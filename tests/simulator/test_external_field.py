#!/usr/bin/env python3
"""
Functional tests for the external magnetic field B0.
"""

import os

import numpy as np

import pyphare.pharein as ph
from pyphare.core.phare_utilities import any_fp_tol
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = (40, 40)
dl = (0.2, 0.2)
time_step = 0.001


def az(x, y):
    """Bilinear vector potential: curl is (y, -x),
    and is exactly computed by the discrete curl."""
    return x * y


def az_twice(x, y):
    return 2.0 * x * y


def az_t(x, y, t):
    return x * y * (1.0 + t)


def dazdt(x, y, t):
    return x * y


def expected_b0(name, x, y, t=0.0):
    """B0 = curl(A0) for the potentials above, evaluated at the node coordinates."""
    return {
        "B0x": x * (1.0 + t),
        "B0y": -y * (1.0 + t),
        "B0z": np.zeros_like(x),
    }[name]


def plasma_model(sim):
    """A pair of current sheets, deliberately not in pressure balance so that B
    evolves fast enough to move the tagged region and force real regrids."""

    def sheets(y):
        Ly = sim.simulation_domain()[1]
        return np.tanh((y - 0.3 * Ly) / 0.5) - np.tanh((y - 0.7 * Ly) / 0.5)

    def density(x, y):
        return 0.4 + 1.0 / np.cosh(sheets(y)) ** 2

    def bx(x, y):
        return -1.0 + 2.0 * sheets(y)

    def zero(x, y):
        return 0.0

    def vth(x, y):
        return 0.3

    ph.MaxwellianFluidModel(
        bx=bx,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": density,
            **{f"vbulk{c}": zero for c in "xyz"},
            **{f"vth{c}": vth for c in "xyz"},
            "nbr_part_per_cell": 25,
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)


class ExternalFieldTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super().tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    # -- helpers ---------------------------------------------------------------

    def config(self, nbr_steps, external_field=None, **kwargs):
        """A short 2D run, dumping EM_B0 at every step plus the initial time when an
        external field is declared."""
        timestamps = np.arange(nbr_steps + 1) * time_step
        opts = dict(
            time_step=time_step,
            time_step_nbr=nbr_steps,
            cells=cells,
            dl=dl,
            diag_options={
                "format": "phareh5",
                "options": {
                    "dir": "phare_outputs/test_external_field",
                    "mode": "overwrite",  # a failed run leaves its diags behind
                },
            },
        )
        opts.update(kwargs)
        if external_field is not None:
            opts["external_field"] = external_field

        sim = self.simulation(**opts)
        plasma_model(sim)
        if external_field is not None:
            ph.ElectromagDiagnostics(quantity="B0", write_timestamps=timestamps)
        return sim, timestamps

    def run_sim(self, sim):
        # NB: not named `run` - unittest.TestCase.run is the test entry point
        Simulator(sim).run().reset()
        return Run(sim.diag_options["options"]["dir"])

    def b0_levels(self, run, time):
        hier = run.GetB0(time, all_primal=False)
        return hier.levels(time)

    def assert_b0_is_exact(self, run, time, t_expected=0.0, min_levels=1):
        """Every B0 node equals curl(A0) evaluated at its own coordinates."""
        levels = self.b0_levels(run, time)
        self.assertGreaterEqual(
            len(levels), min_levels, f"t={time}: only {len(levels)} level(s)"
        )
        nodes = 0
        for ilvl, level in levels.items():
            self.assertTrue(level.patches, f"t={time} level {ilvl} has no patch")
            for patch in level.patches:
                for name, pd in patch.patch_datas.items():
                    x, y = pd.meshgrid()
                    data = pd.dataset[:]
                    self.assertTrue(
                        np.all(np.isfinite(data)),
                        f"t={time} lvl={ilvl} {name}: non finite B0",
                    )
                    expected = expected_b0(name, x, y, t_expected)
                    np.testing.assert_allclose(
                        data,
                        expected,
                        err_msg=f"t={time} lvl={ilvl} {name}",
                        **any_fp_tol(data, expected, atol=1e-12, atol_fp32=1e-6),
                    )
                    nodes += data.size
        self.assertGreater(nodes, 0, f"t={time}: nothing was checked")

    def test_b0_is_exact_on_every_level_at_init(self):
        """Covers the initLevel path: a refined level is never passed to
        model.initialize, so its B0 only exists if the level initializer fills it."""
        sim, timestamps = self.config(
            nbr_steps=1,
            external_field={"type": "user-defined", "potential": (None, None, az)},
            refinement_boxes={"L0": {"B0": [(10, 10), (29, 29)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
        )
        run = self.run_sim(sim)
        # min_levels=2: without it a hierarchy that failed to refine passes vacuously
        self.assert_b0_is_exact(run, timestamps[0], min_levels=2)

        # the primal-interpolated read path, which needs B0x/B0y/B0z in yee_centering
        b0 = run.GetB0(timestamps[0])
        self.assertTrue(b0.levels()[1].patches[0].attrs)

    def test_b0_is_exact_after_regrid(self):
        """Covers the regrid path: patches born from a regrid get their state from
        the messenger, which knows nothing about B0."""
        sim, timestamps = self.config(
            nbr_steps=10,
            external_field={"type": "user-defined", "potential": (None, None, az)},
            refinement="tagging",
            max_nbr_levels=2,
            nesting_buffer=1,
            clustering="tile",
            tag_buffer="1",
        )
        run = self.run_sim(sim)

        for time in timestamps:
            self.assert_b0_is_exact(run, time, min_levels=2)

    def test_time_dependent_b0_follows_the_dump_time(self):
        """Covers the end-of-step update: B0 dumped at t must be B0(t), not B0(0)."""
        sim, timestamps = self.config(
            nbr_steps=5,
            external_field={
                "type": "user-defined",
                "potential": (None, None, az_t),
                "potential_time_derivative": (None, None, dazdt),
            },
            refinement_boxes={"L0": {"B0": [(10, 10), (29, 29)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
        )
        run = self.run_sim(sim)

        for time in timestamps:
            self.assert_b0_is_exact(run, time, t_expected=time, min_levels=2)

        # and it really did change: a frozen B0 would satisfy t=0 only
        first = run.GetB0(timestamps[0], all_primal=False)
        last = run.GetB0(timestamps[-1], all_primal=False)
        self.assertFalse(
            np.array_equal(
                first.levels(timestamps[0])[0].patches[0].patch_datas["B0x"].dataset[:],
                last.levels(timestamps[-1])[0].patches[0].patch_datas["B0x"].dataset[:],
            ),
            "time dependent B0 did not change over the run",
        )

    def test_static_b0_is_untouched_by_advance(self):
        """A time independent field must be computed once and then left alone."""
        sim, timestamps = self.config(
            nbr_steps=5,
            external_field={
                "type": "dipole",
                "position": (-5.0, -5.0),  # outside the domain: no singularity
                "moment": (0.0, 1.0),
                "radius": 0.0,
            },
            refinement_boxes={"L0": {"B0": [(10, 10), (29, 29)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
        )
        run = self.run_sim(sim)

        def b0_of(time):
            out = {}
            for ilvl, level in self.b0_levels(run, time).items():
                for patch in level.patches:
                    for name, pd in patch.patch_datas.items():
                        out[(ilvl, str(patch.box), name)] = pd.dataset[:].copy()
            return out

        first, last = b0_of(timestamps[0]), b0_of(timestamps[-1])
        self.assertTrue(first, "no B0 data was dumped")
        for key, values in first.items():
            self.assertTrue(np.all(np.isfinite(values)), f"{key}: non finite B0")
            self.assertIn(key, last)
            np.testing.assert_array_equal(values, last[key], err_msg=str(key))

    def test_dipole_radius_regularizes_the_center(self):
        """A dipole of finite radius sitting on a node inside the domain: B0 must be
        finite everywhere (a point dipole divides by zero there), and inside the radius
        the uniform field m / (2 pi R^2) of a uniformly magnetized cylinder."""
        center = tuple(0.5 * n * d for n, d in zip(cells, dl))  # a node on every level
        moment = (0.0, 1.0)
        radius = 1.0
        sim, timestamps = self.config(
            nbr_steps=1,
            external_field={
                "type": "dipole",
                "position": center,
                "moment": moment,
                "radius": radius,
            },
            refinement_boxes={"L0": {"B0": [(10, 10), (29, 29)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
        )
        run = self.run_sim(sim)

        uniform = {
            "B0x": moment[0] / (2.0 * np.pi * radius**2),
            "B0y": moment[1] / (2.0 * np.pi * radius**2),
            "B0z": 0.0,
        }
        levels = self.b0_levels(run, timestamps[0])
        self.assertGreaterEqual(len(levels), 2)
        inside = 0
        for ilvl, level in levels.items():
            for patch in level.patches:
                for name, pd in patch.patch_datas.items():
                    data = pd.dataset[:]
                    self.assertTrue(
                        np.all(np.isfinite(data)), f"lvl={ilvl} {name}: non finite B0"
                    )
                    x, y = pd.meshgrid()
                    r = np.hypot(x - center[0], y - center[1])
                    # the curl stencil reaches half a cell away: keep it all inside
                    mask = r + pd.layout.dl[0] < radius
                    expected = np.full_like(data[mask], uniform[name])
                    np.testing.assert_allclose(
                        data[mask],
                        expected,
                        err_msg=f"lvl={ilvl} {name}",
                        **any_fp_tol(data, expected, atol=1e-10, atol_fp32=1e-6),
                    )
                    inside += np.count_nonzero(mask)
        self.assertGreater(inside, 0, "no node was checked inside the radius")

    def restart_pair(self, first_field, restarted_field):
        restart_time = 2 * time_step
        timestamps = [restart_time, 4 * time_step]
        out = "phare_outputs/test_external_field_restart"

        def options(external_field, diag_dir, restart_options):
            opts = dict(
                time_step=time_step,
                time_step_nbr=4,
                cells=cells,
                dl=dl,
                refinement_boxes={"L0": {"B0": [(10, 10), (29, 29)]}},
                smallest_patch_size=10,
                largest_patch_size=20,
                diag_options={
                    "format": "phareh5",
                    "options": {"dir": diag_dir, "mode": "overwrite"},
                },
                restart_options=restart_options,
            )
            if external_field is not None:
                opts["external_field"] = external_field
            return opts

        def setup(sim, external_field):
            plasma_model(sim)
            ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)
            if external_field is not None:
                ph.ElectromagDiagnostics(quantity="B0", write_timestamps=timestamps)

        restart_options = {
            "dir": out,
            "mode": "overwrite",
            "timestamps": [restart_time],
        }
        sim = self.simulation(**options(first_field, out, restart_options))
        setup(sim, first_field)
        Simulator(sim).run().reset()
        first_dir = sim.diag_options["options"]["dir"]
        restarted_dir = first_dir + "_restarted"

        restart_options = {
            "dir": first_dir,
            "mode": "overwrite",
            "restart_time": restart_time,
        }
        ph.global_vars.sim = None
        sim = ph.Simulation(**options(restarted_field, restarted_dir, restart_options))
        setup(sim, restarted_field)
        self.register_diag_dir_for_cleanup(restarted_dir)
        Simulator(sim).run().reset()
        return Run(restarted_dir), timestamps

    def test_restart_declaring_a_field_computes_b0(self):
        run, timestamps = self.restart_pair(
            None, {"type": "user-defined", "potential": (None, None, az)}
        )
        for time in timestamps:
            self.assert_b0_is_exact(run, time, min_levels=2)

    def test_restart_with_other_parameters_takes_the_new_field(self):
        run, timestamps = self.restart_pair(
            {"type": "user-defined", "potential": (None, None, az_twice)},
            {"type": "user-defined", "potential": (None, None, az)},
        )
        for time in timestamps:
            self.assert_b0_is_exact(run, time, min_levels=2)

    def test_restart_dropping_the_field_runs(self):
        run, timestamps = self.restart_pair(
            {"type": "user-defined", "potential": (None, None, az)}, None
        )
        self.assertFalse(os.path.exists(os.path.join(run.path, "EM_B0.h5")))
        for time in timestamps:
            b = run.GetB(time, all_primal=False)
            self.assertGreaterEqual(len(b.levels(time)), 2)

    def restart_entries(self, external_field):
        import glob

        import h5py

        opts = dict(
            time_step=time_step,
            time_step_nbr=2,
            cells=cells,
            dl=dl,
            refinement_boxes={"L0": {"B0": [(10, 10), (29, 29)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
            diag_options={
                "format": "phareh5",
                "options": {
                    "dir": "phare_outputs/test_external_field_entries",
                    "mode": "overwrite",
                },
            },
            restart_options={
                "dir": "phare_outputs/test_external_field_entries",
                "mode": "overwrite",
                "timestamps": [time_step],
            },
        )
        if external_field is not None:
            opts["external_field"] = external_field
        sim = self.simulation(**opts)
        plasma_model(sim)
        Simulator(sim).run().reset()

        entries = set()
        pattern = os.path.join(sim.restart_options["dir"], "**", "proc.*")
        for path in glob.glob(pattern, recursive=True):
            with h5py.File(path, "r") as restart_file:
                restart_file.visit(lambda name: entries.add(name.split("/")[-1]))
        return entries

    def test_restart_files_hold_b0_only_when_declared(self):
        declared = self.restart_entries(
            {"type": "user-defined", "potential": (None, None, az)}
        )
        self.assertIn("HybridModel_B0##default", declared)
        self.assertIn("HybridModel_dB0dt##default", declared)

        undeclared = self.restart_entries(None)
        self.assertTrue(undeclared)
        self.assertFalse([name for name in undeclared if "B0" in name])

    def test_b0_readback_without_field_raises(self):
        sim, _ = self.config(
            nbr_steps=2,
            refinement_boxes={"L0": {"B0": [(10, 10), (29, 29)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
        )
        self.simulator = Simulator(sim)
        self.simulator.initialize()
        self.simulator.advance()
        self.simulator.advance()

        with self.assertRaises(RuntimeError):
            hierarchy_from(simulator=self.simulator, qty="EM_B0_x")


if __name__ == "__main__":
    import unittest

    unittest.main()
