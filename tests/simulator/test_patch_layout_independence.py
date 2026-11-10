#!/usr/bin/env python3
"""
With a seed, the initial particles of a cell do not depend on the patches, so a run must not
depend on the patch layout: the same deck run with small patches and with one patch per level
must agree on every level, to round-off.

- a uniform plasma, two levels (one refinement box), 8-cell patches against one patch per
  level, 1D and 2D, every interpolation order: same number of domain particles per level at
  every step, and same moments, E and B to round-off. Level ghost particles are held by every
  patch whose ghost box covers them: one entering the domain of a same-level patch belongs to
  that patch, and each is deposited by one patch only.
- one level of patches G + 1 and G + 2 cells wide, G the field ghost width: two patches are
  then G + 1 cells apart, which a border sum schedule does not pair, so the outermost ghost
  node of one of them misses the other's deposit (#1575). That node must not reach the domain.
- on every run, the nodes shared by the domains of several patches (patch borders, periodic
  images included) hold the same bits on all of them; ghost nodes may differ.
- a run with tagging and regrids, overlapping patches included: shared domain nodes as above,
  and the number of level 0 particles is constant (periodic domain).
"""

import unittest

import h5py
import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core.box import Box1D, Box2D
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

dl = 0.125
time_step = 0.00125
time_step_nbr = 5
small_patch_size = 8
field_ghost_width = {1: 2, 2: 4, 3: 4}  # hybrid_options.hpp

two_level_cells = {1: 64, 2: 32}
refinement_box = {1: (16, 47), 2: (8, 23)}  # same bounds along every direction
ppc = {1: 128, 2: 64}

cases = [(ndim, interp) for ndim in (1, 2) for interp in (1, 2, 3)]
# interp 1 patches cannot be 4 cells wide (pharein, SAMRAI issue 311): G + 1 only
gap_cases = [
    (ndim, interp, extra)
    for ndim, interp in cases
    for extra in (1, 2)
    if not (interp == 1 and extra == 2)
]

# file, dataset keys
FIELDS = {
    "n": ("ions_charge_density.h5", ["charge_density"]),
    "flux": ("ions_charge_flux.h5", [f"charge_flux_{c}" for c in "xyz"]),
    "E": ("EM_E.h5", [f"EM_E_{c}" for c in "xyz"]),
    "B": ("EM_B.h5", [f"EM_B_{c}" for c in "xyz"]),
}


def const(value):
    return lambda *xyz: value + np.zeros_like(xyz[0])


def timestamps(nbr_steps=time_step_nbr, dt=time_step):
    return np.arange(nbr_steps + 1) * dt


def add_diagnostics(times):
    for quantity in ("E", "B"):
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=times)
    for quantity in ("charge_density", "charge_flux"):
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=times)
    ph.ParticleDiagnostics(quantity="domain", write_timestamps=times, population_name="protons")


def uniform_plasma(ndim, nbr_part_per_cell):
    ph.MaxwellianFluidModel(
        bx=const(0.0),
        by=const(0.0),
        bz=const(1.0),
        protons={
            "charge": 1,
            "mass": 1,
            "density": const(1.0),
            "vbulkx": const(0.0),
            "vbulky": const(0.0),
            "vbulkz": const(0.0),
            "vthx": const(0.3),
            "vthy": const(0.3),
            "vthz": const(0.3),
            "nbr_part_per_cell": nbr_part_per_cell,
            "init": {"seed": 11},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)


def config(test, ndim, interp, cells, patch_size, out, refinement_boxes=None):
    sim = test.simulation(
        cells=[cells] * ndim,
        dl=[dl] * ndim,
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        interp_order=interp,
        smallest_patch_size=patch_size,
        largest_patch_size=patch_size,
        refinement_boxes=refinement_boxes,
        resistivity=0.0,
        hyper_resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": f"phare_outputs/{out}", "mode": "overwrite"},
        },
    )
    uniform_plasma(ndim, ppc[ndim])
    add_diagnostics(timestamps())
    return sim


def time_key(h5, time):
    return next(k for k in h5["t"] if abs(float(k) - time) < 1e-10)


def patches_of(h5, time, ilvl):
    level = h5["t"][time_key(h5, time)].get(f"pl{ilvl}")
    return [] if level is None else list(level.values())


def domain_nodes(dataset, patch):
    """the AMR index of the first domain node and the domain node values of `dataset`"""
    ghosts = int(dataset.attrs["ghosts"])
    lower = np.asarray(patch.attrs["lower"], dtype=int)
    ncells = np.asarray(patch.attrs["nbrCells"], dtype=int)
    nodes = [n - 2 * ghosts for n in dataset.shape]  # nbr of domain nodes, + 1 if primal
    assert all(n in (c, c + 1) for n, c in zip(nodes, ncells)), (dataset.name, nodes, ncells)
    return lower, dataset[...][tuple(slice(ghosts, ghosts + n) for n in nodes)]


class LevelField:
    """the domain nodes of one dataset on one level, from all patches, each node once; checks
    that every node shared by several patches has the same bits on each of them"""

    def __init__(self, h5, time, ilvl, key, period):
        self.mismatches = []
        values = {}
        for patch in patches_of(h5, time, ilvl):
            lower, nodes = domain_nodes(patch[key], patch)
            for index in np.ndindex(nodes.shape):
                node = tuple((lower[d] + index[d]) % period[d] for d in range(len(index)))
                value = nodes[index]
                if node in values:
                    if values[node].tobytes() != value.tobytes():
                        self.mismatches.append((node, values[node], value))
                else:
                    values[node] = value
        self.values = values


def tolerance(dtype):
    return 1e-12 if dtype == np.float64 else 1e-5


@ddt
class PatchLayoutIndependenceTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_sim(self, sim):
        Simulator(sim).run().reset()
        return sim.diag_options["options"]["dir"]

    def check_shared_domain_nodes(self, diag_dir, ndim, cells, times):
        with h5py.File(f"{diag_dir}/{FIELDS['B'][0]}") as h5:
            levels = range(max(len(h5["t"][time_key(h5, t)]) for t in times))
        for name, (filename, keys) in FIELDS.items():
            with h5py.File(f"{diag_dir}/{filename}") as h5:
                for time in times:
                    for ilvl in levels:
                        period = [cells * 2**ilvl] * ndim
                        for key in keys:
                            field = LevelField(h5, time, ilvl, key, period)
                            self.assertEqual(
                                field.mismatches[:3],
                                [],
                                f"{key} t={time} level {ilvl}: {len(field.mismatches)} nodes",
                            )

    def compare_layouts(self, dir_a, dir_b, ndim, cells, times):
        for name, (filename, keys) in FIELDS.items():
            with h5py.File(f"{dir_a}/{filename}") as ha, h5py.File(f"{dir_b}/{filename}") as hb:
                for time in times:
                    for ilvl in range(len(ha["t"][time_key(ha, time)])):
                        period = [cells * 2**ilvl] * ndim
                        for key in keys:
                            a = LevelField(ha, time, ilvl, key, period).values
                            b = LevelField(hb, time, ilvl, key, period).values
                            self.assertEqual(a.keys(), b.keys(), f"{key} level {ilvl}")
                            va = np.array([a[n] for n in sorted(a)])
                            vb = np.array([b[n] for n in sorted(a)])
                            atol = tolerance(va.dtype) * max(1.0, np.abs(va).max())
                            np.testing.assert_allclose(
                                va, vb, rtol=0, atol=atol, err_msg=f"{key} t={time} level {ilvl}"
                            )

    def particle_counts(self, diag_dir, times):
        """number of domain particles per (time, level)"""
        counts = {}
        with h5py.File(f"{diag_dir}/ions_pop_protons_domain.h5") as h5:
            for time in times:
                for ilvl in range(len(h5["t"][time_key(h5, time)])):
                    counts[(time, ilvl)] = sum(
                        patch["weight"].shape[0] for patch in patches_of(h5, time, ilvl)
                    )
        return counts

    def nbr_patches(self, diag_dir, ilvl):
        with h5py.File(f"{diag_dir}/{FIELDS['B'][0]}") as h5:
            return len(patches_of(h5, 0.0, ilvl))

    @data(*cases)
    @unpack
    def test_two_levels(self, ndim, interp):
        print(f"{self._testMethodName}: {ndim}D interp {interp}")
        cells = two_level_cells[ndim]
        lower, upper = refinement_box[ndim]
        boxes = {"L0": {"B0": (Box1D if ndim == 1 else Box2D)(lower, upper)}}
        small = self.run_sim(
            config(self, ndim, interp, cells, small_patch_size, "patch_layout/small", boxes)
        )
        whole = self.run_sim(config(self, ndim, interp, cells, cells, "patch_layout/whole", boxes))

        self.assertEqual(self.nbr_patches(whole, 1), 1)
        self.assertGreater(self.nbr_patches(small, 1), 1)

        small_counts = self.particle_counts(small, timestamps())
        whole_counts = self.particle_counts(whole, timestamps())
        self.assertEqual(small_counts.keys(), whole_counts.keys())
        self.assertIn((timestamps()[-1], 1), small_counts)
        for key, count in whole_counts.items():
            self.assertEqual(small_counts[key], count, f"(time, level) {key}")

        self.check_shared_domain_nodes(small, ndim, cells, timestamps())
        self.compare_layouts(small, whole, ndim, cells, timestamps())

    @data(*gap_cases)
    @unpack
    def test_patches_ghost_width_plus(self, ndim, interp, extra):
        print(f"{self._testMethodName}: {ndim}D interp {interp}, patches G + {extra} wide")
        width = field_ghost_width[interp] + extra
        cells = width * 6
        small = self.run_sim(config(self, ndim, interp, cells, width, "patch_layout/small"))
        whole = self.run_sim(config(self, ndim, interp, cells, cells, "patch_layout/whole"))

        self.assertEqual(self.nbr_patches(whole, 0), 1)
        self.assertEqual(self.nbr_patches(small, 0), 6**ndim)

        self.check_shared_domain_nodes(small, ndim, cells, timestamps())
        self.compare_layouts(small, whole, ndim, cells, timestamps())

    def test_regrids_with_tagging(self):
        ndim, interp, cells, steps = 2, 2, 32, 8
        times = timestamps(steps)
        sim = self.simulation(
            cells=[cells] * ndim,
            dl=[0.2] * ndim,
            time_step=time_step,
            time_step_nbr=steps,
            interp_order=interp,
            refinement="tagging",
            max_nbr_levels=3,
            nesting_buffer=1,
            tagging_threshold=0.1,
            resistivity=0.0,
            # without: E and B copies differ in the last bits where the current is not zero,
            # from the hyper-resistive term, a separate cause
            hyper_resistivity=0.0,
            diag_options={
                "format": "phareh5",
                "options": {"dir": "phare_outputs/patch_layout/tagging", "mode": "overwrite"},
            },
        )

        def bx(x, y):  # a bump, tagged, carried by the flow so that regrids move the levels
            return 1.0 + np.exp(-((x - 3.2) ** 2 + (y - 3.2) ** 2) / 0.8**2)

        ph.MaxwellianFluidModel(
            bx=bx,
            by=const(0.0),
            bz=const(0.0),
            protons={
                "charge": 1,
                "density": const(1.0),
                "vbulkx": const(2.0),
                "vbulky": const(1.0),
                "vbulkz": const(0.0),
                "vthx": const(0.3),
                "vthy": const(0.3),
                "vthz": const(0.3),
                "nbr_part_per_cell": 20,
                "init": {"seed": 1},
            },
        )
        ph.ElectronModel(closure="isothermal", Te=0.0)
        add_diagnostics(times)
        diag_dir = self.run_sim(sim)

        counts = self.particle_counts(diag_dir, times)
        self.assertGreater(max(ilvl for _, ilvl in counts), 0, "no refined level")
        level0 = {counts[(time, 0)] for time in times}
        self.assertEqual(len(level0), 1, f"level 0 particles {level0}")
        self.check_shared_domain_nodes(diag_dir, ndim, cells, times)


if __name__ == "__main__":
    unittest.main()
