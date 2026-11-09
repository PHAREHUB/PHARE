import unittest

import numpy as np
from ddt import data, ddt, unpack

import pyphare.core.box as boxm
import pyphare.pharein as ph
from pyphare.core.box import Box
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest
from tests.simulator.utilities.field_coarsening import coarsen

ph.NO_GUI()

CELLS = 40
DL = 0.1
SENTINEL = -144123


def blast_pressure(x, y):
    r = np.sqrt((x - 2.0) ** 2 + (y - 2.0) ** 2)
    return 0.1 + 4.95 * (1.0 - np.tanh((r - 0.3) / 0.1))


def box_sim(
    levels,
    diag_dir,
    density,
    pressure,
    time_step,
    final_time,
    time_step_nbr=None,
    hall=False,
    extra_diag_options=None,
    write_timestamps=None,
):
    refinement_boxes = {"L0": {"B0": [[10, 10], [29, 29]]}}
    if levels == 3:
        refinement_boxes["L1"] = {"B0": [[30, 30], [49, 49]]}

    diag_options = {"dir": diag_dir, "mode": "overwrite"}
    diag_options.update(extra_diag_options or {})

    time_args = {"time_step": time_step}
    if isinstance(time_step, dict):
        time_args["final_time"] = final_time
    else:
        time_args["time_step_nbr"] = time_step_nbr

    sim = ph.Simulation(
        **time_args,
        cells=(CELLS, CELLS),
        dl=(DL, DL),
        refinement_boxes=refinement_boxes,
        diag_options={"format": "phareh5", "options": diag_options},
        strict=True,
        nesting_buffer=1,
        hyper_mode="spatial",
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        mhd_timestepper="TVDRK3",
        hall=hall,
        model_options=["MHDModel"],
        max_mhd_level=levels,
    )

    def zero(x, y):
        return 0.0

    def bx(x, y):
        return 0.1

    ph.MHDModel(
        density=density, vx=zero, vy=zero, vz=zero, bx=bx, by=zero, bz=zero, p=pressure
    )

    if write_timestamps is None:
        write_timestamps = np.linspace(0, final_time, 6)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=write_timestamps)
    ph.MHDDiagnostics(quantity="rho", write_timestamps=write_timestamps)
    return sim


def fast_magnetosonic(gamma, rho, b_dir, b2, p):
    c02 = gamma * p / rho
    ca2 = b2 / rho
    cadir2 = b_dir * b_dir / rho
    return np.sqrt(
        (c02 + ca2) * 0.5
        + np.sqrt((c02 + ca2) * (c02 + ca2) - 4.0 * c02 * cadir2) * 0.5
    )


def uniform_stable_dt(cfl_wave, gamma, rho, v, b, p, dl):
    b2 = sum(bi * bi for bi in b)
    inv_dt = sum(
        (abs(v[d]) + fast_magnetosonic(gamma, rho, b[d], b2, p)) / dl[d]
        for d in range(len(dl))
    )
    return cfl_wave / inv_dt


def run_steps(sim, final_time):
    simulator = Simulator(sim)
    simulator.initialize()
    steps = []
    while simulator.currentTime() < final_time - 1e-12:
        t_before = simulator.currentTime()
        simulator.advance()
        steps.append((t_before, simulator.currentTime()))
    simulator.reset()
    return steps


def finest_substeps(diag_dir, steps):
    hier = hierarchy_from(h5_filename=diag_dir + "/mhd_rho.h5")
    fine = [float(t) for t in hier.time_hier if 0 not in hier.levelNbrs(t)]
    margin = 1e-9
    return [
        sum(1 for t in fine if t_start + margin < t < t_end - margin) + 1
        for t_start, t_end in steps
    ]


def coarsened(qty, coarse_box, data_box, coarse_field, fine_field):
    after = np.copy(coarse_field.dataset[:])
    boxm.DataSelector(after)[data_box] = SENTINEL
    coarsen(qty, coarse_field, fine_field, coarse_box, fine_field.dataset[:], after)
    return after


def composite_mass(hier, t):
    level_nbrs = sorted(hier.levelNbrs(t))
    filled = {}
    values = {}
    for ilvl in level_nbrs:
        shape = (CELLS * 2**ilvl,) * 2
        filled[ilvl] = np.zeros(shape, dtype=bool)
        values[ilvl] = np.zeros(shape)
        for patch in hier.level(ilvl, t).patches:
            pd = patch["mhdRho"]
            ghosts = int(pd.layout.nbrGhosts(None))
            lo, up = patch.box.lower, patch.box.upper
            window = (slice(lo[0], up[0] + 1), slice(lo[1], up[1] + 1))
            values[ilvl][window] = pd.dataset[ghosts:-ghosts, ghosts:-ghosts]
            filled[ilvl][window] = True

    mass = 0.0
    for ilvl in level_nbrs:
        weights = np.ones_like(values[ilvl])
        if ilvl + 1 in filled:
            weights[filled[ilvl + 1][0::2, 0::2]] = 0.0
        mass += float(np.sum(values[ilvl] * weights)) * (DL / 2**ilvl) ** 2
    return mass


def covered_increment_errors(diag_dir, t_ref, t, qties):
    hier = hierarchy_from(h5_filename=diag_dir + "/mhd_rho.h5")
    hier = hierarchy_from(h5_filename=diag_dir + "/EM_B.h5", hier=hier)

    errors = {}
    level_nbrs = sorted(set(hier.levelNbrs(t)) & set(hier.levelNbrs(t_ref)))
    for ilvl in level_nbrs[:-1]:
        if ilvl + 1 not in level_nbrs:
            continue
        coarse_patches = hier.level(ilvl, t).patches
        fine_patches = hier.level(ilvl + 1, t).patches
        coarse_patches_ref = hier.level(ilvl, t_ref).patches
        fine_patches_ref = hier.level(ilvl + 1, t_ref).patches

        for qty in qties:
            key = f"{qty}_L{ilvl}"
            errors.setdefault(key, 0.0)
            for ic, coarse_patch in enumerate(coarse_patches):
                for jf, fine_patch in enumerate(fine_patches):
                    overlap = boxm.refine(coarse_patch.box, 2) * fine_patch.box
                    if overlap is None:
                        continue
                    coarse_pd = coarse_patch[qty]
                    fine_pd = fine_patch[qty]
                    coarse_pd_ref = coarse_patches_ref[ic][qty]
                    fine_pd_ref = fine_patches_ref[jf][qty]

                    coarse_box = boxm.coarsen(overlap, 2)
                    coarse_offset = coarse_box.lower - coarse_pd.layout.box.lower
                    data_lower = coarse_offset + coarse_pd.layout.nbrGhosts(None)
                    data_box = Box(data_lower, data_lower + coarse_box.shape - 1)

                    d_coarse = coarse_pd.dataset[:] - coarse_pd_ref.dataset[:]
                    d_fine = coarsened(
                        qty, coarse_box, data_box, coarse_pd, fine_pd
                    ) - coarsened(qty, coarse_box, data_box, coarse_pd_ref, fine_pd_ref)
                    diff = np.abs(
                        boxm.DataSelector(d_coarse)[data_box]
                        - boxm.DataSelector(d_fine)[data_box]
                    )
                    errors[key] = float(np.maximum(errors[key], diff.max()))
    return errors


@ddt
class AdaptiveDtMHDTest(SimulatorTest):
    @data((2, "constant"), (3, "constant"), (2, "adaptive"), (3, "adaptive"))
    @unpack
    def test_covered_coarse_follows_fine(self, levels, mode):
        final_time = 0.3
        timestamps = np.linspace(0, final_time, 6)

        if mode == "adaptive":
            time_args = dict(
                time_step={"mode": "adaptive", "cfl_wave": 0.5},
                final_time=final_time,
            )
        else:
            time_args = dict(time_step=0.002, time_step_nbr=150, final_time=final_time)

        def density(x, y):
            return 1.0 + 0.0 * x

        diag_dir = f"phare_outputs/adaptive_dt_mhd/{mode}_{levels}"
        self.register_diag_dir_for_cleanup(diag_dir)
        ph.global_vars.sim = None
        sim = box_sim(
            levels,
            diag_dir,
            density,
            blast_pressure,
            write_timestamps=timestamps,
            **time_args,
        )

        simulator = Simulator(sim)
        simulator.initialize()
        dts = []
        while simulator.currentTime() < final_time - 1e-12:
            t_before = simulator.currentTime()
            simulator.advance()
            dts.append(simulator.currentTime() - t_before)
        simulator.reset()

        if mode == "adaptive":
            self.assertGreater(max(dts) / min(dts), 2.0)

        hier = hierarchy_from(h5_filename=diag_dir + "/mhd_rho.h5")
        times = sorted(hier.time_hier.keys(), key=float)
        self.assertGreater(len(times), 1)
        for t in times[1:]:
            errors = covered_increment_errors(
                diag_dir, times[0], t, ["mhdRho", "Bx", "By", "Bz"]
            )
            self.assertTrue(errors)
            for key, err in errors.items():
                self.assertLess(err, 1e-12, f"{key} at t={t}")

        mass_hier = hierarchy_from(h5_filename=diag_dir + "/mhd_rho.h5")
        mass_ref = composite_mass(mass_hier, times[0])
        for t in times[1:]:
            drift = composite_mass(mass_hier, t) - mass_ref
            self.assertLess(abs(drift), 1e-12 * mass_ref, f"mass drift at t={t}")

    def test_dt0_is_level0_stable_dt(self):
        cfl_wave = 0.5
        gamma = 5.0 / 3.0

        def density(x, y):
            return 1.0 + 0.0 * x

        def pressure(x, y):
            return 1.0 + 0.0 * x

        diag_dir = "phare_outputs/adaptive_dt_mhd/dt0_level0"
        self.register_diag_dir_for_cleanup(diag_dir)
        ph.global_vars.sim = None
        sim = box_sim(
            2,
            diag_dir,
            density,
            pressure,
            time_step={"mode": "adaptive", "cfl_wave": cfl_wave},
            final_time=1.0,
        )

        expected = uniform_stable_dt(
            cfl_wave, gamma, 1.0, (0.0, 0.0, 0.0), (0.1, 0.0, 0.0), 1.0, (DL, DL)
        )

        simulator = Simulator(sim)
        simulator.initialize()
        dt0 = simulator.timeStep()
        dts = []
        for _ in range(5):
            simulator.advance()
            dts.append(simulator.timeStep())
        simulator.reset()

        self.assertLess(abs(dt0 - expected), 1e-12 * expected, f"{dt0} != {expected}")
        for dt in dts:
            self.assertTrue(np.isfinite(dt))
            self.assertGreater(dt, 0.0)

    def test_finest_level_substeps_follow_stable_dt(self):
        def density(x, y):
            return 1.0 + 0.0 * x

        def pressure(x, y):
            return 1.0 + 0.0 * x

        final_time = 0.1
        diag_dir = "phare_outputs/adaptive_dt_mhd/substeps_adaptive"
        self.register_diag_dir_for_cleanup(diag_dir)
        ph.global_vars.sim = None
        sim = box_sim(
            2,
            diag_dir,
            density,
            pressure,
            time_step={"mode": "adaptive", "cfl_wave": 0.5},
            final_time=final_time,
            extra_diag_options={"fine_dump_lvl_max": 10},
        )
        steps = run_steps(sim, final_time)
        full_steps = steps[:-1]
        adaptive_substeps = finest_substeps(diag_dir, full_steps)

        dt = steps[0][1] - steps[0][0]
        time_step_nbr = 4
        diag_dir = "phare_outputs/adaptive_dt_mhd/substeps_constant"
        self.register_diag_dir_for_cleanup(diag_dir)
        ph.global_vars.sim = None
        sim = box_sim(
            2,
            diag_dir,
            density,
            pressure,
            time_step=dt,
            final_time=dt * time_step_nbr,
            time_step_nbr=time_step_nbr,
            extra_diag_options={"fine_dump_lvl_max": 10},
            write_timestamps=dt * np.arange(time_step_nbr + 1),
        )
        constant_steps = run_steps(sim, dt * time_step_nbr)
        constant_substeps = finest_substeps(diag_dir, constant_steps)

        with self.subTest(mode="adaptive"):
            self.assertGreater(len(adaptive_substeps), 2)
            self.assertEqual(
                adaptive_substeps, [3] + [2] * (len(adaptive_substeps) - 1)
            )
        with self.subTest(mode="constant"):
            self.assertEqual(constant_substeps, [4] * time_step_nbr)

    def test_finest_level_substeps_follow_varying_stable_dt(self):
        def density(x, y):
            return 1.0 + 0.0 * x

        final_time = 0.8
        diag_dir = "phare_outputs/adaptive_dt_mhd/substeps_blast"
        self.register_diag_dir_for_cleanup(diag_dir)
        ph.global_vars.sim = None
        sim = box_sim(
            2,
            diag_dir,
            density,
            blast_pressure,
            time_step={"mode": "adaptive", "cfl_wave": 0.5},
            final_time=final_time,
            extra_diag_options={"fine_dump_lvl_max": 10},
        )
        full_steps = run_steps(sim, final_time)[:-1]
        dts = [t_end - t_start for t_start, t_end in full_steps]
        substeps = finest_substeps(diag_dir, full_steps)

        i_min = int(np.argmin(dts))
        self.assertGreater(max(dts[i_min:]) / dts[i_min], 1.5)
        for i, count in enumerate(substeps):
            self.assertGreaterEqual(count, 2, f"step {i}, dt {dts[i]}")
            self.assertLessEqual(count, 3, f"step {i}, dt {dts[i]}")


if __name__ == "__main__":
    unittest.main()
