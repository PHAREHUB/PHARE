import unittest

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator.advance.adaptive_dt_mhd import (
    CELLS,
    DL,
    AdaptiveDtMHDTestBase,
    box_sim,
    composite_mass,
    finest_substeps,
    run_steps,
    uniform_stable_dt,
)

ph.NO_GUI()


@ddt
class AdaptiveDtMHDTest(AdaptiveDtMHDTestBase):
    @data((2, "adaptive", False))
    @unpack
    def test_covered_coarse_follows_fine(self, levels, mode, hall):
        self.check_covered_coarse_follows_fine(levels, mode, hall)

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

    def test_level_created_by_regridding(self):
        lx = CELLS * DL

        def density(x, y):
            return 1.0 + 0.0 * x

        def pressure(x, y):
            return 0.1 + 0.0 * x

        def vx(x, y):
            return np.sin(2.0 * np.pi * x / lx)

        def by(x, y):
            return 0.5 + 0.0 * x

        final_time = 0.45
        diag_dir = "phare_outputs/adaptive_dt_mhd/regrid_created_level"
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
            vx=vx,
            by=by,
            tagging_threshold=0.1,
        )
        full_steps = run_steps(sim, final_time)[:-1]
        substeps = finest_substeps(diag_dir, full_steps)

        hier = hierarchy_from(h5_filename=diag_dir + "/mhd_rho.h5")
        times = sorted(hier.time_hier.keys(), key=float)
        self.assertEqual(sorted(hier.levelNbrs(times[0])), [0])

        created = next(i for i, count in enumerate(substeps) if count > 1)
        self.assertGreater(created, 0)
        self.assertEqual(substeps[created], 3)
        t_start, t_end = full_steps[created]
        first_fine_dump = min(
            float(t)
            for t in times
            if 0 not in hier.levelNbrs(t) and t_start < float(t) < t_end
        )
        self.assertAlmostEqual(
            first_fine_dump - t_start, (t_end - t_start) / 4.0, delta=1e-9
        )
        later = substeps[created + 1 :]
        self.assertGreater(len(later), 2)
        for i, count in enumerate(later):
            self.assertIn(count, (2, 3), f"step {created + 1 + i}")

        coarse_times = [t for t in times if 0 in hier.levelNbrs(t)]
        self.assertEqual(
            {tuple(sorted(hier.levelNbrs(t))) for t in coarse_times}, {(0,), (0, 1)}
        )
        rho_dtype = hier.level(0, times[0]).patches[0]["mhdRho"].dataset[:].dtype
        mass_rtol = 1e-12 if np.dtype(rho_dtype) == np.float64 else 1e-6
        mass_ref = composite_mass(hier, coarse_times[0])
        for t in coarse_times[1:]:
            drift = composite_mass(hier, t) - mass_ref
            self.assertLess(abs(drift), mass_rtol * mass_ref, f"mass drift at t={t}")


if __name__ == "__main__":
    unittest.main()
