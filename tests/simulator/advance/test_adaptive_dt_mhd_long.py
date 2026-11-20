import unittest

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from tests.simulator.advance.adaptive_dt_mhd import (
    AdaptiveDtMHDTestBase,
    blast_pressure,
    box_sim,
    finest_substeps,
    run_steps,
)

ph.NO_GUI()


@ddt
class AdaptiveDtMHDLongTest(AdaptiveDtMHDTestBase):
    @data(
        (2, "constant", False),
        (3, "constant", False),
        (3, "adaptive", False),
        (2, "adaptive", True),
        (3, "adaptive", True),
    )
    @unpack
    def test_covered_coarse_follows_fine(self, levels, mode, hall):
        self.check_covered_coarse_follows_fine(levels, mode, hall)

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
