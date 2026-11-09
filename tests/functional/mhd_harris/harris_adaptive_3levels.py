#!/usr/bin/env python3

from common import HarrisTest, full_cells

from pyphare.simulator.simulator import startMPI


class HarrisAdaptive3LevelsTest(HarrisTest):
    def test_run(self):
        return self.run_case(
            full_cells,
            final_time=50,
            diag_dir="phare_outputs/mhd_harris_adaptive_3levels",
            time_step={"mode": "adaptive", "cfl_wave": 0.8},
            levels=3,
        )


if __name__ == "__main__":
    startMPI()
    HarrisAdaptive3LevelsTest().test_run().tearDown()
