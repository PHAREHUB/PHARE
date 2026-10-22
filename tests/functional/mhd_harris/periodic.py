#!/usr/bin/env python3

from common import HarrisTest, full_cells

from pyphare.simulator.simulator import startMPI


class HarrisPeriodicTest(HarrisTest):
    def test_run(self):
        return self.run_case(
            full_cells, final_time=50, diag_dir="phare_outputs/mhd_harris"
        )


if __name__ == "__main__":
    startMPI()
    HarrisPeriodicTest().test_run().tearDown()
