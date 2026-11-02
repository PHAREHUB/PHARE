#!/usr/bin/env python3

from common import ShockTest, dx, ncells

from pyphare.simulator.simulator import startMPI

final_time = 80


class ShockPeriodicTest(ShockTest):
    def test_run(self):
        return self.run_case(
            time_step=0.2,
            final_time=final_time,
            timestamps=[final_time],
            diag_dir="phare_outputs/shock",
            x0=ncells * dx / 2,
            u=0.0,
        )


if __name__ == "__main__":
    startMPI()
    ShockPeriodicTest().test_run().tearDown()
