#!/usr/bin/env python3

import unittest

import numpy as np

from pyphare.simulator import simulator
from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator.test_adaptive_stable_dt import CFL_WAVE, AdaptiveDtCase, P


class StopOnLevel0DtShrinkTest(AdaptiveDtCase):
    def test_shrinking_level0_dt_stops_the_run(self):
        sim = self.uniform_simulation(
            False,
            (0.3, 0.2),
            smallest_patch_size=5,
            time_step={
                "mode": "adaptive",
                "cfl_wave": CFL_WAVE,
                "level0_dt_shrink_max": 1,
            },
        )

        def pressure(x, y):
            return P + 5.0 * np.exp(-((x - 3.0) ** 2 + (y - 1.0) ** 2) / 0.25)

        self.uniform_model(1.0, pressure)
        self.simulator = Simulator(sim).initialize()
        simulator.exit_on_exception = False
        try:
            with self.assertRaises(RuntimeError) as ctx:
                for _ in range(20):
                    self.simulator.advance()
        finally:
            simulator.exit_on_exception = True
        for text in ("stable dt on level 0", "level0_dt_shrink_max = 1"):
            self.assertIn(text, str(ctx.exception))


if __name__ == "__main__":
    startMPI()
    unittest.main()
