#!/usr/bin/env python3

import unittest

import numpy as np

from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator.test_adaptive_stable_dt import AdaptiveDtCase, P


class StopOnCorruptStateTest(AdaptiveDtCase):
    def test_corrupt_state_stops_the_run(self):
        sim = self.uniform_simulation(False, (0.3, 0.2), smallest_patch_size=5)

        def pressure(x, y):
            return np.where((x < 0.6) & (y < 0.4), -1.0, P) + 0.0 * x

        self.uniform_model(1.0, pressure)
        self.simulator = Simulator(sim)
        with self.assertRaises(ValueError) as ctx:
            self.simulator.initialize()
        cause = ctx.exception.__cause__
        self.assertIsInstance(cause, RuntimeError)
        self.assertIn("non-finite wave speed on level 0", str(cause))


if __name__ == "__main__":
    startMPI()
    unittest.main()
