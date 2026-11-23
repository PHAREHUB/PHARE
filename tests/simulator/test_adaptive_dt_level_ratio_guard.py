#!/usr/bin/env python3

import unittest

from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator.test_adaptive_stable_dt import CFL_WAVE, AdaptiveDtCase


class StopOnFineLevelSubstepsTest(AdaptiveDtCase):
    def test_too_many_fine_level_substeps_stops_the_run(self):
        sim = self.uniform_simulation(
            False,
            (0.3, 0.2),
            max_mhd_level=2,
            smallest_patch_size=5,
            refinement_boxes={"L0": {"B0": [[4, 2], [15, 7]]}},
            time_step={
                "mode": "adaptive",
                "cfl_wave": CFL_WAVE,
                "fine_level_substeps_max": 1,
            },
        )
        self.uniform_model(1.0)
        self.simulator = Simulator(sim)
        with self.assertRaises(ValueError) as ctx:
            self.simulator.initialize()
        cause = ctx.exception.__cause__
        self.assertIsInstance(cause, RuntimeError)
        for text in ("stable dt on level 1", "level 0", "fine_level_substeps_max = 1"):
            self.assertIn(text, str(cause))


if __name__ == "__main__":
    startMPI()
    unittest.main()
