#!/usr/bin/env python3

import unittest

from pyphare.simulator.simulator import startMPI
from tests.simulator.test_restart_time import CONSTANT, RecordingTest


class AdvancePastFinalTimeTest(RecordingTest):
    def test_advance_past_final_time_raises(self):
        self.run_recording({}, advance_past_end=True, **CONSTANT)


if __name__ == "__main__":
    startMPI()
    unittest.main()
