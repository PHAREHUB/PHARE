#!/usr/bin/env python3

import unittest
from unittest import mock

from pyphare.simulator.simulator import startMPI
from tests.simulator.test_restart_time import CONSTANT, RecordingTest


class StartTimeMismatchTest(RecordingTest):
    def test_python_start_time_must_match_the_restored_time(self):
        self.run_recording(dict(timestamps=[0.3]), **CONSTANT)
        with (
            mock.patch(
                "pyphare.pharein.restarts.restored_time", lambda options: 0.3 + 1e-9
            ),
            self.assertRaises((RuntimeError, ValueError)) as cm,
        ):
            self.run_recording(dict(restart_time=0.3), time_step=0.1, time_step_nbr=3)
        message = str(cm.exception) + str(cm.exception.__cause__)
        self.assertIn("differs from the restored time", message)


if __name__ == "__main__":
    startMPI()
    unittest.main()
