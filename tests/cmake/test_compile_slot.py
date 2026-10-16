import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "res/cmake/compile_slot.py"


def launch(slots, lock_dir, code):
    return subprocess.Popen(
        [sys.executable, str(SCRIPT), str(slots), lock_dir, sys.executable, "-c", code],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


class CompileSlotTest(unittest.TestCase):
    def run_second_while_first_holds(self, slots):
        with tempfile.TemporaryDirectory() as lock_dir:
            first = launch(slots, lock_dir, "import time; time.sleep(1.5)")
            try:
                time.sleep(0.3)
                start = time.monotonic()
                second = launch(slots, lock_dir, "pass")
                second.communicate()
                elapsed = time.monotonic() - start
            finally:
                first.communicate()
        return elapsed

    def test_returns_wrapped_exit_code(self):
        with tempfile.TemporaryDirectory() as lock_dir:
            proc = launch(1, lock_dir, "raise SystemExit(3)")
            proc.communicate()
            self.assertEqual(proc.returncode, 3)

    def test_waits_for_the_only_slot(self):
        self.assertGreaterEqual(self.run_second_while_first_holds(1), 1.0)

    def test_two_slots_run_concurrently(self):
        self.assertLess(self.run_second_while_first_holds(2), 1.0)

    def test_rejects_zero_slots(self):
        with tempfile.TemporaryDirectory() as lock_dir:
            proc = launch(0, lock_dir, "pass")
            _, err = proc.communicate()
            self.assertEqual(proc.returncode, 2)
            self.assertIn("slots", err)


if __name__ == "__main__":
    unittest.main()
