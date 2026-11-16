import contextlib
import errno
import importlib.util
import io
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[2] / "res/cmake/compile_slot.py"
TIMEOUT = 120


def load_script():
    spec = importlib.util.spec_from_file_location("compile_slot", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def launch(slots, lock_dir, code):
    return subprocess.Popen(
        [sys.executable, str(SCRIPT), str(slots), lock_dir, sys.executable, "-c", code],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def hold_until(ready, release):
    return (
        "import os, time\n"
        f"open({str(ready)!r}, 'w').close()\n"
        f"deadline = time.monotonic() + {TIMEOUT}\n"
        f"while not os.path.exists({str(release)!r}) and time.monotonic() < deadline:\n"
        "    time.sleep(0.01)\n"
    )


def record_release(ran, release):
    return (
        "import os\n"
        f"open({str(ran)!r}, 'w').write(str(os.path.exists({str(release)!r})))\n"
    )


def wait_for(path):
    deadline = time.monotonic() + TIMEOUT
    while not path.exists():
        if time.monotonic() > deadline:
            raise AssertionError(f"timed out waiting for {path}")
        time.sleep(0.01)


class CompileSlotTest(unittest.TestCase):
    def run_second_while_first_holds(self, slots, settle):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            lock_dir = str(tmp / "locks")
            ready, release, ran = tmp / "ready", tmp / "release", tmp / "ran"
            first = launch(slots, lock_dir, hold_until(ready, release))
            second = None
            try:
                wait_for(ready)
                second = launch(slots, lock_dir, record_release(ran, release))
                settle(second)
            finally:
                release.touch()
                first.communicate()
                if second is not None:
                    try:
                        second.communicate(timeout=TIMEOUT)
                    except subprocess.TimeoutExpired:
                        second.kill()
                        second.communicate()
            self.assertEqual(first.returncode, 0)
            self.assertEqual(second.returncode, 0)
            return ran.read_text()

    def test_returns_wrapped_exit_code(self):
        with tempfile.TemporaryDirectory() as lock_dir:
            proc = launch(1, lock_dir, "raise SystemExit(3)")
            proc.communicate()
            self.assertEqual(proc.returncode, 3)

    def test_waits_for_the_only_slot(self):
        def still_waiting(second):
            time.sleep(0.5)
            self.assertIsNone(second.poll())

        self.assertEqual(self.run_second_while_first_holds(1, still_waiting), "True")

    def test_two_slots_run_concurrently(self):
        def finishes(second):
            second.communicate(timeout=TIMEOUT)

        self.assertEqual(self.run_second_while_first_holds(2, finishes), "False")

    def test_rejects_zero_slots(self):
        with tempfile.TemporaryDirectory() as lock_dir:
            proc = launch(0, lock_dir, "pass")
            _, err = proc.communicate()
            self.assertEqual(proc.returncode, 2)
            self.assertIn("slots", err)


class CompileSlotLockingTest(unittest.TestCase):
    def setUp(self):
        self.script = load_script()

    def test_busy_slot_is_retried(self):
        with tempfile.TemporaryDirectory() as lock_dir:
            with patch.object(
                self.script.fcntl, "flock", side_effect=[BlockingIOError(), None]
            ) as flock:
                fd = self.script.acquire(1, lock_dir)
            self.assertIsNotNone(fd)
            self.assertEqual(flock.call_count, 2)

    def test_unsupported_locking_runs_without_slot(self):
        error = OSError(errno.ENOLCK, "No locks available")
        stderr = io.StringIO()
        with (
            tempfile.TemporaryDirectory() as lock_dir,
            patch.object(self.script.fcntl, "flock", side_effect=error) as flock,
            patch.object(self.script.os, "execvp") as execvp,
            contextlib.redirect_stderr(stderr),
        ):
            self.script.main(["compile_slot.py", "2", lock_dir, "c++", "-c", "a.cpp"])
        self.assertEqual(flock.call_count, 1)
        execvp.assert_called_once_with("c++", ["c++", "-c", "a.cpp"])
        self.assertEqual(
            stderr.getvalue().splitlines(),
            [
                f"compile_slot: file locking unavailable ({error}), running without a slot limit"
            ],
        )


if __name__ == "__main__":
    unittest.main()
