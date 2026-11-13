"""Tests of pharesee.amr.io.rawio: threaded range reads and file identity."""

import os
import subprocess
import sys
import unittest
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.errors import PhareseeFileChanged, PhareseeIOError
from pyphare.pharesee.amr.io import rawio
from pyphare_tests.test_pharesee.amr._common import tmp_dir


class RawIOTest(unittest.TestCase):
    def setUp(self):
        path = tmp_dir(self) / "data.bin"
        self.data = np.random.default_rng(12).integers(0, 256, size=3_000_000, dtype=np.uint8)
        self.data.tofile(path)
        self.fd = rawio.open_readonly(path)
        self.addCleanup(os.close, self.fd)

    def ranges(self):
        offsets = np.array([2_000_000, 7, 0, 1_234_567, 2_999_000], dtype=np.int64)
        lengths = np.array([700_001, 3, 0, 1_000_000, 1000], dtype=np.int64)
        return offsets, lengths

    def expected(self, offsets, lengths):
        return [self.data[o : o + n] for o, n in zip(offsets, lengths)]

    def test_new_buffer_concatenates_ranges_in_given_order(self):
        offsets, lengths = self.ranges()
        for nthreads in (1, 3, 8):
            out = rawio.read_ranges(self.fd, offsets, lengths, nthreads=nthreads, block=4096)
            np.testing.assert_array_equal(out, np.concatenate(self.expected(offsets, lengths)))

    def test_one_buffer_in_place(self):
        offsets, lengths = self.ranges()
        buf = bytearray(int(lengths.sum()) + 10)
        ret = rawio.read_ranges(self.fd, offsets, lengths, buf, block=100_000)
        self.assertIs(ret, buf)
        got = np.frombuffer(buf, np.uint8)[: int(lengths.sum())]
        np.testing.assert_array_equal(got, np.concatenate(self.expected(offsets, lengths)))
        self.assertEqual(bytes(buf[-10:]), bytes(10))

    def test_list_of_buffers_and_typed_arrays(self):
        offsets, lengths = self.ranges()
        outs = [np.zeros(int(n), np.uint8) for n in lengths]
        rawio.read_ranges(self.fd, offsets, lengths, outs, nthreads=4, block=65536)
        for got, ref in zip(outs, self.expected(offsets, lengths)):
            np.testing.assert_array_equal(got, ref)
        # a float64 destination receives the raw bytes
        f64 = np.empty(1000, np.float64)
        rawio.read_ranges(self.fd, [8000], [8000], f64)
        np.testing.assert_array_equal(f64.view(np.uint8), self.data[8000:16000])

    def test_bad_requests_raise(self):
        cases = [  # (error, offsets, lengths, destination, keywords)
            (PhareseeIOError, [2_999_990], [100], None, dict(nthreads=2)),  # past the end
            (PhareseeIOError, [0, 2_999_000], [10, 5000], None, dict(nthreads=4, block=1000)),
            (ValueError, [0], [10], bytes(10), {}),  # read-only
            (ValueError, [0], [10], bytearray(5), {}),  # too small
            (ValueError, [0, 1], [10], None, {}),
            (ValueError, [0], [10], [bytearray(10), bytearray(10)], {}),
        ]
        for error, offsets, lengths, out, kw in cases:
            with self.subTest(offsets=offsets, out=type(out)), self.assertRaises(error):
                rawio.read_ranges(self.fd, offsets, lengths, out, **kw)

    def test_io_threads_env(self):
        for value, expect in (("3", 3), ("0", 1), (" ", 8)):
            with mock.patch.dict(os.environ, {"PHARESEE_IO_THREADS": value}):
                self.assertEqual(rawio.io_threads(), expect)
        for bad in ("x", "abc", "2.5", "8 threads"):
            with mock.patch.dict(os.environ, {"PHARESEE_IO_THREADS": bad}), \
                    self.assertRaisesRegex(ValueError, "PHARESEE_IO_THREADS.*whole number"):
                rawio.io_threads()
        env = dict(os.environ)
        env.pop("PHARESEE_IO_THREADS", None)
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertEqual(rawio.io_threads(), 8)


class EnvTest(unittest.TestCase):
    def test_flag(self):
        for value, expect in (("1", True), ("TRUE", True), ("yes", True), ("On", True),
                              ("0", False), ("false", False), ("No", False), ("off", False),
                              ("", True)):
            with mock.patch.dict(os.environ, {"PHARESEE_X": value}):
                self.assertIs(rawio.env_flag("PHARESEE_X", True), expect, value)
        for bad in ("2", "maybe", "y"):
            with mock.patch.dict(os.environ, {"PHARESEE_X": bad}), \
                    self.assertRaisesRegex(ValueError, "PHARESEE_X.*true.*false"):
                rawio.env_flag("PHARESEE_X", True)
        with mock.patch.dict(os.environ, {"PHARESEE_RAW_IO": "nope"}), \
                self.assertRaisesRegex(ValueError, "PHARESEE_RAW_IO"):
            rawio.raw_io_enabled()

    def test_numbers(self):
        with mock.patch.dict(os.environ, {"PHARESEE_X": "2.5"}):
            self.assertEqual(rawio.env_float("PHARESEE_X", 1.0), 2.5)
            with self.assertRaisesRegex(ValueError, "PHARESEE_X"):
                rawio.env_int("PHARESEE_X", 1)
        with mock.patch.dict(os.environ, {"PHARESEE_X": "lots"}), \
                self.assertRaisesRegex(ValueError, "PHARESEE_X.*number"):
            rawio.env_float("PHARESEE_X", 1.0)
        env = {k: v for k, v in os.environ.items() if k != "PHARESEE_X"}
        with mock.patch.dict(os.environ, env, clear=True):
            self.assertEqual(rawio.env_int("PHARESEE_X", 7), 7)
            self.assertEqual(rawio.env_float("PHARESEE_X", 7.5), 7.5)


class TimeMatchTest(unittest.TestCase):

    def test_nearest_time(self):
        times = np.array([0.0, 0.2, 0.4, 10.0])
        self.assertEqual(rawio.nearest_time(times, np.float32(0.2)), 1)
        self.assertEqual(rawio.nearest_time(times, np.float32(10.0)), 3)
        self.assertEqual(rawio.nearest_time(times, 0.2 + 4e-11), 1)
        self.assertIsNone(rawio.nearest_time(times, float(np.float32(0.2))))
        self.assertIsNone(rawio.nearest_time(times, np.float32(0.3)))
        self.assertIsNone(rawio.nearest_time(times, float("nan")))
        self.assertIsNone(rawio.nearest_time(np.zeros(0), 0.0))
        with self.assertRaises(TypeError):
            rawio.nearest_time(times, "soon")


class IdentityTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tmp_dir(self)
        self.path = self.tmp / "f.bin"
        self.path.write_bytes(b"x" * 1000)

    def test_replaced_file_raises(self):
        ident = rawio.file_identity(self.path)
        other = self.tmp / "g.bin"
        other.write_bytes(b"z" * 2000)
        os.replace(other, self.path)
        with self.assertRaises(PhareseeFileChanged) as cm:
            rawio.check_identity(ident)
        self.assertIn("replaced", str(cm.exception))
        self.assertIsInstance(cm.exception, PhareseeIOError)
        self.assertIsInstance(cm.exception, OSError)

    def test_deleted_file_raises(self):
        ident = rawio.file_identity(self.path)
        os.remove(self.path)
        with self.assertRaises(PhareseeFileChanged):
            rawio.check_identity(ident)


class ImportTest(unittest.TestCase):
    def test_no_thread_created_at_import(self):
        code = (
            "import threading\n"
            "n = threading.active_count()\n"
            "from pyphare.pharesee.amr.io import rawio, btree, vtkhdf\n"
            "print(threading.active_count() - n)\n"
        )
        res = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            env=os.environ,
            check=False,
        )
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(res.stdout.strip(), "0")


if __name__ == "__main__":
    unittest.main()
