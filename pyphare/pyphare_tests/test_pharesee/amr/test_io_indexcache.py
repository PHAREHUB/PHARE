#!/usr/bin/env python3
"""Tests of the on-disk phareh5 index cache (pyphare.pharesee.amr.io.indexcache)
and of its use by H5Backend.index()."""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyphare.pharesee.amr.io import indexcache, phareh5
from pyphare.pharesee.amr.io.phareh5 import H5Backend
from pyphare_tests.test_pharesee.amr.h5_files import tiling, write_phareh5

_ENV = ("PHARESEE_CACHE_DIR", "PHARESEE_CACHE", "PHARESEE_CACHE_MB", "XDG_CACHE_HOME")


class CacheBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self._env = {k: os.environ.get(k) for k in _ENV}
        for k in _ENV:
            os.environ.pop(k, None)
        self.cache = self.tmp / "cache"
        os.environ["PHARESEE_CACHE_DIR"] = str(self.cache)
        self.data = self.tmp / "data"
        self.data.mkdir()

    def tearDown(self):
        for k, v in self._env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        self._tmp.cleanup()

    def write(self, name="EM_B.h5", times=(0.0, 1.0)):
        patches = {t: {0: tiling(2, 16, 8), 1: tiling(2, 8, 4, first_id=10)} for t in times}
        return str(write_phareh5(self.data / name, dim=2, patches=patches, domain_cells=[16, 16]))

    def entries(self):
        return sorted(self.cache.glob("*.npz")) if self.cache.exists() else []

    def assert_same_index(self, a, b):
        self.assertEqual(len(a), len(b))
        for x, y in zip(a, b):
            self.assertEqual(x.level, y.level)
            np.testing.assert_array_equal(x.lower, y.lower)
            np.testing.assert_array_equal(x.upper, y.upper)
            for field in ("rank", "localid", "ghosts", "shape", "offset", "nbytes"):
                np.testing.assert_array_equal(getattr(x.payload, field), getattr(y.payload, field))
                self.assertEqual(getattr(y.payload, field).dtype, np.int64)
            self.assertEqual(x.payload.datasets, y.payload.datasets)


class TestModule(CacheBase):
    def test_location(self):
        self.assertEqual(indexcache.cache_dir(), str(self.cache))
        os.environ.pop("PHARESEE_CACHE_DIR")
        os.environ["XDG_CACHE_HOME"] = str(self.tmp / "xdg")
        self.assertEqual(indexcache.cache_dir(), str(self.tmp / "xdg" / "pharesee"))
        os.environ.pop("XDG_CACHE_HOME")
        self.assertEqual(
            indexcache.cache_dir(), os.path.join(os.path.expanduser("~"), ".cache", "pharesee")
        )

    def test_roundtrip_and_key(self):
        arrays = {"a": np.arange(5, dtype=np.int64), "names": np.array(["EM_B_x", "EM_B_y"])}
        self.assertTrue(indexcache.store("/x/f.h5", 7, "0.0000000000", 96, "d1", arrays))
        got = indexcache.load("/x/f.h5", 7, "0.0000000000", 96, "d1")
        np.testing.assert_array_equal(got["a"], arrays["a"])
        np.testing.assert_array_equal(got["names"], arrays["names"])
        for other in (
            ("/x/g.h5", 7, "0.0000000000", 96, "d1"),
            ("/x/f.h5", 8, "0.0000000000", 96, "d1"),
            ("/x/f.h5", 7, "1.0000000000", 96, "d1"),
            ("/x/f.h5", 7, "0.0000000000", 97, "d1"),
            ("/x/f.h5", 7, "0.0000000000", 96, "d2"),  # header changed
        ):
            self.assertIsNone(indexcache.load(*other), other)

    def test_lru_cap(self):
        big = {"a": np.zeros(40_000, dtype=np.int64)}  # ~320 kB per entry
        os.environ["PHARESEE_CACHE_MB"] = "0.5"  # room for one entry only
        indexcache.store("/f", 1, "t0", 1, "d", big)
        first = self.entries()
        self.assertEqual(len(first), 1)
        indexcache.store("/f", 1, "t1", 2, "d", big)
        self.assertEqual(len(self.entries()), 1)
        self.assertIsNone(indexcache.load("/f", 1, "t0", 1, "d"))
        self.assertIsNotNone(indexcache.load("/f", 1, "t1", 2, "d"))

    def test_bad_environment_values(self):
        os.environ["PHARESEE_CACHE"] = "sometimes"
        with self.assertRaisesRegex(ValueError, "PHARESEE_CACHE="):
            indexcache.enabled()
        os.environ["PHARESEE_CACHE"] = "off"
        self.assertFalse(indexcache.enabled())
        os.environ.pop("PHARESEE_CACHE")
        os.environ["PHARESEE_CACHE_MB"] = "lots"
        with self.assertRaisesRegex(ValueError, "PHARESEE_CACHE_MB"):
            indexcache.cap_bytes()
        with self.assertRaisesRegex(ValueError, "PHARESEE_CACHE_MB"):
            indexcache.store("/f", 1, "t", 1, "d", {"a": np.zeros(3)})

    def test_disabled(self):
        os.environ["PHARESEE_CACHE"] = "0"
        self.assertFalse(indexcache.store("/f", 1, "t", 1, "d", {"a": np.zeros(3)}))
        self.assertFalse(self.cache.exists())
        self.assertIsNone(indexcache.load("/f", 1, "t", 1, "d"))


class TestBackendCache(CacheBase):
    def test_second_backend_does_not_walk(self):
        path = self.write()
        first = H5Backend(path).index(1.0)
        walks = phareh5.STATS["walks"]
        self.assertEqual(len(self.entries()), 1)
        second = H5Backend(path).index(1.0)
        self.assertEqual(phareh5.STATS["walks"], walks)
        self.assert_same_index(first, second)
        # the cached index reads the same data
        b = H5Backend(path)
        li = b.index(1.0)[1]
        nn = li.upper - li.lower + 2
        size = np.prod(nn, axis=1)
        out = {n: np.empty(size.sum(), np.float32) for n in "xyz"}
        b.read_into(1.0, 1, np.arange(len(size)), out, np.cumsum(size) - size)
        self.assertEqual(phareh5.STATS["walks"], walks)
        self.assertFalse(np.isnan(out["x"]).any())
        # nothing written next to the data
        self.assertEqual(sorted(p.name for p in self.data.iterdir()), ["EM_B.h5"])

    def test_corrupted_entry_is_a_miss(self):
        path = self.write()
        ref = H5Backend(path).index(0.0)
        (entry,) = self.entries()
        for damage in (b"", b"not a zip file" * 10, entry.read_bytes()[:200]):
            entry.write_bytes(damage)
            walks = phareh5.STATS["walks"]
            got = H5Backend(path).index(0.0)
            self.assertEqual(phareh5.STATS["walks"], walks + 1)
            self.assert_same_index(ref, got)

    def test_changed_time_group_is_a_miss(self):
        import h5py

        path = self.write()
        H5Backend(path).index(0.0)
        inode = os.stat(path).st_ino
        # rewrite the same file in place: the t/<time> header changes, not the inode
        with h5py.File(path, "r+") as f:
            f["t/0.0000000000"].attrs["rewritten"] = np.arange(64)
            del f["t/0.0000000000/pl1/p0#13"]
        self.assertEqual(os.stat(path).st_ino, inode)
        walks = phareh5.STATS["walks"]
        got = H5Backend(path).index(0.0)
        self.assertEqual(phareh5.STATS["walks"], walks + 1)
        self.assertEqual(got[1].lower.shape[0], 3)

    def test_rewritten_file_same_header_address(self):
        """Same path, same inode, same structure, other boxes: the header hash differs."""
        path = self.write()
        a = H5Backend(path).index(0.0)
        with open(path, "r+b") as fh:  # keep the inode
            fh.truncate(0)
        patches = {t: {0: tiling(2, 16, 4), 1: tiling(2, 8, 4, first_id=10)} for t in (0.0, 1.0)}
        tmp = write_phareh5(self.tmp / "new.h5", dim=2, patches=patches, domain_cells=[16, 16])
        with open(path, "r+b") as fh:
            fh.write(Path(tmp).read_bytes())
        b = H5Backend(path).index(0.0)
        self.assertEqual(a[0].lower.shape[0], 4)
        self.assertEqual(b[0].lower.shape[0], 16)

    def test_unwritable_cache_dir(self):
        blocker = self.tmp / "file"
        blocker.write_text("x")
        os.environ["PHARESEE_CACHE_DIR"] = str(blocker / "cache")  # parent is a file
        path = self.write()
        self.assertEqual(H5Backend(path).index(0.0)[0].lower.shape[0], 4)


if __name__ == "__main__":
    unittest.main()
