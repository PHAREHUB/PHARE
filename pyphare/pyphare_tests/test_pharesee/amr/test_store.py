#!/usr/bin/env python3
"""Tests of pharesee.amr.store: lazy buffers and LevelData reading through a backend."""

import gc
import os
import pickle
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.errors import PhareseeIOError
from pyphare.pharesee.amr.geometry import LevelGeometry
from pyphare.pharesee.amr.io import open_backend
from pyphare.pharesee.amr.store import (
    LevelData,
    lazy_empty,
    open_from_spec,
)
from pyphare_tests.test_pharesee.amr.h5_files import write_phareh5
from pyphare_tests.test_pharesee.amr.testing import FakeBackend, fill, make_amr
from pyphare_tests.test_pharesee.amr.vtkhdf_files import (
    tiling,
    write_vtkhdf,
)

PAGE = os.sysconf("SC_PAGE_SIZE")


def rss():
    with open("/proc/self/statm") as f:
        return int(f.read().split()[1]) * PAGE


@unittest.skipUnless(sys.platform.startswith("linux"), "needs /proc")
class TestLazyEmpty(unittest.TestCase):
    def test_scattered_patches_cost_their_bytes(self):
        n = 2**28  # 2 GiB of float64, address space only
        patch = 6900
        starts = np.random.default_rng(0).choice(n // patch, 300, replace=False) * patch
        gc.collect()
        r0 = rss()
        a = lazy_empty(n, np.float64)
        self.assertLess(rss() - r0, 4 * 2**20, "an untouched buffer must cost ~nothing")
        for s in starts:
            a[s : s + patch] = 1.0
        touched = 300 * patch * 8
        grown = rss() - r0
        self.assertGreaterEqual(grown, touched * 0.9)
        self.assertLess(grown, 1.5 * touched + 4 * 2**20, f"{grown} bytes for {touched}")
        del a
        gc.collect()

    def test_properties(self):
        a = lazy_empty(1000, np.float32)
        self.assertEqual((a.shape, a.dtype), ((1000,), np.float32))
        self.assertTrue(a.flags.writeable)
        a[:] = 3
        self.assertIsNotNone(a.base)  # keeps the mapping alive
        self.assertEqual(lazy_empty(0, np.float64).shape, (0,))


class Setup:
    def make(self, dim=2, delay=0.0, nlevels=2, seed=3):
        self.hier = make_amr(dim, nlevels, seed=seed, two_node=True)
        f = {"x": lambda *c: sum(c) + 1.0, "y": lambda *c: c[0] * 2.0}
        self.src = fill(self.hier, f)
        self.fb = FakeBackend(self.src, delay=delay)
        self.lds = [LevelData.from_backend(self.fb, g, 0.0)
                    for g in self.hier.levels]


class TestLevelData(unittest.TestCase, Setup):
    def test_reads_exactly_missing_patches_once(self):
        self.make()
        ld, src, g = self.lds[0], self.src[0], self.hier.levels[0]
        self.assertFalse(ld.loaded.any())
        self.assertEqual(self.fb.calls, 0)
        first = [0, 3, g.npatch - 1]
        ld.ensure(first)
        self.assertEqual(self.fb.calls, 1)
        np.testing.assert_array_equal(self.fb.requests[0][2], np.sort(g.file_index[first]))
        np.testing.assert_array_equal(np.nonzero(ld.loaded)[0], first)
        ld.ensure(first)
        self.assertEqual(self.fb.calls, 1)
        ld.ensure()
        self.assertEqual(self.fb.calls, 2)
        rest = np.setdiff1d(np.arange(g.npatch), first)
        np.testing.assert_array_equal(self.fb.requests[1][2], np.sort(g.file_index[rest]))
        np.testing.assert_array_equal(self.fb.reads[(0.0, 0)], 1)
        for name in ("x", "y"):
            np.testing.assert_array_equal(ld.buffer(name), src.buffer(name))
        self.assertEqual(self.fb.calls, 2)
        self.assertEqual(ld.loaded_nbytes(), ld.nbytes())

    def test_view_reads_one_patch(self):
        self.make(dim=3)
        ld, src = self.lds[1], self.src[1]
        v = ld.view(2, "y")
        np.testing.assert_array_equal(v, src.view(2, "y"))
        self.assertEqual(int(ld.loaded.sum()), 0)  # x and z not read: no patch complete
        self.assertEqual(len(ld.missing(None, "y")[0]), ld.geom.npatch - 1)
        self.assertEqual(self.fb.names_read, [("y",)])
        self.assertEqual(v.shape, tuple(self.hier.levels[1].shape[2]))

    def test_failed_read_leaves_loaded_false(self):
        self.make()
        ld = self.lds[0]
        self.fb.fail_next()
        with self.assertRaises(PhareseeIOError):
            ld.ensure()
        self.assertFalse(ld.loaded.any())
        ld.ensure()
        self.assertTrue(ld.loaded.all())
        np.testing.assert_array_equal(self.fb.reads[(0.0, 0)], 2)
        np.testing.assert_array_equal(ld.buffer("x"), self.src[0].buffer("x"))

    def test_unload(self):
        self.make()
        ld = self.lds[1]
        held = ld.view(0, "x")
        ld.unload()
        self.assertFalse(ld.loaded.any())
        np.testing.assert_array_equal(held, self.src[1].view(0, "x"))  # user's view still valid
        ld.view(0, "x")
        self.assertEqual(self.fb.calls, 2)
        self.src[1].unload()  # in memory: nothing happens
        self.assertTrue(self.src[1].loaded.all())

    def test_arguments(self):
        self.make()
        with self.assertRaises(NotImplementedError):
            self.lds[0].ensure(slab=(0, 1))
        with self.assertRaises(KeyError):
            self.lds[0].ensure(names=["z"])
        with self.assertRaises(KeyError):
            self.lds[0].view(0, "z")
        with self.assertRaises(IndexError):
            self.lds[0].ensure([10**6])
        g = self.hier.levels[0]
        with self.assertRaises(ValueError):
            LevelData.from_arrays(g, 0.0, {"x": np.zeros(g.nnodes + 1)})

    def test_unload_and_partial_read_after_a_read_keep_the_buffer_whole(self):
        """Another thread unloads and reads one patch between this thread's read and the
        return of its buffer: buffer() still returns the buffer it filled, every patch
        valid (not the new, partly filled one)."""
        self.make()
        ld, src = self.lds[0], self.src[0]
        ensure = ld.ensure

        def ensure_then_other_thread(slots=None, names=None, slab=None):
            got = ensure(slots, names, slab)
            ld.ensure = ensure  # the other thread runs once
            ld.unload()
            ld.ensure([0])
            return got

        ld.ensure = ensure_then_other_thread
        np.testing.assert_array_equal(ld.buffer("x"), src.buffer("x"))
        ld.ensure = ensure_then_other_thread
        np.testing.assert_array_equal(ld.view(1, "x"), src.view(1, "x"))

    def test_read_only_arrays_cannot_be_made_writeable(self):
        self.make()
        ld = self.lds[0]
        for a in (ld.buffer("x"), ld.view(0, "y"), ld.raw("x"), ld.loaded):
            self.assertFalse(a.flags.writeable)
            with self.assertRaises(ValueError):
                a.flags.writeable = True

    def test_from_arrays_shares_the_given_arrays(self):
        """Documented choice: no copy (in-memory component views share their parent's
        buffers, and every operation result would be copied once more)."""
        self.make()
        g = self.hier.levels[0]
        a = np.zeros(g.nnodes)
        ld = LevelData.from_arrays(g, 0.0, {"value": a})
        a[0] = 7.0
        self.assertEqual(ld.buffer("value")[0], 7.0)

    def test_memory_error_names_the_level_reservation(self):
        self.make()
        ld = self.lds[1]

        def no_memory(n, dtype):
            raise MemoryError("cannot reserve")

        with mock.patch("pyphare.pharesee.amr.store.lazy_empty", no_memory):
            for call in (lambda: ld.ensure([0]), lambda: ld.raw("x")):
                with self.assertRaises(MemoryError) as cm:
                    call()
                msg = str(cm.exception)
                self.assertIn("level 1", msg)
                self.assertIn("x, y", msg)
                self.assertIn(f"{ld.nbytes() / 2**30:.2f} GiB", msg)
                self.assertIn("box=", msg)
        self.assertFalse(ld.loaded.any())
        ld.ensure([0])  # memory back: reads normally
        np.testing.assert_array_equal(ld.view(0, "x"), self.src[1].view(0, "x"))

    def test_concurrent_ensure_reads_each_patch_once(self):
        self.make(delay=0.02)
        ld = self.lds[0]
        n = self.hier.levels[0].npatch
        rng = np.random.default_rng(0)
        subsets = [rng.choice(n, n // 2, replace=False) for _ in range(7)] + [None]
        errors = []

        def work(s):
            try:
                ld.ensure(s)
            except Exception as e:  # noqa: BLE001  (reported below)
                errors.append(e)

        threads = [threading.Thread(target=work, args=(s,)) for s in subsets]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertEqual(errors, [])
        self.assertTrue(ld.loaded.all())
        self.assertEqual(int(self.fb.reads[(0.0, 0)].max()), 1)
        np.testing.assert_array_equal(ld.buffer("x"), self.src[0].buffer("x"))


class TestRaw(unittest.TestCase, Setup):
    def test_raw_never_reads(self):
        self.make()
        ld, g = self.lds[0], self.hier.levels[0]
        r = ld.raw("x")  # nothing read yet: lazy buffer allocated, no backend call
        self.assertEqual(self.fb.calls, 0)
        self.assertEqual(r.shape, (g.nnodes,))
        self.assertFalse(r.flags.writeable)
        ld.ensure([2])
        self.assertEqual(self.fb.calls, 1)
        r2 = ld.raw("x")
        self.assertTrue(np.shares_memory(r, r2))  # ensure filled the buffer raw handed out
        a, b = g.offset[2], g.offset[3]
        np.testing.assert_array_equal(r[a:b], self.src[0].buffer("x")[a:b])
        self.assertEqual(int(ld.loaded.sum()), 1)
        self.assertEqual(self.fb.calls, 1)
        with self.assertRaises(KeyError):
            ld.raw("z")
        mem = self.src[0].raw("y")  # in memory: the array itself, writeable
        self.assertIs(mem, self.src[0].buffer("y"))
        self.assertTrue(mem.flags.writeable)


class _DerivedFake(FakeBackend):
    """FakeBackend with read_derived: copies every component (records the call)."""

    def read_derived(self, kind, time, level, file_ids, out, dst_start):
        self.derived_calls = getattr(self, "derived_calls", 0) + 1
        self.read_into(time, level, file_ids, out, dst_start)


class TestComponentReads(unittest.TestCase):
    """Only the components asked for are read; loaded state is per (component, patch)."""

    def make(self, backend=FakeBackend, derived=None):
        self.hier = make_amr(2, 1, seed=5, two_node=True)
        f = {"x": lambda *c: sum(c) + 1.0, "y": lambda *c: c[0] * 2.0,
             "z": lambda *c: c[1] - 3.0}
        self.src = fill(self.hier, f)[0]
        self.fb = backend(fill(self.hier, f))
        return LevelData.from_backend(self.fb, self.hier.levels[0], 0.0,
                                      derived=derived)

    def check_values(self, ld, names):
        for n in names:
            np.testing.assert_array_equal(ld.buffer(n), self.src.buffer(n))

    def test_z_then_x_then_rest(self):
        ld = self.make()
        ld.ensure(names="z")
        self.assertEqual(self.fb.names_read, [("z",)])
        self.assertFalse(ld.loaded.any())  # no patch holds every component
        self.assertEqual(len(ld.missing(None, "z")[0]), 0)
        self.assertEqual(ld.missing(None, ["x", "z"])[1], ("x",))
        self.assertEqual(ld.loaded_nbytes(), ld.nbytes() // 3)
        ld.ensure(names=["x"])
        self.assertEqual(self.fb.names_read, [("z",), ("x",)])
        ld.ensure()  # every component: only y is still missing
        self.assertEqual(self.fb.names_read, [("z",), ("x",), ("y",)])
        self.assertTrue(ld.loaded.all())
        ld.ensure(names=["x", "y", "z"])
        self.assertEqual(self.fb.calls, 3)
        self.check_values(ld, "xyz")

    def test_buffer_and_view_read_their_component(self):
        ld = self.make()
        np.testing.assert_array_equal(ld.view(1, "y"), self.src.view(1, "y"))
        self.assertEqual(self.fb.names_read, [("y",)])
        np.testing.assert_array_equal(self.fb.requests[0][2], [ld.geom.file_index[1]])
        ld.buffer("z")
        self.assertEqual(self.fb.names_read[-1], ("z",))
        ld.raw("x")  # never reads
        self.assertEqual(self.fb.calls, 2)

    def test_buffers_reads_once(self):
        ld = self.make()
        ld.ensure([0], "y")
        bufs = ld.buffers()
        self.assertEqual(list(bufs), ["x", "y", "z"])
        self.assertEqual(self.fb.names_read, [("y",), ("x", "y", "z")])
        self.assertEqual(list(ld.buffers(["z", "x"])), ["z", "x"])
        self.assertEqual(self.fb.calls, 2)
        self.check_values(ld, "xyz")

    def test_mixed_missing_patches_read_in_one_call(self):
        ld = self.make()
        ld.ensure([0], "x")
        ld.ensure([1], "y")
        ld.ensure([0, 1, 2], ["x", "y"])
        self.assertEqual(self.fb.calls, 3)
        self.assertEqual(self.fb.names_read[-1], ("x", "y"))
        np.testing.assert_array_equal(self.fb.requests[-1][2],
                                      np.sort(ld.geom.file_index[[0, 1, 2]]))
        for slot in (0, 1, 2):
            for n in "xy":
                np.testing.assert_array_equal(ld.view(slot, n), self.src.view(slot, n))
        self.assertEqual(self.fb.calls, 3)

    def test_component_level_reads_its_component(self):
        from pyphare.pharesee.amr.field import ComponentLevel, _missing_bytes

        ld = self.make()
        cz = ComponentLevel(ld, "z")
        self.assertEqual(_missing_bytes(cz, None), ld.nbytes() // 3)
        self.assertEqual(cz.read_nbytes(), ld.read_nbytes(None, "z"))
        np.testing.assert_array_equal(cz.buffer("value"), self.src.buffer("z"))
        self.assertEqual(self.fb.names_read, [("z",)])
        self.assertTrue(cz.loaded.all())
        self.assertFalse(ld.loaded.any())
        self.assertEqual(_missing_bytes(cz, None), 0)
        self.assertEqual(cz.loaded_nbytes(), cz.nbytes())
        cx = ComponentLevel(ld, "x")
        self.assertFalse(cx.loaded.any())
        cx.view(0, "value")
        self.assertEqual(self.fb.names_read, [("z",), ("x",)])
        np.testing.assert_array_equal(cx._loaded, np.arange(ld.geom.npatch) == 0)
        ld.unload()  # clears every component
        self.assertFalse(cz.loaded.any())
        self.assertEqual(ld.loaded_nbytes(), 0)
        ld.ensure()
        self.assertEqual(self.fb.names_read[-1], ("x", "y", "z"))

    def test_read_nbytes_counts_wanted_components(self):
        ld = self.make()
        one = ld.read_nbytes(None, "z")
        self.assertEqual(one * 3, ld.read_nbytes())
        ld.ensure(names="z")
        self.assertEqual(ld.read_nbytes(), 2 * one)

    def test_derived_reads_every_component(self):
        ld = self.make(backend=_DerivedFake, derived="J_yee")
        ld.ensure(names="z")
        self.assertEqual(self.fb.derived_calls, 1)
        self.assertEqual(self.fb.names_read, [("x", "y", "z")])
        self.assertTrue(ld.loaded.all())
        self.assertEqual(ld.read_nbytes(None, "x"), 0)


class TestRealBackendPickle(unittest.TestCase):
    """File-backed LevelData of real vtkhdf / phareh5 files survive a pickle round trip."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        env = {"PHARESEE_CACHE_DIR": str(self.tmp / "cache")}
        self._env = mock.patch.dict(os.environ, env)
        self._env.start()

    def tearDown(self):
        self._env.stop()
        self._tmp.cleanup()

    def level(self, backend, file_ids=None):
        """File-backed LevelData of level 0 at the first time, optionally on a patch subset."""
        meta = backend.meta()
        t = float(backend.times()[0])
        idx = next(li for li in backend.index(t) if li.level == 0)
        lower, upper = idx.lower, idx.upper
        if file_ids is not None:
            lower, upper = lower[file_ids], upper[file_ids]
        geom = LevelGeometry.for_level(0, lower, upper, meta.dl0, meta.domain_cells0,
                                       meta.periodic)
        return LevelData.from_backend(backend, geom, t, file_ids=file_ids)

    def check_round_trip(self, path, fmt):
        backend = open_backend(path)
        spec = backend.spec()
        self.assertEqual(spec[0], fmt)
        self.assertEqual(type(open_from_spec(spec)), type(backend))
        whole = self.level(backend)
        self.assertGreater(len(set(map(tuple, whole.geom.shape))), 1, "want unequal patches")
        values = {n: np.array(whole.buffer(n)) for n in whole.names}
        back = pickle.loads(pickle.dumps(whole))
        self.assertIsNotNone(back.backend)
        self.assertFalse(back.loaded.any())
        self.assertEqual(back.geom, whole.geom)
        for n in whole.names:
            np.testing.assert_array_equal(back.buffer(n), values[n])
        # a subset of the patches (file ids not in increasing store order on purpose)
        ids = np.array([1, whole.geom.npatch - 1, 3], dtype=np.int64)
        sub = self.level(backend, file_ids=ids)
        sub_back = pickle.loads(pickle.dumps(sub))
        np.testing.assert_array_equal(sub_back.file_ids, ids)
        self.assertEqual(sub_back.geom, sub.geom)
        for k, fid in enumerate(ids):
            for n in whole.names:
                np.testing.assert_array_equal(
                    sub_back.view(sub.geom.store_index[k], n),
                    whole.view(whole.geom.store_index[fid], n),
                )
        self.assertEqual(int(sub_back.loaded.sum()), len(ids))

    def test_vtkhdf(self):
        path = write_vtkhdf(self.tmp / "f.vtkhdf", [[tiling([14, 10], [4, 4])]],
                            domain_cells=[14, 10])
        self.check_round_trip(str(path), "vtkhdf")

    def test_phareh5(self):
        boxes = tiling([14, 10], [4, 4])
        patches = {0.0: {0: [(0, i, lo, up) for i, (lo, up) in enumerate(boxes)]}}
        path = write_phareh5(self.tmp / "f.h5", dim=2, patches=patches, domain_cells=[14, 10],
                             dtype=np.float64)
        self.check_round_trip(str(path), "phareh5")

    def test_bad_spec(self):
        with self.assertRaises(ValueError):
            open_from_spec(())
        with self.assertRaises(ValueError):
            open_from_spec(("nonsense", 1, 2, 3))


if __name__ == "__main__":
    unittest.main()
