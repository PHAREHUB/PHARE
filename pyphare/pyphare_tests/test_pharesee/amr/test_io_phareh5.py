#!/usr/bin/env python3
"""Tests of the phareh5 reader (pyphare.pharesee.amr.io.phareh5.H5Backend).

Synthetic files come from h5_files.write_phareh5 (structure of real files,
linear dyadic fields: node values must be exact). Real-data tests use
build/sample_data (an error when absent; the 3-rank run: test_mpi3.py): h5 vs vtkhdf of the same run, and the
new reader vs master's Run (``legacy=True``) on the same h5 files.
"""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.errors import (
    FormatError,
    PhareseeFileChanged,
    PhareseeIOError,
    TimeNotFoundError,
)
from pyphare.pharesee.amr.io import phareh5, rawio
from pyphare.pharesee.amr.io.phareh5 import H5Backend
from pyphare_tests.test_pharesee.amr.h5_files import (
    DATASETS,
    linear_coeffs,
    linear_value,
    node_coords,
    tiling,
    write_phareh5,
)
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES, require

SAMPLE_CASES = ("1d", "2d_static", "2d_tagged", "3d")


def layout_offsets(lower, upper):
    """Node count of each patch and its start in a flat buffer (patches back to back)."""
    nn = upper - lower + 2
    size = np.prod(nn, axis=1)
    return nn, np.cumsum(size) - size, int(size.sum())


def read_level(backend, time, li, dtype=None, ids=None):
    """Whole level (or patches `ids`) read back to back: {name: flat}, nn, dst."""
    ids = np.arange(li.lower.shape[0]) if ids is None else np.asarray(ids)
    nn, dst, total = layout_offsets(li.lower[ids], li.upper[ids])
    dtype = backend.meta().dtype if dtype is None else dtype
    out = {n: np.full(total, np.nan, dtype=dtype) for n in backend.meta().names}
    backend.read_into(time, li.level, ids, out, dst)
    return out, nn, dst


def patch_views(flat, nn, dst):
    return [flat[d : d + np.prod(n)].reshape(tuple(n), order="F") for n, d in zip(nn, dst)]


def two_level_patches(dim):
    """Level 0 tiled by 3 patches (one of them 1 cell wide = 2 nodes), level 1 by 2."""
    c = 16
    if dim == 1:
        l0 = [(0, 0, (0,), (6,)), (0, 1, (7,), (7,)), (1, 2, (8,), (c - 1,))]
        l1 = [(0, 3, (4,), (9,)), (1, 4, (10,), (10,))]
    elif dim == 2:
        l0 = [
            (0, 0, (0, 0), (6, c - 1)),
            (0, 1, (7, 0), (7, c - 1)),
            (1, 2, (8, 0), (c - 1, c - 1)),
        ]
        l1 = [(0, 3, (4, 4), (9, 11)), (1, 4, (10, 4), (10, 9))]
    else:
        l0 = [
            (0, 0, (0, 0, 0), (6, c - 1, 5)),
            (0, 1, (7, 0, 0), (7, c - 1, 5)),
            (1, 2, (8, 0, 0), (c - 1, c - 1, 5)),
        ]
        l1 = [(0, 3, (4, 4, 2), (9, 7, 5)), (1, 4, (10, 4, 2), (10, 5, 3))]
    return {0: l0, 1: l1}


def expected_nodes(dataset, lower, upper, dl):
    coords = np.meshgrid(*node_coords(lower, upper, dl), indexing="ij")
    return np.broadcast_to(linear_value(dataset, tuple(coords)), coords[0].shape)


def ulp_close(mine, ref):
    """|mine - ref| <= 1 float32 ulp of ref, NaN where ref is NaN."""
    ref = np.asarray(ref, dtype=np.float64)
    tol = np.spacing(np.abs(ref).astype(np.float32)).astype(np.float64)
    return np.abs(np.asarray(mine, np.float64) - ref) <= tol


class SyntheticBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self._env = {k: os.environ.get(k) for k in ("PHARESEE_CACHE_DIR", "PHARESEE_CACHE")}
        os.environ["PHARESEE_CACHE_DIR"] = str(self.tmp / "cache")
        os.environ.pop("PHARESEE_CACHE", None)

    def tearDown(self):
        for k, v in self._env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        self._tmp.cleanup()

    def write(self, name="f.h5", **kw):
        return str(write_phareh5(self.tmp / name, **kw))


class TestYeeToNodes(SyntheticBase):
    """Linear fields: every 2-point average is exact, node values must be exact."""

    QUANTITIES = (
        ("EM_B", "hybrid"),
        ("EM_E", "hybrid"),
        ("ions_bulkVelocity", "hybrid"),
        ("ions_charge_density", "hybrid"),
        ("ions_momentum_tensor", "hybrid"),
        ("mhd_rho", "mhd"),
        ("mhd_V", "mhd"),
    )

    def test_linear_fields_exact_every_centering(self):
        for dim in (1, 2, 3):
            for quantity, model in self.QUANTITIES:
                with self.subTest(dim=dim, quantity=quantity):
                    path = self.write(
                        f"{quantity}_{dim}.h5", dim=dim, quantity=quantity, model=model,
                        patches={0.0: two_level_patches(dim)}, domain_cells=[16] * dim,
                        nan_outer_ghost=True,
                    )
                    b = H5Backend(path)
                    meta = b.meta()
                    self.assertEqual(meta.model, model)
                    self.assertEqual(meta.interp_order, 1 if model == "hybrid" else None)
                    datasets = DATASETS[quantity]
                    expected_names = {1: ("value",), 3: ("x", "y", "z"), 6: phareh5._TENSOR}
                    self.assertEqual(meta.names, expected_names[len(datasets)])
                    for li in b.index(0.0):
                        out, nn, dst = read_level(b, 0.0, li)
                        dl = meta.dl0 / 2**li.level
                        for name, dname in zip(meta.names, datasets):
                            views = patch_views(out[name], nn, dst)
                            for p, v in enumerate(views):
                                ref = expected_nodes(dname, li.lower[p], li.upper[p], dl)
                                np.testing.assert_array_equal(v, ref, err_msg=f"{dname} patch {p}")

    def test_float64_file_and_float64_out(self):
        path = self.write(dim=2, patches={0.0: two_level_patches(2)}, dtype=np.float64,
                          domain_cells=[16, 16])
        b = H5Backend(path)
        self.assertEqual(b.meta().dtype, np.float64)
        li = b.index(0.0)[1]
        out, nn, dst = read_level(b, 0.0, li)
        for name, dname in zip("xyz", DATASETS["EM_B"]):
            for p, v in enumerate(patch_views(out[name], nn, dst)):
                ref = expected_nodes(dname, li.lower[p], li.upper[p], b.meta().dl0 / 2)
                np.testing.assert_array_equal(v, ref)

    def test_non_dyadic_rounding(self):
        """Arbitrary values: one dual axis = float32 mean, two = float64 mean rounded once."""
        rng = np.random.default_rng(3)
        table = {}

        def noise(dname, coords):
            v = rng.standard_normal(coords[0].shape)
            table[dname] = v.astype(np.float32)
            return v

        path = self.write(dim=2, patches={0.0: {0: [(0, 0, (0, 0), (7, 5))]}}, func=noise,
                          domain_cells=[8, 6])
        b = H5Backend(path)
        out, nn, dst = read_level(b, 0.0, b.index(0.0)[0])
        n = nn[0]
        bz = table["EM_B_z"].astype(np.float64)  # dual, dual, ghosts 2
        x = bz[1 : 1 + n[0] + 1, 1 : 1 + n[1] + 1]
        x = 0.5 * (x[:-1] + x[1:])
        x = 0.5 * (x[:, :-1] + x[:, 1:])
        np.testing.assert_array_equal(patch_views(out["z"], nn, dst)[0], x.astype(np.float32))
        bx = table["EM_B_x"]  # primal, dual
        x = bx[2 : 2 + n[0], 1 : 1 + n[1] + 1]
        x = np.float32(0.5) * (x[:, :-1] + x[:, 1:])
        self.assertEqual(x.dtype, np.float32)
        np.testing.assert_array_equal(patch_views(out["x"], nn, dst)[0], x)


class TestIndex(SyntheticBase):
    def test_patch_order_is_numeric_rank_localid(self):
        # created in an order that is neither numeric nor alphabetical
        plist = [
            (1, 0, (12, 0), (15, 7)),
            (0, 10, (4, 0), (11, 7)),
            (0, 2, (0, 0), (3, 7)),
            (1, 11, (0, 8), (15, 15)),
        ]
        path = self.write(dim=2, patches={0.0: {0: plist}}, domain_cells=[16, 16])
        b = H5Backend(path)
        li = b.index(0.0)[0]
        np.testing.assert_array_equal(li.payload.rank, [0, 0, 1, 1])
        np.testing.assert_array_equal(li.payload.localid, [2, 10, 0, 11])
        np.testing.assert_array_equal(li.lower, [[0, 0], [4, 0], [12, 0], [0, 8]])
        np.testing.assert_array_equal(li.upper, [[3, 7], [11, 7], [15, 7], [15, 15]])
        self.assertEqual(li.lower.dtype, np.int64)
        self.assertEqual(li.payload.offset.dtype, np.int64)
        # data follows the same order
        out, nn, dst = read_level(b, 0.0, li)
        for p, v in enumerate(patch_views(out["x"], nn, dst)):
            np.testing.assert_array_equal(
                v, expected_nodes("EM_B_x", li.lower[p], li.upper[p], b.meta().dl0)
            )

    def test_time_without_level0(self):
        patches = {0.0: two_level_patches(2), 0.5: {1: two_level_patches(2)[1]}}
        path = self.write(dim=2, patches=patches, domain_cells=[16, 16])
        b = H5Backend(path)
        np.testing.assert_array_equal(b.times(), [0.0, 0.5])
        with self.assertRaises(FormatError) as cm:
            b.index(0.5)
        self.assertIn("level 0", str(cm.exception))
        out = {n: np.zeros(1000, np.float32) for n in "xyz"}
        with self.assertRaises(FormatError):
            b.read_into(0.5, 1, [0], out, [0])
        self.assertEqual(len(b.index(0.0)), 2)

    def test_time_matching(self):
        times = [0.0, 0.1, 0.30000000000000004, 12.5]
        path = self.write(dim=1, patches={t: {0: [(0, 0, (0,), (15,))]} for t in times},
                          domain_cells=[16])
        b = H5Backend(path)
        np.testing.assert_array_equal(b.times(), [0.0, 0.1, 0.3, 12.5])
        self.assertEqual(b.time_key(0.3), "0.3000000000")
        self.assertEqual(b.time_key(0.1 + 4e-11), "0.1000000000")
        self.assertEqual(b.time_key("12.5000000000"), "12.5000000000")
        with self.assertRaises(TimeNotFoundError) as cm:
            b.index(0.2)
        self.assertIsInstance(cm.exception, KeyError)
        self.assertIsInstance(cm.exception, ValueError)
        msg = str(cm.exception)
        for near in ("0.1", "0.3", "0.0"):
            self.assertIn(near, msg)
        self.assertNotIn("12.5", msg)
        with self.assertRaises(TimeNotFoundError):
            b.time_key(0.1 + 2e-10)

    def test_shape_mismatch_is_format_error(self):
        path = self.write(dim=2, patches={0.0: {0: [(0, 0, (0, 0), (7, 7))]}}, domain_cells=[8, 8])
        import h5py

        with h5py.File(path, "r+") as f:
            f["t/0.0000000000/pl0/p0#0"].attrs["upper"] = np.array([8, 7], np.int32)
        with self.assertRaises(FormatError):
            H5Backend(path).index(0.0)

    def test_one_ghost_layer_is_refused_for_dual(self):
        path = self.write(dim=1, ghosts=1, patches={0.0: {0: [(0, 0, (0,), (7,))]}},
                          domain_cells=[8])
        b = H5Backend(path)
        out = {n: np.zeros(9, np.float32) for n in "xyz"}
        with self.assertRaises(FormatError):
            b.read_into(0.0, 0, [0], out, [0])
        # primal-only quantities need no ghost at all
        path = self.write("n.h5", dim=1, ghosts=0, quantity="ions_charge_density",
                          patches={0.0: {0: [(0, 0, (0,), (7,))]}}, domain_cells=[8])
        b = H5Backend(path)
        out = {"value": np.zeros(9, np.float32)}
        b.read_into(0.0, 0, [0], out, [0])
        np.testing.assert_array_equal(
            out["value"], expected_nodes("charge_density", [0], [7], [0.5])
        )

    def test_unknown_dataset_name(self):
        path = self.write(dim=1, datasets=("mystery",), patches={0.0: {0: [(0, 0, (0,), (7,))]}},
                          domain_cells=[8])
        with self.assertRaises(FormatError):
            H5Backend(path)


class TestReads(SyntheticBase):
    def make(self, **kw):
        plist = [(0, i, (4 * i, 0), (4 * i + 3, 3 + i)) for i in range(6)]
        return self.write(dim=2, patches={0.0: {0: plist}}, domain_cells=[24, 9], **kw)

    def test_subset_scattered_arbitrary_dst(self):
        b = H5Backend(self.make())
        li = b.index(0.0)[0]
        full, _, dst_all = read_level(b, 0.0, li)
        ids = np.array([4, 1, 3])
        nn = li.upper[ids] - li.lower[ids] + 2
        sizes = np.prod(nn, axis=1)
        # gaps between patches and patches in decreasing position
        dst = np.array([7, 7 + sizes[0] + sizes[2] + 20, 7 + sizes[0] + 5])
        total = int((dst + sizes).max()) + 3
        out = {n: np.full(total, -999.0, np.float32) for n in "xyz"}
        b.read_into(0.0, 0, ids, out, dst)
        written = np.zeros(total, bool)
        for k, i in enumerate(ids):
            written[dst[k] : dst[k] + sizes[k]] = True
            for n in "xyz":
                np.testing.assert_array_equal(
                    out[n][dst[k] : dst[k] + sizes[k]],
                    full[n][dst_all[i] : dst_all[i] + sizes[k]],
                )
        for n in "xyz":
            self.assertTrue((out[n][~written] == -999.0).all())

    def test_component_subsets(self):
        """Only the datasets of the components in `out` are read; same bits as a full read."""
        path = self.make()
        for raw_io in (True, False):
            b = H5Backend(path, raw_io=raw_io)
            li = b.index(0.0)[0]
            full, _, dst = read_level(b, 0.0, li)
            ids = np.arange(li.lower.shape[0])
            total = len(full["x"])
            for names in (("z",), ("x", "z"), ("y",), ("z", "x", "y")):
                out = {n: np.full(total, np.nan, np.float32) for n in names}
                seen = []
                real = b._read_patches

                def spy(key, li_, ids_, comps, convert, real=real, seen=seen):
                    seen.append(tuple(comps))
                    return real(key, li_, ids_, comps, convert)

                with mock.patch.object(b, "_read_patches", side_effect=spy):
                    b.read_into(0.0, 0, ids, out, dst)
                self.assertEqual(seen, [tuple(sorted("xyz".index(n) for n in names))])
                for n in names:
                    np.testing.assert_array_equal(out[n].view(np.uint32),
                                                  full[n].view(np.uint32))

    def test_nbytes_per_component(self):
        b = H5Backend(self.make())
        li = b.index(0.0)[0]
        ids = np.array([0, 2, 5])
        full = b.nbytes(0.0, 0, ids)
        parts = [b.nbytes(0.0, 0, ids, names=(n,)) for n in "xyz"]
        self.assertEqual(sum(parts), full)
        self.assertEqual(parts[2], int(li.payload.nbytes[2, ids].sum()))
        self.assertEqual(b.nbytes(0.0, 0, ids, names=["z", "x"]), parts[0] + parts[2])
        with self.assertRaises(KeyError):
            b.nbytes(0.0, 0, ids, names=("w",))

    def test_empty_selection(self):
        b = H5Backend(self.make())
        out = {n: np.full(4, 5.0, np.float32) for n in "xyz"}
        b.read_into(0.0, 0, np.zeros(0, np.int64), out, np.zeros(0, np.int64))
        self.assertTrue(all((v == 5.0).all() for v in out.values()))
        self.assertEqual(b.nbytes(0.0, 0, []), 0)

    def test_bad_requests(self):
        b = H5Backend(self.make())
        out = {n: np.zeros(1000, np.float32) for n in "xyz"}
        with self.assertRaises(IndexError):
            b.read_into(0.0, 0, [6], out, [0])
        with self.assertRaises(ValueError):
            b.read_into(0.0, 0, [0], {n: np.zeros(3, np.float32) for n in "xyz"}, [0])
        with self.assertRaises(KeyError):
            b.read_into(0.0, 0, [0], {"x": out["x"], "w": out["y"]}, [0])  # unknown name
        with self.assertRaises(ValueError):
            b.read_into(0.0, 0, [0], {}, [0])  # no component
        with self.assertRaises(ValueError):
            b.read_into(0.0, 3, [0], out, [0])

    def test_raw_read_equals_read_direct(self):
        path = self.make()
        rng = np.random.default_rng(0)
        chunked = self.write(
            "chunked.h5", dim=2, chunked=True, domain_cells=[24, 9],
            patches={0.0: {0: [(0, i, (4 * i, 0), (4 * i + 3, 3 + i)) for i in range(6)]}},
        )
        raw = H5Backend(path)
        direct = H5Backend(path, raw_io=False)
        fallback = H5Backend(chunked)
        self.assertTrue((fallback.index(0.0)[0].payload.offset == -1).all())
        li = raw.index(0.0)[0]
        ids = np.sort(rng.choice(6, 4, replace=False))
        a, _, _ = read_level(raw, 0.0, li, ids=ids)
        b, _, _ = read_level(direct, 0.0, li, ids=ids)
        c, _, _ = read_level(fallback, 0.0, fallback.index(0.0)[0], ids=ids)
        for n in "xyz":
            self.assertEqual(a[n].tobytes(), b[n].tobytes())
            self.assertEqual(a[n].tobytes(), c[n].tobytes())

    def test_small_batches_and_prefetch(self):
        """Many batches (pipeline) and the prefetch give the same bytes."""
        path = self.make()
        b = H5Backend(path, prefetch="never")
        li = b.index(0.0)[0]
        ref, _, _ = read_level(b, 0.0, li)
        old = phareh5.BATCH_BYTES
        try:
            phareh5.BATCH_BYTES = 1000  # one patch per batch
            b2 = H5Backend(path)  # auto prefetch: whole time requested
            got, _, _ = read_level(b2, 0.0, li)
            self.assertIn("0.0000000000", b2._prefetched)
        finally:
            phareh5.BATCH_BYTES = old
        for n in "xyz":
            np.testing.assert_array_equal(got[n], ref[n])
        self.assertGreaterEqual(b.prefetch(0.0), 0.0)

    def test_prefetch_before_walk(self):
        path = self.make()
        os.environ["PHARESEE_CACHE"] = "0"
        a = H5Backend(path, prefetch="walk")
        ref, _, _ = read_level(H5Backend(path, prefetch="never"), 0.0, a.index(0.0)[0])
        self.assertIn("0.0000000000", a._prefetched)
        got, _, _ = read_level(a, 0.0, a.index(0.0)[0])
        for n in "xyz":
            np.testing.assert_array_equal(got[n], ref[n])
        with self.assertRaises(ValueError):
            H5Backend(path, prefetch="sometimes")

    def test_file_replaced(self):
        path = self.make()
        b = H5Backend(path)
        li = b.index(0.0)[0]
        other = self.write("other.h5", dim=2, domain_cells=[24, 9],
                           patches={0.0: {0: [(0, 0, (0, 0), (3, 3))]}})
        os.replace(other, path)
        out = {n: np.zeros(1000, np.float32) for n in "xyz"}
        with self.assertRaises(PhareseeFileChanged):
            b.read_into(0.0, 0, [0], out, [0])
        self.assertEqual(li.level, 0)


def rewrite_in_place(path, source):
    """Overwrite `path` with the bytes of `source`, keeping the inode."""
    with open(source, "rb") as src, open(path, "r+b") as dst:
        data = src.read()
        dst.truncate(0)
        dst.write(data)


def other_values(dataset, coords):
    """Another quantity on the same layout."""
    return 2.0 * linear_value(dataset, coords) + 7.0


class TestRewriteInPlace(SyntheticBase):
    """A1: an index made before an in-place rewrite (same inode) is never used on the new
    content; a pure append keeps working."""

    def patches(self, times, size=8):
        return {t: {0: tiling(2, [32, 16], [size, 8]), 1: tiling(2, [16, 8], [8, 8])}
                for t in times}

    def check_stale_raises(self, b, times):
        for t in times:
            with self.assertRaisesRegex(PhareseeFileChanged, "rewritten.*new Run"):
                b.index(t)
            out = {n: np.zeros(10**4, np.float32) for n in "xyz"}
            with self.assertRaises(PhareseeFileChanged):
                b.read_into(t, 0, [0], out, [0])
            self.assertFalse(out["x"].any())  # nothing was written
        self.assertEqual(b._indexes, {})

    def write_pair(self, new_kw, cache):
        os.environ["PHARESEE_CACHE"] = "1" if cache else "0"
        path = self.write("q.h5", dim=2, patches=self.patches([0.0, 0.5]),
                          domain_cells=[32, 16])
        new = self.write("new.h5", dim=2, domain_cells=[32, 16], **new_kw)
        return path, new

    def test_other_boxes(self):
        for cache in (False, True):
            path, new = self.write_pair(dict(patches=self.patches([0.0, 0.5], size=4)), cache)
            b = H5Backend(path)
            b.index(0.5)  # a field was created: index cached, nothing read
            rewrite_in_place(path, new)
            self.check_stale_raises(b, (0.5, 0.0))
            fresh = H5Backend(path)
            li = fresh.index(0.5)[0]
            self.assertEqual(len(li.lower), len(tiling(2, [32, 16], [4, 8])))
            got, _, _ = read_level(fresh, 0.5, li)
            ref, _, _ = read_level(H5Backend(new), 0.5, H5Backend(new).index(0.5)[0])
            for n in "xyz":
                np.testing.assert_array_equal(got[n], ref[n])

    def test_same_size_other_data(self):
        for cache, warm in ((False, False), (True, True), (True, False)):
            path, new = self.write_pair(
                dict(patches=self.patches([0.0, 0.5]), func=other_values), cache)
            self.assertEqual(os.path.getsize(path), os.path.getsize(new))
            b = H5Backend(path)
            li = b.index(0.5)[0]
            if warm:
                read_level(b, 0.5, li)
            os.utime(path, ns=(0, 10**18))  # make sure the mtime differs
            rewrite_in_place(path, new)
            self.check_stale_raises(b, (0.5,))
            fresh = H5Backend(path)
            got, _, _ = read_level(fresh, 0.5, fresh.index(0.5)[1])
            want, _, _ = read_level(H5Backend(new), 0.5, H5Backend(new).index(0.5)[1])
            for n in "xyz":
                np.testing.assert_array_equal(got[n], want[n])

    def test_append_keeps_working(self):
        os.environ["PHARESEE_CACHE"] = "0"
        path = self.write("a.h5", dim=2, patches=self.patches([0.0, 0.5]), domain_cells=[32, 16])
        ref = self.write("ref.h5", dim=2, patches=self.patches([0.0, 0.5, 1.0]),
                         domain_cells=[32, 16])
        b = H5Backend(path)
        before = {t: read_level(b, t, b.index(t)[0])[0] for t in (0.5,)}
        identity = rawio.file_identity(path)
        import h5py

        with h5py.File(path, "a") as f, h5py.File(ref, "r") as g:  # the simulation dumps again
            g.copy(g["t/1.0000000000"], f["t"], name="1.0000000000")
        now = rawio.file_identity(path)
        self.assertTrue(rawio.modified_since(identity, now))
        self.assertEqual(now.st_ino, identity.st_ino)
        for t in (0.5, 0.0):  # cached, and indexed for the first time after the append
            for li in b.index(t):
                got, _, _ = read_level(b, t, li)
                want, _, _ = read_level(H5Backend(ref), t, H5Backend(ref).index(t)[li.level])
                for n in "xyz":
                    np.testing.assert_array_equal(got[n], want[n])
        for n in "xyz":
            np.testing.assert_array_equal(read_level(b, 0.5, b.index(0.5)[0])[0][n],
                                          before[0.5][n])
        self.assertEqual(b.times().tolist(), [0.0, 0.5])
        self.assertEqual(H5Backend(path).times().tolist(), [0.0, 0.5, 1.0])


class TestDamagedFiles(SyntheticBase):
    """A9: empty or truncated files give PhareseeIOError / FormatError naming the file."""

    def test_truncated_at_open(self):
        full = self.write("full.h5", dim=2, patches={0.0: {0: tiling(2, [16, 8], [4, 4])}},
                          domain_cells=[16, 8])
        data = Path(full).read_bytes()
        size = len(data)
        for n in (0, 8, 100, 2048, size // 10, size // 2, 9 * size // 10, size - 1):
            p = str(self.tmp / f"t{n}.h5")
            Path(p).write_bytes(data[:n])
            with self.subTest(n=n):
                with self.assertRaises((PhareseeIOError, FormatError)) as cm:
                    b = H5Backend(p)
                    for li in b.index(0.0):
                        read_level(b, 0.0, li)
                self.assertIn(p, str(cm.exception))
                if n == 0:
                    self.assertIn("empty", str(cm.exception))

    def test_file_held_open_by_another_reader(self):
        """Another reader of this session (master's pharesee, plain h5py) keeps the file
        open with h5py's default lock settings: the backend still opens and reads it."""
        import h5py

        p = self.write("held.h5", dim=2, patches={0.0: {0: tiling(2, [16, 8], [4, 4])}},
                       domain_cells=[16, 8])
        b = H5Backend(p)
        expected = read_level(b, 0.0, b.index(0.0)[0])[0]
        with h5py.File(p, "r"):
            b = H5Backend(p)
            got = read_level(b, 0.0, b.index(0.0)[0])[0]
        for n in expected:
            np.testing.assert_array_equal(got[n], expected[n])

    def test_damaged_file_is_not_retried_as_lock_clash(self):
        """Only the lock clash is retried: a file that is not HDF5 keeps its message."""
        p = str(self.tmp / "junk.h5")
        Path(p).write_bytes(b"not an hdf5 file at all" * 100)
        with self.assertRaisesRegex(PhareseeIOError, "truncated or not an HDF5 file"):
            H5Backend(p)

    def test_truncated_after_open(self):
        p = self.write("ta.h5", dim=2, patches={0.0: {0: tiling(2, [16, 8], [4, 4])}},
                       domain_cells=[16, 8])
        b = H5Backend(p)
        li = b.index(0.0)[0]
        os.truncate(p, os.path.getsize(p) // 2)
        with self.assertRaisesRegex(PhareseeFileChanged, "shrank"):
            read_level(b, 0.0, li)

    def test_bad_environment_values(self):
        p = self.write("env.h5", dim=1, patches={0.0: {0: [(0, 0, (0,), (7,))]}},
                       domain_cells=[8])
        for name, value in (("PHARESEE_IO_THREADS", "abc"), ("PHARESEE_RAW_IO", "maybe"),
                            ("PHARESEE_IO_THREADS", "2.5")):
            with mock.patch.dict(os.environ, {name: value}), \
                    self.assertRaisesRegex(ValueError, name):
                H5Backend(p)
        with mock.patch.dict(os.environ, {"PHARESEE_IO_THREADS": "3", "PHARESEE_RAW_IO": "no"}):
            b = H5Backend(p)
        self.assertEqual(b._nthreads, 3)
        self.assertFalse(b._raw_io)
        with mock.patch.dict(os.environ, {"PHARESEE_CACHE": "perhaps"}), \
                self.assertRaisesRegex(ValueError, "PHARESEE_CACHE"):
            b.index(0.0)


class TestFloat32Time(SyntheticBase):
    def test_float32_request(self):
        path = self.write(dim=1, patches={t: {0: [(0, 0, (0,), (7,))]} for t in (0.0, 0.2, 0.4)},
                          domain_cells=[8])
        b = H5Backend(path)
        self.assertEqual(b.time_key(np.float32(0.2)), "0.2000000000")
        self.assertEqual(b.time_key(np.array(0.4, dtype=np.float32)), "0.4000000000")
        self.assertEqual(len(b.index(np.float32(0.2))), 1)
        with self.assertRaises(TimeNotFoundError):  # the float64 value of float32(0.2)
            b.time_key(float(np.float32(0.2)))
        with self.assertRaises(TimeNotFoundError):
            b.time_key(np.float32(0.2) + np.float32(1e-6))
        self.assertFalse(b.reads_whole_level(0.2, 0, [0]))


class TestDerived(SyntheticBase):
    def test_linear_B_curl_and_div_exact(self):
        for dim in (1, 2, 3):
            with self.subTest(dim=dim):
                path = self.write(f"b{dim}.h5", dim=dim, patches={0.0: two_level_patches(dim)},
                                  domain_cells=[16] * dim, nan_outer_ghost=True)
                b = H5Backend(path)
                coef = [linear_coeffs(d)[1] for d in DATASETS["EM_B"]]  # coef[c][a] = dB_c/dx_a

                def d(c, a, coef=coef, dim=dim):
                    return coef[c][a] if a < dim else 0.0

                J = [d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)]
                div = sum(d(a, a) for a in range(dim))
                for li in b.index(0.0):
                    nn, dst, total = layout_offsets(li.lower, li.upper)
                    ids = np.arange(len(nn))
                    out = {n: np.full(total, np.nan, np.float32) for n in "xyz"}
                    b.read_derived("J_yee", 0.0, li.level, ids, out, dst)
                    for a, n in enumerate("xyz"):
                        np.testing.assert_array_equal(out[n], np.float32(J[a]))
                    out = {"value": np.full(total, np.nan)}
                    b.read_derived("divB_yee", 0.0, li.level, ids, out, dst)
                    np.testing.assert_array_equal(out["value"], div)

    def test_derived_needs_B(self):
        path = self.write(dim=1, quantity="EM_E", patches={0.0: {0: [(0, 0, (0,), (7,))]}},
                          domain_cells=[8])
        b = H5Backend(path)
        out = {n: np.zeros(9) for n in "xyz"}
        with self.assertRaises(ValueError):
            b.read_derived("J_yee", 0.0, 0, [0], out, [0])
        with self.assertRaises(ValueError):
            H5Backend(self.write("b.h5", dim=1, patches={0.0: {0: [(0, 0, (0,), (7,))]}},
                                 domain_cells=[8])).read_derived("curl", 0.0, 0, [0], out, [0])


# ----------------------------------------------------------------------------- real data


def vtk_level(path, time, level):
    """Naive vtkhdf reader: per patch (lower, upper, (n_nodes, ncomp) values), file order.

    Returns None when the level does not exist at that time.
    """
    import h5py

    with h5py.File(path, "r") as f:
        v = f["VTKHDF"]
        values = v["Steps/Values"][()]
        s = int(np.argmin(np.abs(values - time)))
        assert abs(values[s] - time) < 1e-6
        if f"Level{level}" not in v:
            return None
        st = v[f"Steps/Level{level}"]
        nb = int(st["NumberOfAMRBox"][s])
        if nb == 0:
            return None
        bo = int(st["AMRBoxOffset"][s])
        po = int(st["PointDataOffset/data"][s])
        dim = int(f.attrs["dimension"])
        boxes = v[f"Level{level}/AMRBox"][bo : bo + nb].astype(np.int64)
        data = v[f"Level{level}/PointData/data"]
        copies = {1: 4, 2: 2, 3: 1}[dim]
        out = []
        off = po
        for box in boxes:
            lo, up = box[0::2][:dim], box[1::2][:dim]
            n = int(np.prod(up - lo + 2))
            out.append((lo, up, data[off : off + n]))
            off += n * copies
        return out


class TestRealData(unittest.TestCase):
    ROOT, CASES = SAMPLES, SAMPLE_CASES  # test_mpi3.py reruns these tests on the 3-rank run
    MIN_COMPARED = 50  # (quantity, time, level) value checks in test_h5_equals_vtkhdf

    @classmethod
    def setUpClass(cls):
        require(cls.ROOT, cls.CASES)
        cls.runs = [Path(cls.ROOT) / c for c in cls.CASES]

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._old = os.environ.get("PHARESEE_CACHE_DIR")
        os.environ["PHARESEE_CACHE_DIR"] = self._tmp.name

    def tearDown(self):
        if self._old is None:
            os.environ.pop("PHARESEE_CACHE_DIR", None)
        else:
            os.environ["PHARESEE_CACHE_DIR"] = self._old
        self._tmp.cleanup()

    def test_h5_equals_vtkhdf(self):
        """Stored quantities, every time and level: same boxes in the same order,
        values equal to float32 rounding.

        The 3-rank sample (sample_data_mpi3) wrote its h5 and vtkhdf files in two
        separate runs, and runs on 3+ MPI ranks are not bitwise reproducible
        (issue #1526): there, values are compared at t=0 only (initial
        condition), boxes at every time.
        """
        compared = 0
        for run in self.runs:
            values_up_to = 0.0 if "mpi3" in str(run) else np.inf
            for h5 in sorted((run / "h5").glob("*.h5")):
                vtk = run / "vtk" / (h5.stem + ".vtkhdf")
                if h5.stem == "tags" or not vtk.exists():
                    continue
                b = H5Backend(str(h5))
                for t in b.times():
                    for li in b.index(t):
                        ref = vtk_level(str(vtk), t, li.level)
                        if ref is None:  # vtkhdf writer does not write this quantity
                            continue
                        with self.subTest(run=str(run), q=h5.stem, t=t, level=li.level):
                            np.testing.assert_array_equal(li.lower, [r[0] for r in ref])
                            np.testing.assert_array_equal(li.upper, [r[1] for r in ref])
                            if t > values_up_to:
                                continue
                            out, _, _ = read_level(b, t, li)
                            for c, name in enumerate(b.meta().names):
                                vals = np.concatenate([r[2][:, c] for r in ref])
                                scale = max(float(np.abs(vals).max()), 1e-30)
                                err = np.abs(out[name].astype(np.float64) - vals).max()
                                self.assertLessEqual(err, 2e-7 * scale, f"{name}")
                            compared += 1
        self.assertGreater(compared, self.MIN_COMPARED)

    def _master_pairs(self, hier, t, li):
        """Master's patches of level li.level keyed by lower corner."""
        lvl = hier.levels(t)[li.level]
        return {tuple(p.box.lower.tolist()): p for p in lvl.patches}

    def test_J_yee_equals_master(self):
        from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from
        from pyphare.pharesee.run import Run
        from pyphare.pharesee.run.utils import _compute_current

        for run in self.runs:
            b = H5Backend(str(run / "h5" / "EM_B.h5"))
            dim = b.meta().dim
            if dim == 3:
                continue  # master has no 3D current
            r = Run(str(run / "h5"), legacy=True)  # master's reader
            for t in b.times():
                if dim == 1:  # master's GetJ fails in 1D (asks for Jx): use its raw Yee J
                    J = compute_hier_from(_compute_current, r.GetB(t, all_primal=False))
                    names = {"y": "Jy", "z": "Jz"}
                else:
                    J = r.GetJ(t)
                    names = {"x": "x", "y": "y", "z": "z"}
                for li in b.index(t):
                    with self.subTest(run=run.name, t=t, level=li.level):
                        nn, dst, total = layout_offsets(li.lower, li.upper)
                        out = {n: np.full(total, np.nan, np.float32) for n in "xyz"}
                        b.read_derived("J_yee", t, li.level, np.arange(len(nn)), out, dst)
                        master = self._master_pairs(J, t, li)
                        for p in range(len(nn)):
                            mp = master[tuple(li.lower[p].tolist())]
                            if dim == 1:
                                self.assertTrue((patch_views(out["x"], nn, dst)[p] == 0).all())
                            for n, mn in names.items():
                                pd = mp.patch_datas[mn]
                                g = int(pd.ghosts_nbr[0])
                                ref = pd.dataset[(slice(g, -g),) * dim]
                                mine = patch_views(out[n], nn, dst)[p]
                                self.assertEqual(mine.shape, ref.shape)
                                self.assertTrue(ulp_close(mine, ref).all(), f"J{n} patch {p}")

    def test_divB_yee_is_centred_node_divergence(self):
        for run in self.runs:
            b = H5Backend(str(run / "h5" / "EM_B.h5"))
            dim = b.meta().dim
            for t in b.times():
                for li in b.index(t):
                    with self.subTest(run=run.name, t=t, level=li.level):
                        dl = b.meta().dl0 / 2**li.level
                        Bn, nn, dst = read_level(b, t, li, dtype=np.float64)
                        out = {"value": np.full(Bn["x"].size, np.nan)}
                        b.read_derived("divB_yee", t, li.level, np.arange(len(nn)), out, dst)
                        scale = max(np.abs(Bn[n]).max() for n in "xyz") / dl.min()
                        for p in range(len(nn)):
                            if (nn[p] < 3).any():
                                continue
                            inner = (slice(1, -1),) * dim
                            div = np.zeros(tuple(nn[p] - 2))
                            for a, n in enumerate("xyz"[:dim]):
                                v = patch_views(Bn[n], nn, dst)[p]
                                hi = [slice(1, -1)] * dim
                                lo = [slice(1, -1)] * dim
                                hi[a], lo[a] = slice(2, None), slice(0, -2)
                                div += (v[tuple(hi)] - v[tuple(lo)]) / (2 * dl[a])
                            got = patch_views(out["value"], nn, dst)[p][inner]
                            np.testing.assert_allclose(got, div, rtol=0, atol=1e-6 * scale)


if __name__ == "__main__":
    unittest.main()
