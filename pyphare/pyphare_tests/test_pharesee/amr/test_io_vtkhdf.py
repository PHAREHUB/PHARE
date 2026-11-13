"""Tests of pharesee.amr.io.vtkhdf.VtkhdfBackend (DESIGN sections 4.1, 5, 5.1, 11.2)."""

import glob
import os
import tempfile
import unittest
import warnings
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.errors import (
    FormatError,
    PhareseeFileChanged,
    PhareseeIOError,
    TimeNotFoundError,
)
from pyphare.pharesee.amr.io import btree, open_backend, rawio
from pyphare.pharesee.amr.io import vtkhdf as vtk
from pyphare.pharesee.amr.io.vtkhdf import VtkhdfBackend
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES as SAMPLE_ROOT
from pyphare_tests.test_pharesee.amr.sample_data import require

try:
    from .vtkhdf_files import (
        append_vtkhdf,
        default_fields,
        naive_read,
        patch_values,
        tiling,
        write_vtkhdf,
    )
except ImportError:  # run as a script
    from vtkhdf_files import (
        append_vtkhdf,
        default_fields,
        naive_read,
        patch_values,
        tiling,
        write_vtkhdf,
    )

SAMPLE_CASES = ("1d", "2d_static", "2d_tagged", "3d")


def hierarchy(dim, nsteps=2):
    """Boxes per step per level: L0 tiling with uneven patches + L1 with 2-node patches."""
    if dim == 1:
        l0 = tiling([60], [13])
        l1 = [([20], [20]), ([21], [40]), ([60], [61]), ([80], [99])]
    elif dim == 2:
        l0 = tiling([30, 20], [11, 7])
        l1 = [([10, 10], [10, 10]), ([11, 10], [24, 17]), ([30, 4], [41, 9])]
    else:
        l0 = tiling([10, 8, 6], [5, 4, 3])
        l1 = [([4, 4, 4], [4, 4, 4]), ([5, 4, 4], [11, 9, 7])]
    steps = [[l0, l1] for _ in range(nsteps)]
    steps[-1] = [l0, l1[:-1]]  # patch count changes with time
    return steps


def read(backend, time, level, ids=None, dst=None, dtype=None, fill=np.nan, names=None):
    """read_into with buffers pre-filled with `fill` for components `names` (None: all);
    returns (out, ids, dst, lix)."""
    lix = next(x for x in backend.index(time) if x.level == level)
    npts = lix.payload.npts
    ids = np.arange(len(npts)) if ids is None else np.asarray(ids, np.int64)
    if dst is None:
        dst = np.concatenate([[0], np.cumsum(npts[ids])[:-1]]).astype(np.int64)
    size = int((dst + npts[ids]).max()) + 5 if len(ids) else 5
    dtype = backend.meta().dtype if dtype is None else dtype
    names = backend.meta().names if names is None else names
    out = {n: np.full(size, fill, dtype=dtype) for n in names}
    backend.read_into(time, level, ids, out, dst)
    return out, ids, dst, lix


def assert_exact(test, backend, out, ids, dst, lix, fields, t, level, dl0, dtype=None):
    """Each requested patch holds the analytic values, nothing else was written."""
    dtype = backend.meta().dtype if dtype is None else dtype
    file_dtype = backend.meta().dtype
    names = backend.meta().names
    written = np.zeros(len(out[names[0]]), bool)
    for i, d in zip(ids.tolist(), dst.tolist()):
        n = int(lix.payload.npts[i])
        ref = patch_values(fields, t, level, lix.lower[i], lix.upper[i], dl0, file_dtype)
        for c, name in enumerate(names):
            got = out[name][d : d + n]
            np.testing.assert_array_equal(got, ref[:, c].astype(dtype), err_msg=f"patch {i}")
        written[d : d + n] = True
    for name in names:
        test.assertTrue(np.isnan(out[name][~written]).all(), "values written outside patches")


def _h5py():
    import h5py

    return h5py


class TmpDirTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, steps, **kw):
        path = os.path.join(self.tmp, name)
        write_vtkhdf(path, steps, **kw)
        return path


class IndexTest(TmpDirTest):
    def test_index_equals_written(self):
        for dim in (1, 2, 3):
            for style in ("production", "master"):
                steps = hierarchy(dim, 3)
                steps[1] = [steps[1][0], []]  # level 1 empty at step 1
                times = [0.0, 0.25, 1.5]
                p = self.write(f"{dim}{style}.vtkhdf", steps, times=times, style=style)
                b = VtkhdfBackend(p)
                np.testing.assert_array_equal(b.times(), times)
                for s, t in enumerate(times):
                    ix = b.index(t)
                    want = [lv for lv in range(2) if len(steps[s][lv])]
                    self.assertEqual([x.level for x in ix], want)
                    for lix in ix:
                        boxes = np.asarray(steps[s][lix.level], np.int64).reshape(-1, 2, dim)
                        np.testing.assert_array_equal(lix.lower, boxes[:, 0])
                        np.testing.assert_array_equal(lix.upper, boxes[:, 1])
                        self.assertEqual(lix.lower.dtype, np.int64)
                        npts = np.prod(boxes[:, 1] - boxes[:, 0] + 2, axis=1)
                        np.testing.assert_array_equal(lix.payload.npts, npts)
                    self.assertIs(b.index(t), ix)  # cached

    def test_meta(self):
        steps = hierarchy(2)
        p = self.write(
            "ions_pop_alpha_flux.vtkhdf", steps, style="master", dl0=[0.4, 0.5],
            domain_cells=[100, 80], periodic=[True, False],
            attrs={"pop_mass": 4.0, "pop_charge": 2.0},
        )
        m = VtkhdfBackend(p).meta()
        self.assertEqual(m.format, "vtkhdf")
        self.assertEqual(m.quantity, "ions_pop_alpha_flux")
        self.assertEqual(m.dim, 2)
        np.testing.assert_array_equal(m.dl0, [0.4, 0.5])  # float32 0.4 rounded back
        self.assertEqual(m.dl0.dtype, np.float64)
        np.testing.assert_array_equal(m.domain_cells0, [100, 80])
        np.testing.assert_array_equal(m.periodic, [True, False])
        self.assertEqual((m.interp_order, m.model), (1, "hybrid"))
        self.assertEqual(m.names, ("x", "y", "z"))
        self.assertEqual(m.dtype, np.float32)
        self.assertTrue(m.colocated)
        self.assertEqual(m.attrs, {"pop_mass": 4.0, "pop_charge": 2.0})
        self.assertIn("(87381, 3)", m.writer_signature)
        self.assertIn("S1 off", m.writer_signature)
        p = self.write("prod.vtkhdf", steps)
        self.assertIn("S1) on", VtkhdfBackend(p).meta().writer_signature)

    def test_component_names(self):
        for ncomp, names in [
            (1, ("value",)), (3, ("x", "y", "z")),
            (6, ("xx", "xy", "xz", "yy", "yz", "zz")),
            (9, ("xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz")),
        ]:
            fields = default_fields(ncomp)
            p = self.write(f"c{ncomp}.vtkhdf", hierarchy(2), fields=fields)
            b = VtkhdfBackend(p)
            self.assertEqual(b.meta().names, names)
            out, ids, dst, lix = read(b, 0.0, 1)
            assert_exact(self, b, out, ids, dst, lix, fields, 0.0, 1, [0.4, 0.4])

    def test_mhd_file_without_interp_order(self):
        p = self.write("mhd.vtkhdf", hierarchy(2), interp_order=None)
        m = VtkhdfBackend(p).meta()
        self.assertEqual((m.interp_order, m.model), (None, "mhd"))

    def test_old_boundary_conditions_attribute(self):
        p = self.write("old.vtkhdf", hierarchy(2), boundary_conditions=True,
                       periodic=[False, True])
        np.testing.assert_array_equal(VtkhdfBackend(p).meta().periodic, [False, True])


class TimeTest(TmpDirTest):
    def test_time_lookup_and_tolerance(self):
        times = [0.0, 0.007, 0.014, 0.02, 0.5]
        p = self.write("t.vtkhdf", [hierarchy(1, 1)[0]] * 5, times=times, style="master")
        b = VtkhdfBackend(p)
        # float32 times come back as the shortest decimal
        self.assertEqual(b.times().tolist(), times)
        self.assertEqual(b.index(0.007)[0].payload.step, 1)
        self.assertEqual(b.index(float(np.float32(0.007)))[0].payload.step, 1)
        self.assertEqual(b.index(0.007 * (1 + 1e-8))[0].payload.step, 1)
        with self.assertRaises(TimeNotFoundError) as cm:
            b.index(0.0071)
        e = cm.exception
        self.assertIsInstance(e, KeyError)
        self.assertIsInstance(e, ValueError)
        msg = str(e)
        for near in ("0.007", "0.014", "0.0"):
            self.assertIn(near, msg)
        self.assertNotIn("0.5", msg)
        with self.assertRaises(KeyError):
            b.read_into(3.0, 0, [0], {}, [0])
        with self.assertRaises(ValueError):
            b.index(-1.0)


class ReadTest(TmpDirTest):
    def check_all_strategies(self, dim, style, **kw):
        steps = hierarchy(dim)
        fields = default_fields(3)
        p = self.write(f"r{dim}{style}.vtkhdf", steps, style=style, fields=fields, **kw)
        dl0 = [0.4] * dim
        for strategy in ("auto", "s1", "s2", "s3"):
            b = VtkhdfBackend(p, strategy=strategy)
            for t in b.times():
                for lix in b.index(t):
                    lv = lix.level
                    n = len(lix.lower)
                    # all patches
                    out, ids, dst, _ = read(b, t, lv)
                    assert_exact(self, b, out, ids, dst, lix, fields, t, lv, dl0)
                    # subset
                    out, ids, dst, _ = read(b, t, lv, ids=np.arange(0, n, 3))
                    assert_exact(self, b, out, ids, dst, lix, fields, t, lv, dl0)
                    # scattered single patches, destinations not in file order, with gaps
                    rng = np.random.default_rng(n)
                    ids = np.sort(rng.choice(n, size=max(1, n // 2), replace=False))
                    npts = lix.payload.npts[ids]
                    slots = rng.permutation(len(ids))
                    sizes = npts[slots] + 3
                    starts = np.concatenate([[1], 1 + np.cumsum(sizes)[:-1]])
                    dst = np.empty(len(ids), np.int64)
                    dst[slots] = starts
                    out, ids, dst, _ = read(b, t, lv, ids=ids, dst=dst)
                    assert_exact(self, b, out, ids, dst, lix, fields, t, lv, dl0)
                    for i in (0, n - 1):
                        out, ids, dst, _ = read(b, t, lv, ids=[i])
                        assert_exact(self, b, out, ids, dst, lix, fields, t, lv, dl0)

    def test_production_1d_2d_3d(self):
        for dim in (1, 2, 3):
            self.check_all_strategies(dim, "production")

    def test_master_1d_2d_3d(self):
        for dim in (1, 2, 3):
            self.check_all_strategies(dim, "master")

    def test_old_variants(self):
        self.check_all_strategies(2, "production", offsets_dtype=np.int32, amrbox_cols=4)

    def test_output_dtype_cast(self):
        p = self.write("cast.vtkhdf", hierarchy(2))
        b = VtkhdfBackend(p)
        out, ids, dst, lix = read(b, 0.0, 0, dtype=np.float32)
        assert_exact(self, b, out, ids, dst, lix, default_fields(), 0.0, 0, [0.4, 0.4],
                     dtype=np.float32)

    def test_bad_arguments(self):
        p = self.write("args.vtkhdf", hierarchy(2))
        b = VtkhdfBackend(p)
        lix = b.index(0.0)[0]
        n = int(lix.payload.npts[:3].sum())
        out = {k: np.zeros(n) for k in "xyz"}
        with self.assertRaises(ValueError):
            b.read_into(0.0, 0, [1, 0], out, [0, 10])
        with self.assertRaises(ValueError):
            b.read_into(0.0, 0, [0, 0], out, [0, 10])
        with self.assertRaises(ValueError):
            b.read_into(0.0, 0, [0], {}, [0])  # no component
        with self.assertRaises(ValueError):
            b.read_into(0.0, 0, [0], {"x": out["x"], "w": out["y"]}, [0])  # unknown name
        with self.assertRaises(ValueError):
            b.read_into(0.0, 0, [0], out, [n])  # does not fit
        with self.assertRaises(ValueError):
            b.read_into(0.0, 7, [0], out, [0])  # no such level
        b.read_into(0.0, 0, [], out, [])  # empty selection: nothing to do


class StrategyTest(TmpDirTest):
    """S1 == S2 == S3, choice rule, fallbacks."""

    @classmethod
    def setUpClass(cls):
        # production file, 2 steps, ~2.1 M rows at level 0 in the second (~60 MB):
        # a chunk B-tree with 3 levels; read with small windows to exercise window cuts
        cls._big = tempfile.TemporaryDirectory()
        cls.big = os.path.join(cls._big.name, "big.vtkhdf")
        l0 = tiling([1200, 800], [30, 20])
        l1 = tiling([1800, 1200], [32, 24], origin=[400, 300], stop=[1400, 900])
        write_vtkhdf(cls.big, [[l0[:50], l1[:10]], [l0, l1]], times=[0.0, 1.0])

    @classmethod
    def tearDownClass(cls):
        cls._big.cleanup()

    def read_bits(self, b, t, level, ids=None, names=None):
        out = read(b, t, level, ids=ids, names=names)[0]
        return {k: v.view(np.uint64).copy() for k, v in out.items()}, dict(b.last_read)

    def test_s1_s2_s3_bit_identical(self):
        for window in (vtk.WINDOW_BYTES, 3 << 20, 40 << 10):
            for level in (0, 1):
                ref, info = self.read_bits(VtkhdfBackend(self.big, strategy="s3"), 1.0, level)
                self.assertEqual(info["strategy"], "S3")
                b1 = VtkhdfBackend(self.big, strategy="s1", window_bytes=window)
                got, info = self.read_bits(b1, 1.0, level)
                self.assertEqual(info["strategy"], "S1", info)
                if window == 40 << 10 and level == 0:
                    self.assertGreater(info["windows"], 10)
                    self.assertGreater(info["strays"], 0)  # forced cuts split chunks
                for k in ref:
                    np.testing.assert_array_equal(got[k], ref[k])
                got2, info = self.read_bits(VtkhdfBackend(self.big, strategy="s2"), 1.0, level)
                for k in ref:
                    np.testing.assert_array_equal(got2[k], ref[k])
        # analytic check of S1 on the big file
        b = VtkhdfBackend(self.big, strategy="s1", window_bytes=1 << 20)
        out, ids, dst, lix = read(b, 1.0, 0)
        assert_exact(self, b, out, ids, dst, lix, default_fields(), 1.0, 0, [0.4, 0.4])

    def test_component_subsets_bit_identical(self):
        """Reading some columns gives the bits of the full read, per strategy, and only
        those columns are asked for."""
        n = len(VtkhdfBackend(self.big).index(1.0)[0].lower)
        cases = [(None, 0), (np.r_[np.arange(5, 200, 7), np.arange(n - 40, n)], 0),
                 (np.arange(0, 6), 0), (None, 1)]
        for ids, level in cases:
            ref, _ = self.read_bits(VtkhdfBackend(self.big, strategy="s3"), 1.0, level, ids)
            for names in (("z",), ("x", "z"), ("y",), ("x", "y", "z")):
                for strategy, window in (("s1", vtk.WINDOW_BYTES), ("s1", 40 << 10),
                                         ("s2", None), ("s3", None), ("s3", 40 << 10)):
                    b = VtkhdfBackend(self.big, strategy=strategy, window_bytes=window)
                    got, info = self.read_bits(b, 1.0, level, ids, names=names)
                    what = (ids is not None, level, names, strategy, window, info)
                    self.assertEqual(info["strategy"], strategy.upper(), what)
                    self.assertIsNone(info["fallback"], what)
                    cols = tuple(sorted("xyz".index(c) for c in names))
                    self.assertEqual(info["columns"], cols, what)
                    self.assertEqual(sorted(got), sorted(names))
                    for k in names:
                        np.testing.assert_array_equal(got[k], ref[k], err_msg=str(what))

    def test_component_subset_s1_skips_other_columns(self):
        b = VtkhdfBackend(self.big, strategy="s1", window_bytes=40 << 10)
        read(b, 1.0, 0)
        full = dict(b.last_read)
        read(b, 1.0, 0, names=("z",))
        one = dict(b.last_read)
        self.assertLess(one["strays"], full["strays"])  # unwanted strays not fetched
        with mock.patch.object(vtk._S1Reader, "_decode_cluster_cols",
                               side_effect=AssertionError("subset path used")):
            read(b, 1.0, 0)  # every column: the all-columns decoder only

    def test_nbytes_per_component(self):
        b = VtkhdfBackend(self.big)
        ids = np.arange(10)
        full = b.nbytes(1.0, 0, ids)
        self.assertEqual(b.nbytes(1.0, 0, ids, names=("z",)) * 3, full)
        self.assertEqual(b.nbytes(1.0, 0, ids, names=["x", "y"]) * 3, 2 * full)
        self.assertEqual(b.nbytes(1.0, 0, ids, names=("x", "y", "z")), full)
        with self.assertRaises(ValueError):
            b.nbytes(1.0, 0, ids, names=("w",))

    def test_choice_rule(self):
        b = VtkhdfBackend(self.big)
        lix = b.index(1.0)[0]
        npts = lix.payload.npts
        self.assertEqual(b.choose_strategy(0, lix, npts), "S1")
        self.assertEqual(b.choose_strategy(0, lix, npts[: len(npts) // 20]), "S2")
        self.assertEqual(b.choose_strategy(0, lix, npts[: len(npts) // 10]), "S1")
        with mock.patch.dict(os.environ, {"PHARESEE_RAW_IO": "0"}):
            self.assertEqual(b.choose_strategy(0, lix, npts), "S3")
            read(b, 1.0, 1)
            self.assertEqual(b.last_read["strategy"], "S3")
        read(b, 1.0, 1, ids=[0, 1])
        self.assertEqual(b.last_read["strategy"], "S2")
        p = self.write("m.vtkhdf", hierarchy(2), style="master")
        bm = VtkhdfBackend(p)
        lix = bm.index(0.0)[0]
        self.assertEqual(bm.choose_strategy(0, lix, lix.payload.npts), "S3")

    def check_fallback(self, path, reason):
        b = VtkhdfBackend(path, strategy="s1")
        for t in b.times():
            for lix in b.index(t):
                out, ids, dst, _ = read(b, t, lix.level)
                self.assertEqual(b.last_read["strategy"], "S3")
                self.assertRegex(b.last_read["fallback"], reason)
                assert_exact(self, b, out, ids, dst, lix, default_fields(), t, lix.level,
                             [0.4, 0.4])

    def test_fallback_on_unsupported_layouts(self):
        """master chunks, compression, libver latest: the chunk walker refuses, S3 reads."""
        for name, kw, reason in (("m.vtkhdf", dict(style="master"), "chunk shape"),
                                 ("z.vtkhdf", dict(compression="gzip"), "filter"),
                                 ("l.vtkhdf", dict(libver="latest"), "version")):
            with self.subTest(reason=reason):
                self.check_fallback(self.write(name, hierarchy(2), **kw), reason)

    def test_fallback_on_missing_chunk(self):
        real = btree.ChunkBTree.parse_leaves

        def drop_one(self_, buf, rel):
            rows, cols, addrs = real(self_, buf, rel)
            return rows[1:], cols[1:], addrs[1:]

        b = VtkhdfBackend(self.big, strategy="s1", window_bytes=4 << 20)
        with mock.patch.object(btree.ChunkBTree, "parse_leaves", drop_one):
            out, ids, dst, lix = read(b, 1.0, 1)
        self.assertEqual(b.last_read["strategy"], "S3")
        self.assertIn("inconsistent", b.last_read["fallback"])
        assert_exact(self, b, out, ids, dst, lix, default_fields(), 1.0, 1, [0.4, 0.4])

    def test_fallback_when_file_changes_during_read(self):
        counter = iter(range(1000))
        b = VtkhdfBackend(self.big, strategy="s1")

        def moving_mtime(fd):  # same size, new mtime at every call
            return os.fstat(fd).st_size, next(counter)

        with mock.patch.object(rawio, "stat_signature", moving_mtime):
            out, ids, dst, lix = read(b, 1.0, 1)
        self.assertEqual(b.last_read["strategy"], "S3")
        self.assertIn("modified", b.last_read["fallback"])
        assert_exact(self, b, out, ids, dst, lix, default_fields(), 1.0, 1, [0.4, 0.4])


class RobustnessTest(TmpDirTest):
    """DESIGN 5.1 odd and broken files."""

    def steps(self, n=3, dim=2):
        return [hierarchy(dim, 1)[0]] * n

    def test_incomplete_last_values(self):
        p = self.write("a.vtkhdf", self.steps(), times=[0, 1, 2], incomplete_last="values")
        with self.assertWarnsRegex(RuntimeWarning, "last dump incomplete, ignored"):
            b = VtkhdfBackend(p)
        self.assertEqual(b.times().tolist(), [0, 1, 2])

    def test_incomplete_last_level0(self):
        p = self.write("b.vtkhdf", self.steps(), times=[0, 1, 2], incomplete_last="level0")
        with self.assertWarnsRegex(RuntimeWarning, "last dump incomplete"):
            b = VtkhdfBackend(p)
        self.assertEqual(b.times().tolist(), [0, 1, 2])
        out, ids, dst, lix = read(b, 2.0, 0)
        assert_exact(self, b, out, ids, dst, lix, default_fields(), 2.0, 0, [0.4, 0.4])

    def test_emergency_dump(self):
        p = self.write("c.vtkhdf", self.steps(), times=[0, 1, 2], emergency_dump=True)
        with _h5py().File(p) as f:
            self.assertEqual(len(f["VTKHDF/Steps/Values"]), 4)  # 2 levels: one extra time
        with self.assertWarnsRegex(RuntimeWarning, "incomplete"):
            b = VtkhdfBackend(p)
        self.assertEqual(b.times().tolist(), [0, 1, 2])
        for lv in (0, 1):
            out, ids, dst, lix = read(b, 2.0, lv)
            assert_exact(self, b, out, ids, dst, lix, default_fields(), 2.0, lv, [0.4, 0.4])

    def test_normal_file_does_not_warn(self):
        p = self.write("ok.vtkhdf", self.steps())
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            VtkhdfBackend(p)

    def test_per_level_dumps_are_refused(self):
        p = self.write("d.vtkhdf", self.steps(), per_level_misalign=True)
        with self.assertRaisesRegex(FormatError, "per-level"):
            VtkhdfBackend(p)

    def test_restart_duplicate_times_last_copy_wins(self):
        base = hierarchy(2, 1)[0]
        other = [base[0][:3], base[1][:1]]
        steps = [base, base, base, other, base]
        p = self.write("e.vtkhdf", steps, times=[0, 0.5, 1.0, 1.0, 1.5])
        b = VtkhdfBackend(p)
        self.assertEqual(b.times().tolist(), [0, 0.5, 1.0, 1.5])
        ix = b.index(1.0)
        self.assertEqual(ix[0].payload.step, 3)
        self.assertEqual(len(ix[0].lower), 3)
        out, ids, dst, _ = read(b, 1.0, 0)
        ref = naive_read(p, 3, 0)
        for i, d in zip(ids, dst):
            n = len(ref[i][2])
            np.testing.assert_array_equal(out["x"][d : d + n], ref[i][2][:, 0])
        # a restart from an earlier time: repeated times, the later copy wins
        p = self.write("f.vtkhdf", [base] * 6, times=[0, 1, 2, 1, 2, 3])
        b = VtkhdfBackend(p)
        self.assertEqual(b.times().tolist(), [0, 1, 2, 3])
        self.assertEqual([b.index(t)[0].payload.step for t in b.times()], [0, 3, 4, 5])

    def test_int32_saturation(self):
        p = self.write("g.vtkhdf", self.steps(4), saturate_int32=True)
        b = VtkhdfBackend(p)
        b.index(0.0)
        b.index(0.5)
        for t in (1.0, 1.5):  # the step whose next offset saturates, and the last one
            with self.assertRaisesRegex(FormatError, "2\\*\\*31-1"):
                b.index(t)

    def test_offsets_inconsistent(self):
        p = self.write("h.vtkhdf", self.steps())
        with _h5py().File(p, "a") as f:
            f["VTKHDF/Steps/Level0/PointDataOffset/data"][1] += 2
        b = VtkhdfBackend(p)
        with self.assertRaisesRegex(FormatError, "PointDataOffset"):
            b.index(0.0)
        p = self.write("i.vtkhdf", self.steps())
        with _h5py().File(p, "a") as f:
            f["VTKHDF/Steps/Level1/AMRBoxOffset"][2] -= 1
        with self.assertRaisesRegex(FormatError, "AMRBoxOffset"):
            VtkhdfBackend(p).index(0.5)

    def test_truncated_data(self):
        p = self.write("j.vtkhdf", self.steps())
        with _h5py().File(p, "a") as f:
            d = f["VTKHDF/Level0/PointData/data"]
            d.resize((d.shape[0] - 10, 3))
        with self.assertRaisesRegex(FormatError, "beyond"):
            VtkhdfBackend(p).index(1.0)

    def test_format_checks(self):
        p = self.write("m.vtkhdf", self.steps())
        with _h5py().File(p, "a") as f:
            f["VTKHDF"].attrs["Type"] = "ImageData"
        with self.assertRaisesRegex(FormatError, "Type"):
            VtkhdfBackend(p)
        p = self.write("n.vtkhdf", self.steps(), fields=default_fields(2))
        with self.assertRaisesRegex(FormatError, "PointData/data"):
            VtkhdfBackend(p)
        p = self.write("o.vtkhdf", self.steps())
        with _h5py().File(p, "a") as f:
            f.attrs["domain_box"] = np.array([10, 10, 10], np.int32)
        with self.assertRaisesRegex(FormatError, "domain_box"):
            VtkhdfBackend(p)
        p = self.write("q.vtkhdf", self.steps())
        with _h5py().File(p, "a") as f:
            del f["VTKHDF/Steps/Level1/NumberOfAMRBox"]
        with self.assertRaisesRegex(FormatError, "NumberOfAMRBox"):
            VtkhdfBackend(p)


class IdentityTest(TmpDirTest):
    def test_append_accepted_replace_refused(self):
        steps = hierarchy(2)
        p = self.write("id.vtkhdf", steps)
        b = VtkhdfBackend(p)
        b.index(0.5)
        with open(p, "ab") as f:  # a running simulation appends
            f.write(b"\0" * 100_000)
        out, ids, dst, lix = read(b, 0.5, 1)
        assert_exact(self, b, out, ids, dst, lix, default_fields(), 0.5, 1, [0.4, 0.4])
        spec = b.spec()
        new = os.path.join(self.tmp, "new.vtkhdf")
        write_vtkhdf(new, steps)
        os.replace(new, p)
        with self.assertRaises(PhareseeFileChanged):
            read(b, 0.5, 1)
        with self.assertRaises(PhareseeFileChanged):
            open_backend(spec[1], **spec[2])
        self.assertTrue(issubclass(PhareseeFileChanged, PhareseeIOError))

    def test_read_error_is_wrapped(self):
        p = self.write("err.vtkhdf", hierarchy(2))
        b = VtkhdfBackend(p, strategy="s3")
        b.index(0.0)
        with (
            mock.patch.object(VtkhdfBackend, "_read_s3", side_effect=OSError("disk")),
            self.assertRaises(PhareseeIOError) as cm,
        ):
            read(b, 0.0, 0)
        msg = str(cm.exception)
        for part in ("err.vtkhdf", "time 0.0", "level 0", "patches"):
            self.assertIn(part, msg)


def rewrite_in_place(path, source):
    """Overwrite `path` with the bytes of `source`, keeping the inode (a re-run into the
    same folder: open, truncate, write)."""
    with open(source, "rb") as src, open(path, "r+b") as dst:
        data = src.read()
        dst.truncate(0)
        dst.write(data)


def shifted_fields(ncomp=3):
    """Another quantity on the same layout: different values everywhere."""
    base = default_fields(ncomp)
    return [lambda t, level, *c, f=f: 2.0 * f(t, level, *c) + 7.0 for f in base]


class RewriteInPlaceTest(TmpDirTest):
    """A1: a backend opened before its file is rewritten in place (same inode) never
    returns numbers from the new content; a pure append keeps working."""

    def check_stale_raises(self, b, times):
        for t in times:
            with self.assertRaisesRegex(PhareseeFileChanged, "rewritten.*new Run"):
                b.index(t)
            with self.assertRaises(PhareseeFileChanged):
                read(b, t, 0)
        self.assertEqual(len(b._index_cache), 0)  # nothing cached is served any more

    def test_other_boxes(self):
        p = self.write("q.vtkhdf", hierarchy(2))
        b = VtkhdfBackend(p)
        b.index(0.5)  # a field was created: index cached, nothing read
        other = [[tiling([30, 20], [6, 5]), []], [tiling([30, 20], [10, 4]), []]]
        new = self.write("new.vtkhdf", other)
        rewrite_in_place(p, new)
        self.check_stale_raises(b, (0.5, 0.0))
        fresh = VtkhdfBackend(p)  # a new backend reads the new content
        out, ids, dst, lix = read(fresh, 0.5, 0)
        assert_exact(self, fresh, out, ids, dst, lix, default_fields(), 0.5, 0, [0.4, 0.4])

    def test_same_size_other_data(self):
        p = self.write("q.vtkhdf", hierarchy(2))
        new = self.write("new.vtkhdf", hierarchy(2), fields=shifted_fields())
        self.assertEqual(os.path.getsize(p), os.path.getsize(new))
        for warm in (False, True):
            p = self.write("q.vtkhdf", hierarchy(2))
            b = VtkhdfBackend(p)
            if warm:
                read(b, 0.5, 0)
            else:
                b.index(0.5)
            os.utime(p, ns=(0, 10**18))  # an mtime surely different from the write's
            rewrite_in_place(p, new)
            self.assertEqual(os.stat(p).st_size, os.path.getsize(new))
            self.check_stale_raises(b, (0.5,))
        fresh = VtkhdfBackend(p)
        out, ids, dst, lix = read(fresh, 0.5, 1)
        assert_exact(self, fresh, out, ids, dst, lix, shifted_fields(), 0.5, 1, [0.4, 0.4])

    def test_append_keeps_working(self):
        steps = hierarchy(2)
        for style in ("production", "master"):
            p = self.write(f"a_{style}.vtkhdf", steps, style=style)
            b = VtkhdfBackend(p)
            out, ids, dst, lix = read(b, 0.5, 1)
            before = rawio.file_identity(p)
            append_vtkhdf(p, [steps[0]], times=[1.0])  # a running simulation dumps again
            self.assertTrue(rawio.modified_since(before, rawio.file_identity(p)))
            for t in (0.5, 0.0):  # cached and not yet indexed
                for lv in (0, 1):
                    out, ids, dst, lix = read(b, t, lv)
                    assert_exact(self, b, out, ids, dst, lix, default_fields(), t, lv,
                                 [0.4, 0.4])
            self.assertEqual(b.times().tolist(), [0.0, 0.5])  # times seen at open
            self.assertEqual(tuple(b.spec()[2]["identity"]), tuple(rawio.file_identity(p)))
            fresh = VtkhdfBackend(p)
            self.assertEqual(fresh.times().tolist(), [0.0, 0.5, 1.0])
            out, ids, dst, lix = read(fresh, 1.0, 0)
            assert_exact(self, fresh, out, ids, dst, lix, default_fields(), 1.0, 0, [0.4, 0.4])

    def test_rewrite_after_append_is_detected(self):
        steps = hierarchy(2)
        p = self.write("ra.vtkhdf", steps)
        b = VtkhdfBackend(p)
        b.index(0.0)
        append_vtkhdf(p, [steps[0]], times=[1.0])
        read(b, 0.0, 0)  # accepted, identity refreshed
        bigger = self.write("big.vtkhdf", steps + [steps[0]], fields=shifted_fields(),
                            times=[0.0, 0.5, 1.0, 1.5])
        rewrite_in_place(p, bigger)
        self.check_stale_raises(b, (0.0,))


class DamagedFileTest(TmpDirTest):
    """A9: empty or truncated files give PhareseeIOError / FormatError naming the file."""

    def test_truncated_at_open(self):
        full = self.write("full.vtkhdf", hierarchy(2))
        size = os.path.getsize(full)
        with open(full, "rb") as f:
            data = f.read()
        for n in (0, 8, 100, 2048, size // 100, size // 10, size // 2, 9 * size // 10, size - 1):
            p = os.path.join(self.tmp, f"t{n}.vtkhdf")
            with open(p, "wb") as f:
                f.write(data[:n])
            with self.subTest(n=n):
                with self.assertRaises((PhareseeIOError, FormatError)) as cm:
                    b = VtkhdfBackend(p)
                    for t in b.times():
                        read(b, t, 0)
                self.assertIn(p, str(cm.exception))
                if n == 0:
                    self.assertIn("empty", str(cm.exception))

    def test_file_held_open_by_another_reader(self):
        """Another reader of this session (Ulysse's vtkamr, plain h5py) keeps the file open
        with h5py's default lock settings: the backend still opens and reads it."""
        import h5py

        p = self.write("held.vtkhdf", hierarchy(2))
        expected = read(VtkhdfBackend(p), 0.5, 0)[0]
        with h5py.File(p, "r"):
            if os.environ.get("HDF5_USE_FILE_LOCKING", "").upper() != "FALSE":
                with self.assertRaisesRegex(OSError, "file locking"):  # the clash is real
                    h5py.File(p, "r", locking=False)
            got = read(VtkhdfBackend(p), 0.5, 0)[0]
        for n in expected:
            np.testing.assert_array_equal(got[n], expected[n])

    def test_lock_clash_message(self):
        """When the retry with default lock settings fails too, the message names the
        clash, not a damaged file."""
        import h5py

        from pyphare.pharesee.amr.io import h5open

        p = self.write("clash.vtkhdf", hierarchy(2))
        clash = OSError("Unable to synchronously open file (file locking flag values don't match)")
        with mock.patch.object(h5py, "File", side_effect=clash), \
                self.assertRaisesRegex(PhareseeIOError, "already open elsewhere") as cm:
            h5open.open_h5(p)
        self.assertIn(p, str(cm.exception))
        self.assertNotIn("damaged", str(cm.exception))

    def test_truncated_after_open(self):
        p = self.write("ta.vtkhdf", hierarchy(2))
        b = VtkhdfBackend(p)
        b.index(0.5)
        os.truncate(p, os.path.getsize(p) // 2)
        with self.assertRaisesRegex(PhareseeFileChanged, "shrank"):
            read(b, 0.5, 0)

    def test_bad_environment_values(self):
        p = self.write("env.vtkhdf", hierarchy(2))
        b = VtkhdfBackend(p, strategy="s1")
        for name, value in (("PHARESEE_RAW_IO", "maybe"), ("PHARESEE_IO_THREADS", "abc")):
            with mock.patch.dict(os.environ, {name: value}), \
                    self.assertRaisesRegex(ValueError, name):
                read(b, 0.5, 0)
        with mock.patch.dict(os.environ, {"PHARESEE_RAW_IO": "off"}):
            read(b, 0.5, 0)
            self.assertEqual(b.last_read["strategy"], "S3")


class Float32TimeTest(TmpDirTest):
    """A10: a time passed as np.float32 finds the dump, whatever the stored time dtype."""

    def test_float32_request(self):
        for style in ("production", "master"):  # times stored float64 / float32
            p = self.write(f"f{style}.vtkhdf", [hierarchy(1, 1)[0]] * 3, style=style,
                           times=[0.0, 0.2, 0.4])
            b = VtkhdfBackend(p)
            self.assertEqual(b.index(np.float32(0.2))[0].payload.step, 1)
            self.assertEqual(b.index(np.array(0.4, dtype=np.float32))[0].payload.step, 2)
            self.assertEqual(b.index(np.float16(0.2))[0].payload.step, 1)
            with self.assertRaises(TimeNotFoundError):
                b.index(np.float32(0.2) + np.float32(1e-6))
        p = self.write("f64.vtkhdf", [hierarchy(1, 1)[0]] * 3, times=[0.0, 0.2, 0.4])
        with self.assertRaises(TimeNotFoundError):  # a float64 request keeps the tight rule
            VtkhdfBackend(p).index(float(np.float32(0.2)))


class ReadsWholeLevelTest(TmpDirTest):
    """F2: reads_whole_level is True exactly when the request would be read with S1."""

    def test_rule(self):
        l0 = tiling([300, 200], [10, 10])
        p = self.write("w.vtkhdf", [[l0]], times=[0.0])
        b = VtkhdfBackend(p)
        lix = b.index(0.0)[0]
        n = len(lix.lower)
        for ids in (np.arange(n), np.arange(n // 3), np.arange(5), np.array([0, n - 1])):
            expect = b.choose_strategy(0, lix, lix.payload.npts[ids]) == "S1"
            self.assertEqual(b.reads_whole_level(0.0, 0, ids), expect)
        self.assertTrue(b.reads_whole_level(0.0, 0, np.arange(n // 3)))
        self.assertFalse(b.reads_whole_level(0.0, 0, np.arange(5)))
        self.assertFalse(b.reads_whole_level(0.0, 0, np.zeros(0, np.int64)))
        with mock.patch.dict(os.environ, {"PHARESEE_RAW_IO": "0"}):
            self.assertFalse(b.reads_whole_level(0.0, 0, np.arange(n)))
        pm = self.write("wm.vtkhdf", [[l0]], times=[0.0], style="master")  # S1 refused
        self.assertFalse(VtkhdfBackend(pm).reads_whole_level(0.0, 0, np.arange(n)))


class RealFilesTest(unittest.TestCase):
    """Real PHARE output: every quantity, time and level equals a naive h5py read."""

    ROOT, CASES = SAMPLE_ROOT, SAMPLE_CASES  # test_mpi3.py reruns it on the 3-rank run

    @classmethod
    def setUpClass(cls):
        require(cls.ROOT, cls.CASES)
        cls.samples = sorted(glob.glob(os.path.join(cls.ROOT, "*", "vtk", "*.vtkhdf")))

    def test_samples(self):
        checked = 0
        for path in self.samples:
            with _h5py().File(path, "r") as f:
                has_data = "Level0" in f["VTKHDF"]
                values = f["VTKHDF/Steps/Values"][()]
            if not has_data:
                with self.assertRaisesRegex(FormatError, "no data"):
                    VtkhdfBackend(path)
                continue
            for strategy in ("auto", "s2"):
                b = VtkhdfBackend(path, strategy=strategy)
                self.assertEqual(len(b.times()), len(values))
                for t in b.times():
                    for lix in b.index(t):
                        out, ids, dst, _ = read(b, t, lix.level)
                        ref = naive_read(path, lix.payload.step, lix.level)
                        self.assertEqual(len(ref), len(ids))
                        for i, d in zip(ids.tolist(), dst.tolist()):
                            lower, upper, rows = ref[i]
                            np.testing.assert_array_equal(lix.lower[i], lower)
                            np.testing.assert_array_equal(lix.upper[i], upper)
                            n = len(rows)
                            for c, name in enumerate(b.meta().names):
                                np.testing.assert_array_equal(
                                    out[name][d : d + n], rows[:, c], err_msg=path
                                )
                        checked += 1
        self.assertGreater(checked, 0)


if __name__ == "__main__":
    unittest.main()
