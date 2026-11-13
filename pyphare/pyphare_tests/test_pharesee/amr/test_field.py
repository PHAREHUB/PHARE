#!/usr/bin/env python3
"""Tests of pharesee.amr.field: lazy construction, times, components, load/copy/pickle, repr.

Synthetic data: ``testing.make_amr`` layouts filled with analytic functions and served by
``testing.FakeBackend`` (which counts every read). Real data (skipped when absent):
``build/sample_data/*/{vtk,h5}``, the same runs in both formats.
"""

import copy
import pickle
import unittest
from pathlib import Path

import numpy as np

from pyphare.pharesee.amr.boxes import PhysBox
from pyphare.pharesee.amr.errors import TimeNotFoundError
from pyphare.pharesee.amr.field import (
    AMRField,
    ScalarField,
    TensorField,
    VectorField,
    _missing_bytes,
    make_field,
)
from pyphare.pharesee.amr.geometry import GeometryCache
from pyphare.pharesee.amr.io import open_backend
from pyphare.pharesee.amr.store import Snapshot
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES

SAMPLE = Path(SAMPLES)
RUNS = [r for r in ("1d", "2d_static", "2d_tagged", "3d")
        if (SAMPLE / r / "vtk").is_dir() and (SAMPLE / r / "h5").is_dir()]


# --------------------------------------------------------------------------- helpers


def vec_funcs(t=0.0):
    """Analytic, distinct components (linear in space so probes/derivatives are exact)."""
    def fx(x, y=0.0, z=0.0):
        return 1.0 + 0.5 * x - 0.25 * np.asarray(y) + 0.125 * np.asarray(z) + t + 0 * x

    def fy(x, y=0.0, z=0.0):
        return -2.0 + 0.3 * x + 0.7 * np.asarray(y) + 10 * t + 0 * x

    def fz(x, y=0.0, z=0.0):
        return 3.0 - 0.1 * x + 0.2 * np.asarray(z) + 100 * t + 0 * x

    return {"x": fx, "y": fy, "z": fz}


def coords_of(geom):
    return [geom.node_coords(a) for a in range(geom.dim)]


def fake_vector(dim=2, nlevels=3, seed=1, times=(0.0,), dtype=np.float64, **kw):
    """FakeBackend serving a vector field on a make_amr layout at the given times."""
    hier = make_amr(dim, nlevels, seed, **kw)
    data = {float(t): fill(hier, vec_funcs(t), dtype=dtype, time=t) for t in times}
    return FakeBackend(data, quantity="EM_B"), hier


def fake_scalar(dim=2, nlevels=2, seed=2, times=(0.0,), f=None, names=("value",), **kw):
    hier = make_amr(dim, nlevels, seed, **kw)
    f = f or (lambda x, y=0.0, z=0.0: 2.0 + 0.5 * x + 0 * np.asarray(y))
    if len(names) == 1:
        data = {float(t): fill(hier, f, names=names, time=t) for t in times}
    else:
        data = {float(t): fill(hier, {n: (lambda *c, k=k: f(*c) + k) for k, n in
                                      enumerate(names)}, time=t) for t in times}
    return FakeBackend(data, quantity="fake"), hier


def counting(backend):
    """Wrap ``backend.read_into`` (instance attribute) to record (time, level, file_ids)."""
    calls = []
    original = backend.read_into

    def read_into(time, level, file_ids, out, dst_start):
        calls.append((float(time), level, np.asarray(file_ids).copy()))
        return original(time, level, file_ids, out, dst_start)

    backend.read_into = read_into
    return calls


def expected_value(ld, name, t=0.0):
    return vec_funcs(t)[name](*coords_of(ld.geom))


# --------------------------------------------------------------------------- construction


class TestLazyConstruction(unittest.TestCase):
    def test_nothing_read_until_data_is_used(self):
        backend, _ = fake_vector(times=(0.0, 0.5))
        B = AMRField.from_backend(backend)
        self.assertIsInstance(B, VectorField)
        text = repr(B)
        bx = B.x
        _ = B["By"], repr(bx), B.at(0.5), B.times(), B.hier(0.0), B.levels_data(0.5)
        _ = B.nbytes(), B.loaded_fraction, B.extent, B.dl, B.domain, B.ndim
        self.assertEqual(backend.calls, 0)
        self.assertIn("level 0:", text)
        self.assertGreater(B.nbytes(), 0)
        self.assertEqual(B.loaded_nbytes, 0)
        self.assertEqual(B.loaded_fraction, 0.0)
        # first use reads, then nothing more
        ld = bx.levels_data(0.0)[1]
        np.testing.assert_array_equal(ld.buffer("value"), expected_value(ld, "x"))
        self.assertEqual(backend.calls, 1)
        ld.buffer("value")
        self.assertEqual(backend.names_read, [("x",)])  # only the component used
        B.levels_data(0.0)[1].buffer("x")  # shares the read of the parent
        self.assertEqual(backend.calls, 1)
        B.levels_data(0.0)[1].buffer("y")  # y was not read yet
        self.assertEqual(backend.names_read, [("x",), ("y",)])

    def test_components_read_on_demand(self):
        """B.z reads z only, B.x then x only, a vector operation then only y."""
        backend, _ = fake_vector(nlevels=2, seed=6)
        B = AMRField.from_backend(backend)
        box = dict(x=(1.0, 4.0), y=(0.5, 3.0))
        gz = B.z.finest(box=box)
        self.assertEqual(set(backend.names_read), {("z",)})
        n = len(backend.names_read)
        gx = B.x.finest()  # whole domain: every patch of x
        self.assertEqual(set(backend.names_read[n:]), {("x",)})
        n = len(backend.names_read)
        B.z.load()  # the z patches outside the box only
        self.assertEqual(set(backend.names_read[n:]), {("z",)})
        n = len(backend.names_read)
        C = B * 2.0
        self.assertEqual(set(backend.names_read[n:]), {("y",)})
        ref = AMRField.from_backend(fake_vector(nlevels=2, seed=6)[0])
        np.testing.assert_array_equal(gz.data, ref.z.finest(box=box).data)
        np.testing.assert_array_equal(gx.data, ref.x.finest().data)
        for name in "xyz":
            for ld, lr in zip(C.levels_data(), ref.levels_data()):
                np.testing.assert_array_equal(ld.buffer(name), 2.0 * lr.buffer(name))
        self.assertEqual(B.loaded_fraction, 1.0)

    def test_component_load_estimate(self):
        backend, _ = fake_vector(nlevels=2, seed=6)
        B = AMRField.from_backend(backend)
        bz = B.z
        self.assertEqual(sum(_missing_bytes(ld, None) for ld in bz.levels_data()),
                         bz.nbytes())
        self.assertEqual(3 * bz.nbytes(), B.nbytes())
        bz.load()
        self.assertEqual(set(backend.names_read), {("z",)})
        self.assertEqual(bz.loaded_fraction, 1.0)
        self.assertAlmostEqual(B.loaded_fraction, 1 / 3)

    def test_class_from_components(self):
        backend, _ = fake_scalar()
        self.assertIsInstance(AMRField.from_backend(backend), ScalarField)
        backend, _ = fake_scalar(names=("xx", "xy", "xz", "yy", "yz", "zz"))
        P = AMRField.from_backend(backend)
        self.assertIsInstance(P, TensorField)
        with self.assertRaises(TypeError):
            VectorField.from_backend(backend)
        bv, _ = fake_vector()
        snaps = AMRField.from_backend(bv).snapshots()
        with self.assertRaisesRegex(TypeError, "make_field"):
            ScalarField(snaps)
        with self.assertRaises(TypeError):
            TensorField(snaps)
        scal = AMRField.from_backend(fake_scalar()[0]).snapshots()
        with self.assertRaises(TypeError):
            VectorField(scal)
        two, _ = fake_scalar(names=("a", "b"))
        f = AMRField.from_backend(two)
        self.assertIs(type(f), AMRField)
        self.assertIs(type(make_field(snaps)), VectorField)

    def test_geometry_shared_through_cache(self):
        backend, _ = fake_vector(times=(0.0, 1.0))
        cache = GeometryCache()
        B = AMRField.from_backend(backend, cache=cache)
        for g0, g1 in zip(B.hier(0.0).levels, B.hier(1.0).levels):
            self.assertIs(g0, g1)
        self.assertEqual(len(cache), len(B.hier(0.0).levels))


class TestExtent(unittest.TestCase):
    def check_extent(self, dim, seed, box):
        backend, hier = fake_vector(dim=dim, nlevels=3, seed=seed)
        pbox = PhysBox.parse(box, dim)
        B = AMRField.from_backend(backend, extent=box)
        self.assertEqual(B.extent, pbox)
        self.assertEqual(backend.calls, 0)
        # expected: per level, file-order patches selected by the PhysBox rule
        expected = []
        for g in hier.levels:
            lo, hi = g.file_order_boxes()
            sel = np.nonzero(pbox.selects(lo, hi, g.dl, g.domain_cells))[0]
            if len(sel) == 0:
                break
            expected.append((lo[sel], hi[sel], sel))
        sub = B.hier()
        self.assertEqual(len(sub.levels), len(expected))
        for g, (lo, hi, _sel) in zip(sub.levels, expected):
            glo, ghi = g.file_order_boxes()
            np.testing.assert_array_equal(glo, lo)
            np.testing.assert_array_equal(ghi, hi)
        B.load()
        requests = sorted(backend.requests, key=lambda r: r[1])
        for ilvl, ((_, level, ids), (_, _, sel)) in enumerate(zip(requests, expected)):
            self.assertEqual(level, ilvl)
            np.testing.assert_array_equal(ids, sel)
        self.assertEqual(len(backend.requests), len(expected))
        for ld in B.levels_data():
            for c in "xyz":
                np.testing.assert_allclose(ld.buffer(c), expected_value(ld, c), atol=1e-12)
        return B

    def test_extent_1d_2d_3d(self):
        self.check_extent(1, 3, dict(x=(3.1, 6.0)))
        self.check_extent(2, 5, dict(x=(2.0, 5.3), y=(1.0, 3.0)))
        self.check_extent(3, 6, dict(x=(0.5, 1.5), z=(0.4, 1.2)))

    def test_hierarchy_stops_at_first_empty_level(self):
        _, hier = fake_vector(dim=2, nlevels=3, seed=5)
        g0, g1 = hier.levels[0], hier.levels[1]
        covered = np.zeros(tuple(g1.domain_cells), bool)
        for lo, hi in zip(g1.lower, g1.upper):
            covered[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1] = True
        free = np.argwhere(~covered)
        # a box strictly inside a level-1 cell that no level-1 patch holds
        i, j = free[len(free) // 2]
        d = g1.dl
        box = dict(x=((i + 0.25) * d[0], (i + 0.75) * d[0]),
                   y=((j + 0.25) * d[1], (j + 0.75) * d[1]))
        B = self.check_extent(2, 5, box)
        self.assertEqual(B.hier().level_numbers(), (0,))
        self.assertGreater(g0.npatch, B.hier().levels[0].npatch)


class TestLoad(unittest.TestCase):

    def test_unload_reads_again(self):
        backend, _ = fake_vector(dim=1, nlevels=2)
        B = AMRField.from_backend(backend).load()
        n = backend.calls
        B.unload()
        self.assertEqual(B.loaded_nbytes, 0)
        B.levels_data()[0].buffer("x")
        self.assertEqual(backend.calls, n + 1)

    def test_file_backed_buffers_read_only_copy_writeable(self):
        backend, _ = fake_vector(dim=2, nlevels=2)
        B = AMRField.from_backend(backend)
        buf = B.levels_data()[0].buffer("x")
        with self.assertRaises(ValueError):
            buf[0] = 1.0
        with self.assertRaises(ValueError):
            B.x.levels_data()[0].buffer("value")[0] = 1.0
        C = B.copy()
        self.assertIsInstance(C, VectorField)
        cbuf = C.levels_data()[0].buffer("x")
        before = float(buf[0])
        cbuf[0] = before + 42.0
        self.assertEqual(float(B.levels_data()[0].buffer("x")[0]), before)
        self.assertFalse(np.shares_memory(cbuf, buf))
        D = copy.deepcopy(B)
        self.assertFalse(np.shares_memory(D.levels_data()[0].buffer("x"), buf))
        E = copy.copy(B)
        self.assertIsNot(E, B)
        self.assertIs(E.levels_data()[1], B.levels_data()[1])


class TestTimes(unittest.TestCase):
    def setUp(self):
        self.backend, _ = fake_vector(dim=1, nlevels=2, times=(0.0, 0.25, 0.5))

    def test_caller_order_and_at(self):
        B = AMRField.from_backend(self.backend, times=[0.5, 0.0])
        np.testing.assert_array_equal(B.times(), [0.5, 0.0])
        self.assertEqual([s.time for s in B.snapshots()], [0.5, 0.0])
        b0 = B.at(0.0 + 1e-13)  # within the tolerance
        self.assertIsInstance(b0, VectorField)
        np.testing.assert_array_equal(b0.times(), [0.0])
        self.assertIs(b0.snapshot(), B.snapshots()[1])
        ld = B.levels_data(0.5)[0]
        np.testing.assert_allclose(ld.buffer("y"), expected_value(ld, "y", 0.5))

    def test_time_not_found(self):
        B = AMRField.from_backend(self.backend)
        with self.assertRaises(TimeNotFoundError) as cm:
            B.at(0.3)
        self.assertIn("0.25", str(cm.exception))
        self.assertIsInstance(cm.exception, KeyError)
        self.assertIsInstance(cm.exception, ValueError)
        with self.assertRaises(TimeNotFoundError):
            AMRField.from_backend(self.backend, times=0.7)
        with self.assertRaisesRegex(ValueError, "twice"):
            AMRField.from_backend(self.backend, times=[0.0, 0.0])

    def test_single_time_rule(self):
        B = AMRField.from_backend(self.backend)
        for call in (B.snapshot, B.hier, B.levels_data, B.finest, B.mean, B.min, B.max,
                     B.integral, B.check_shared_nodes, lambda: B(x=1.0)):
            with self.assertRaises(ValueError) as cm:
                call()
            self.assertIn("0.25", str(cm.exception))
        self.assertEqual(self.backend.calls, 0)


class TestComponents(unittest.TestCase):
    def test_vector_components_share_buffers(self):
        backend, _ = fake_vector(dim=2, nlevels=2)
        B = AMRField.from_backend(backend)
        for key in ("x", "Bx", "EM_Bx"):
            c = B[key]
            self.assertIsInstance(c, ScalarField)
            self.assertEqual(c.quantity, "EM_B_x")
            self.assertEqual(c.names, ("value",))
        self.assertIsInstance(B.y, ScalarField)
        self.assertEqual(backend.calls, 0)
        bz = B.z
        a = bz.levels_data()[1].buffer("value")
        b = B.levels_data()[1].buffer("z")
        self.assertTrue(np.shares_memory(a, b))
        self.assertEqual(backend.calls, 1)
        self.assertGreater(B.loaded_nbytes, 0)
        self.assertEqual(bz.nbytes() * 3, B.nbytes())
        self.assertIs(bz["value"], bz)

    def test_in_memory_components_share_and_write(self):
        hier = make_amr(2, 2, 3)
        B = make_field([Snapshot(hier, tuple(fill(hier, vec_funcs())))], quantity="B")
        bx = B.x
        bx.levels_data()[0].buffer("value")[0] = -99.0
        self.assertEqual(B.levels_data()[0].buffer("x")[0], -99.0)

    def test_unknown_components(self):
        backend, _ = fake_vector(dim=1, nlevels=1)
        B = AMRField.from_backend(backend)
        for bad in ("w", "1x", "Bw", "", "_x", 3):
            with self.assertRaises(KeyError) as cm:
                B[bad]
            self.assertIn("('x', 'y', 'z')", str(cm.exception))
        with self.assertRaisesRegex(AttributeError, "components"):
            _ = B.w
        with self.assertRaises(KeyError):
            B.x["x"]

    def test_tensor_components(self):
        names = ("xx", "xy", "xz", "yy", "yz", "zz")
        backend, _ = fake_scalar(names=names)
        P = AMRField.from_backend(backend, quantity="protons_P")
        for key in ("xx", "Pxx", "protons_Pxx"):
            self.assertEqual(P[key].quantity, "protons_P_xx")
        self.assertEqual(P.yz.quantity, "protons_P_yz")
        with self.assertRaises(KeyError):
            P["Px"]
        ld = P["Pxy"].levels_data()[0]
        np.testing.assert_allclose(ld.buffer("value"),
                                   P.levels_data()[0].buffer("xy"))


class TestPickleRepr(unittest.TestCase):
    def test_fake_backed_pickle_without_data(self):
        backend, _ = fake_vector(dim=2, nlevels=3)
        B = AMRField.from_backend(backend, extent=dict(x=(1.0, 4.0)))
        blob = pickle.dumps(B.y)
        self.assertEqual(backend.calls, 0)
        C = pickle.loads(blob)
        self.assertIsInstance(C, ScalarField)
        self.assertEqual(C.extent, B.extent)
        self.assertEqual(C.loaded_nbytes, 0)
        for la, lb in zip(C.levels_data(), B.y.levels_data()):
            np.testing.assert_array_equal(la.buffer("value"), lb.buffer("value"))

    def test_in_memory_pickles_arrays(self):
        hier = make_amr(3, 2, 9, extent=dict(x=(0.2, 2.0)))
        B = make_field([Snapshot(hier, tuple(fill(hier, vec_funcs())))], quantity="B")
        C = pickle.loads(pickle.dumps(B))
        self.assertIsInstance(C, VectorField)
        self.assertEqual(C.extent, B.extent)
        for la, lb in zip(C.levels_data(), B.levels_data()):
            self.assertEqual(la.geom.key(), lb.geom.key())
            np.testing.assert_array_equal(la.buffer("z"), lb.buffer("z"))

    def test_repr(self):
        backend, _ = fake_vector(dim=2, nlevels=2, times=(0.0, 0.5))
        B = AMRField.from_backend(backend, stencil="node")
        text = repr(B)
        for part in ("VectorField", "EM_B", "('x', 'y', 'z')", "float64", "2 times",
                     "level 0:", "level 1:", "patches", "nodes", "extent whole domain",
                     "loaded 0.0 %", "fake.vtkhdf", "[fake]", "writer fake", "stencil node",
                     "dl (0.2, 0.2)"):
            self.assertIn(part, text)
        g = B.hier(0.0).levels[1]
        self.assertIn(f"level 1: {g.npatch} patches", text)
        self.assertIn("time 0.5", repr(B.at(0.5)))
        self.assertEqual(backend.calls, 0)
        self.assertIn("in memory", repr(B.at(0.5) + 1.0))


class TestDelegations(unittest.TestCase):
    def setUp(self):
        hier = make_amr(2, 2, 11)
        self.hier = hier
        self.B = make_field([Snapshot(hier, tuple(fill(hier, vec_funcs())))], quantity="B")

    def test_finest_and_reductions(self):
        g = self.B.y.finest()
        X, Y = np.meshgrid(g.x, g.y, indexing="ij")
        np.testing.assert_allclose(g.data, vec_funcs()["y"](X, Y), atol=1e-12)
        try:
            from pyphare.pharesee.amr.ops import reduce  # noqa: F401
        except ImportError:
            self.skipTest("ops/reduce.py not written yet")
        one = self.B.x * 0.0 + 2.0
        area = float(np.prod(self.B.domain))
        self.assertAlmostEqual(one.integral(), 2.0 * area, places=9)
        self.assertAlmostEqual(one.mean(), 2.0, places=12)
        lo, hi = self.B.x.min(), self.B.x.max()
        allx = np.concatenate([ld.buffer("value") for ld in self.B.x.levels_data()])
        self.assertAlmostEqual(lo, allx.min())
        self.assertAlmostEqual(hi, allx.max())
        per = self.B.x.mean(per_level=True)
        self.assertEqual(sorted(per), [0, 1])

    def test_check_shared_nodes(self):
        for periodic in (True, False):
            hier = make_amr(1, 1, 12, periodic=periodic)
            # periodic images hold equal values only for periodic data (domain 64 * 0.2)
            lds = fill(hier, lambda x: np.sin(2 * np.pi * x / 12.8))
            f = make_field([Snapshot(hier, tuple(lds))])
            self.assertLess(f.check_shared_nodes()[0], 1e-14)
            g = hier.levels[0]
            first = int(g.store_index[np.argmin(g.file_order_boxes()[0][:, 0])])
            c = f.copy()
            c.levels_data()[0].buffer("value")[g.offset[first]] += 0.5  # node 0
            want = 0.5 if periodic else 0.0
            self.assertAlmostEqual(c.check_shared_nodes()[0], want, places=12)
            inner = int(np.argmax(g.lower[:, 0]))  # its node `lower` is shared with a patch
            c2 = f.copy()
            c2.levels_data()[0].buffer("value")[g.offset[inner]] -= 0.25
            self.assertAlmostEqual(c2.check_shared_nodes()[0], 0.25, places=12)

    def test_check_shared_nodes_2d_corners(self):
        hier = make_amr(2, 2, 13, periodic=True)
        Lx, Ly = 48 * 0.2, 40 * 0.2

        def per(x, y):
            return np.cos(2 * np.pi * x / Lx) * np.sin(2 * np.pi * y / Ly) + x * 0

        f = make_field([Snapshot(hier, tuple(fill(hier, per)))])
        spread = f.check_shared_nodes()
        self.assertEqual(sorted(spread), [0, 1])
        self.assertLess(max(spread.values()), 1e-14)
        c = f.copy()
        g = hier.levels[1]
        c.levels_data()[1].buffer("value")[g.offset[1] - 1] += 0.75  # a corner node
        self.assertAlmostEqual(c.check_shared_nodes()[1], 0.75)


# --------------------------------------------------------------------------- files


@unittest.skipUnless(RUNS, "build/sample_data/*/{vtk,h5} not found")
class TestRealSamples(unittest.TestCase):
    def test_pickle_lazy_real_field(self):
        run = "2d_tagged" if "2d_tagged" in RUNS else RUNS[0]
        for fmt, ext in (("vtk", "vtkhdf"), ("h5", "h5")):
            with self.subTest(fmt=fmt):
                backend = open_backend(str(SAMPLE / run / fmt / f"EM_B.{ext}"))
                t = float(backend.times()[-1])
                B = AMRField.from_backend(backend, times=t)
                blob = pickle.dumps(B)
                self.assertLess(len(blob), 50_000)
                C = pickle.loads(blob)
                self.assertEqual(C.loaded_nbytes, 0)
                self.assertEqual(len(C.hier().levels), len(B.hier().levels))
                for la, lb in zip(C.levels_data(), B.levels_data()):
                    np.testing.assert_array_equal(la.buffer("y"), lb.buffer("y"))
                bz = pickle.loads(pickle.dumps(B.z))
                np.testing.assert_array_equal(bz.levels_data()[-1].buffer("value"),
                                              B.levels_data()[-1].buffer("z"))

    def test_real_extent_values(self):
        run = "2d_tagged" if "2d_tagged" in RUNS else RUNS[0]
        path = str(SAMPLE / run / "vtk" / "EM_B.vtkhdf")
        backend = open_backend(path)
        t = float(backend.times()[-1])
        calls = counting(backend)
        dom = AMRField.from_backend(open_backend(path), times=t)
        box = {a: (0.3 * d, 0.45 * d) for a, d in zip("xyz", dom.domain)}
        B = AMRField.from_backend(backend, times=t, extent=box)
        B.load()
        pbox = PhysBox.parse(box, B.ndim)
        for g_sub, ld_full in zip(B.hier().levels, dom.levels_data()):
            gf = ld_full.geom
            lo, hi = gf.file_order_boxes()
            sel = np.nonzero(pbox.selects(lo, hi, gf.dl, gf.domain_cells))[0]
            got = np.concatenate([ids for _, lvl, ids in calls if lvl == g_sub.level])
            np.testing.assert_array_equal(got, sel)
            ld_sub = B.levels_data()[g_sub.level]
            for slot in range(g_sub.npatch):
                fid = sel[g_sub.file_index[slot]]
                np.testing.assert_array_equal(ld_sub.view(slot, "x"),
                                              ld_full.view(gf.store_index[fid], "x"))


if __name__ == "__main__":
    unittest.main()
