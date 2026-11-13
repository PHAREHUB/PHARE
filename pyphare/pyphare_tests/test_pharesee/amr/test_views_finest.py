#!/usr/bin/env python3
"""Tests of pharesee.amr.views.finest: regular-grid view of an AMR snapshot."""

import gc
import itertools
import os
import tracemalloc
import unittest
import warnings
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.boxes import PhysBox
from pyphare.pharesee.amr.errors import PhareseeMemoryError
from pyphare.pharesee.amr.store import LevelData, Snapshot
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr
from pyphare.pharesee.amr.views.finest import NO_DATA, FinestGrid, finest


def linear(*c):
    return 1.0 + sum((a + 1.5) * x for a, x in enumerate(c))


def rough(*c):
    """Deterministic 'noisy' field with a jump: any interpolation overshoot shows."""
    v = np.sin(37.0 * c[0] + 1.0) * 3.0 + (c[0] > 3.1) * 5.0
    for a, x in enumerate(c[1:]):
        v = v + np.cos((53.0 + 7 * a) * x * (1 + 0.3 * c[0]))
    return v


def snapshot(hier, f=linear, dtype=np.float64, lazy=False, names=None):
    """Snapshot of ``f``; lazy=True -> file-backed through a FakeBackend (returned too)."""
    src = fill(hier, f, names=names, dtype=dtype)
    if not lazy:
        return Snapshot(hier, tuple(src)), None
    fb = FakeBackend(src)
    lds = tuple(LevelData.from_backend(fb, g, 0.0) for g in hier.levels)
    return Snapshot(hier, lds), fb


def node_range(box, hier, L):
    """Expected output node range at level L (the documented rule, recomputed here)."""
    dl = hier.levels[0].dl / 2**L
    N = hier.levels[0].domain_cells * 2**L
    lo, hi = PhysBox.parse(box, hier.dim).clipped(N * dl) if box is not None else (
        np.zeros(hier.dim), N * dl)
    q0, q1 = lo / dl, hi / dl
    q0 = np.where(np.abs(q0 - np.round(q0)) < 1e-9, np.round(q0), q0)
    q1 = np.where(np.abs(q1 - np.round(q1)) < 1e-9, np.round(q1), q1)
    return np.floor(q0).astype(np.int64), np.ceil(q1).astype(np.int64)


def reference(snap, nlo, nhi, L, name=None):
    """Brute force, patch by patch: value and level of every output node that is a node of
    a level <= L holding it in a patch (finest such level wins); NaN / -1 elsewhere."""
    shape = tuple(nhi - nlo + 1)
    ref = np.full(shape, np.nan)
    lev = np.full(shape, -1, np.int64)
    for geom, ld in zip(snap.hier.levels, snap.levels):
        l = geom.level
        if l > L:
            continue
        r = 2 ** (L - l)
        name_ = name or ld.names[0]
        for slot in range(geom.npatch):  # reference loop over patches (test only)
            idx = []
            for a in range(geom.dim):
                k0 = max(geom.lower[slot, a] * r, -(-nlo[a] // r) * r)
                k1 = min((geom.upper[slot, a] + 1) * r, nhi[a])
                idx.append(np.arange(k0, k1 + 1, r))
            if any(len(i) == 0 for i in idx):
                continue
            view = ld.view(slot, name_)
            local = np.ix_(*[i // r - geom.lower[slot, a] for a, i in enumerate(idx)])
            out = np.ix_(*[i - nlo[a] for a, i in enumerate(idx)])
            ref[out] = view[local]
            lev[out] = l
    return ref, lev


CASES = [  # (dim, nlevels, seed, periodic)
    (1, 1, 0, False), (1, 2, 1, True), (1, 3, 2, False),
    (2, 1, 3, True), (2, 2, 4, False), (2, 3, 5, True), (2, 3, 6, False),
    (3, 2, 7, False), (3, 3, 8, True),
]


def boxes_for(hier, rng):
    """Whole domain, a random box, a box narrower than a cell, a box on the domain edge."""
    dom = hier.domain
    dim = hier.dim
    dl = hier.levels[-1].dl
    out = [None]
    a = rng.uniform(0.1, 0.5, dim) * dom
    b = a + rng.uniform(0.2, 0.45, dim) * dom
    out.append({k: (a[i], b[i]) for i, k in enumerate("xyz"[:dim])})
    c = rng.uniform(0.2, 0.8, dim) * dom
    out.append({k: (c[i], c[i] + 0.3 * dl[i]) for i, k in enumerate("xyz"[:dim])})
    out.append({k: (0.0, 0.3 * dom[i]) if i == 0 else (0.6 * dom[i], dom[i])
                for i, k in enumerate("xyz"[:dim])})
    return out


class TestFinestValues(unittest.TestCase):
    def test_linear_field_exact_everywhere_cubic(self):
        rng = np.random.default_rng(0)
        for dim, nl, seed, per in CASES:
            hier = make_amr(dim, nl, seed=seed, periodic=per, two_node=True)
            snap, _ = snapshot(hier)
            for box in boxes_for(hier, rng):
                with self.subTest(dim=dim, nl=nl, box=box):
                    g = finest(snap, box=box)
                    X = np.meshgrid(*g.coords, indexing="ij")
                    np.testing.assert_allclose(g.data, linear(*X), rtol=0, atol=1e-12)

    def test_fine_values_bit_identical_and_level_map(self):
        rng = np.random.default_rng(1)
        for (dim, nl, seed, per), fill_ in itertools.product(CASES, ("cubic", "nearest")):
            hier = make_amr(dim, nl, seed=seed, periodic=per, two_node=True)
            for dtype in (np.float64, np.float32):
                snap, _ = snapshot(hier, rough, dtype=dtype)
                for box in boxes_for(hier, rng):
                    for max_level in (None, 0, 1):
                        L = hier.finest_level if max_level is None else min(max_level, nl - 1)
                        with self.subTest(dim=dim, nl=nl, box=box, ml=max_level, fill=fill_):
                            g = finest(snap, box=box, max_level=max_level, fill=fill_)
                            nlo, nhi = node_range(box, hier, L)
                            ref, lev = reference(snap, nlo, nhi, L)
                            self.assertEqual(g.level, L)
                            self.assertEqual(g.data.dtype, dtype)
                            self.assertEqual(g.data.shape, ref.shape)
                            self.assertEqual(g.level_map.dtype, np.uint8)
                            held = lev >= 0
                            # bit-identical: compare the raw bytes
                            np.testing.assert_array_equal(
                                g.data[held].view(np.uint8), ref[held].astype(dtype).view(np.uint8))
                            np.testing.assert_array_equal(g.level_map[held], lev[held])
                            gap = g.level_map[~held]
                            self.assertTrue((gap < L).all() if L > 0 else gap.size == 0)
                            self.assertFalse(np.isnan(g.data).any())
                            self.assertFalse((g.level_map == NO_DATA).any())

    def coarse_nodes(self, hier, snap):
        """Level-0 values on the whole level-0 node grid (level 0 tiles the domain)."""
        N0 = hier.levels[0].domain_cells
        ref, lev = reference(snap, np.zeros(hier.dim, np.int64), N0, 0)
        self.assertTrue((lev == 0).all())
        return ref

    def test_no_new_extrema_cubic(self):
        """Two levels: every gap value lies within [min, max] of the level-0 nodes of its
        enclosing level-0 cell (1D/2D/3D, a field with a jump)."""
        for dim, seed in ((1, 1), (2, 2), (2, 9), (3, 3)):
            hier = make_amr(dim, 2, seed=seed)
            snap, _ = snapshot(hier, rough)
            c0 = self.coarse_nodes(hier, snap)
            g = finest(snap, fill="cubic")
            nlo, nhi = node_range(None, hier, 1)
            _, lev = reference(snap, nlo, nhi, 1)
            k = np.nonzero(lev < 0)  # gap nodes (fine indices, nlo = 0 here)
            corners = []
            for bits in itertools.product((0, 1), repeat=dim):
                idx = tuple(-(-ki // 2) if bit else ki // 2 for ki, bit in zip(k, bits))
                corners.append(c0[idx])
            lo, hi = np.min(corners, axis=0), np.max(corners, axis=0)
            v = g.data[k]
            self.assertGreater(len(v), 0)
            self.assertTrue(((v >= lo) & (v <= hi)).all(), f"{dim}D: overshoot")
            # and the interpolation is not just a copy: some gap values are strictly inside
            self.assertTrue(((v > lo) & (v < hi)).any())

    def test_no_new_extrema_three_levels(self):
        hier = make_amr(2, 3, seed=5)
        snap, _ = snapshot(hier, rough)
        allv = np.concatenate([ld.buffer("value") for ld in snap.levels])
        g = finest(snap)
        self.assertGreaterEqual(g.data.min(), allv.min())
        self.assertLessEqual(g.data.max(), allv.max())

    def test_nearest_takes_lower_coarse_node(self):
        for dim, seed in ((1, 1), (2, 2), (3, 3)):
            hier = make_amr(dim, 2, seed=seed)
            snap, _ = snapshot(hier, rough)
            c0 = self.coarse_nodes(hier, snap)
            g = finest(snap, fill="nearest")
            nlo, nhi = node_range(None, hier, 1)
            _, lev = reference(snap, nlo, nhi, 1)
            k = np.nonzero(lev < 0)
            np.testing.assert_array_equal(g.data[k], c0[tuple(ki // 2 for ki in k)])
            # 3 levels: nearest never creates a value that is not in the data
            hier3 = make_amr(dim, 3, seed=seed)
            snap3, _ = snapshot(hier3, rough)
            allv = np.concatenate([ld.buffer("value") for ld in snap3.levels])
            self.assertTrue(np.isin(finest(snap3, fill="nearest").data, allv).all())

    def test_box_geometry(self):
        hier = make_amr(2, 2, seed=4, dl0=0.2)  # level 1: dl = 0.1, 96 x 80 cells
        snap, _ = snapshot(hier)
        g = finest(snap, box=dict(x=(1.23, 1.25), y=(2.0, 2.0)))
        np.testing.assert_allclose(g.x, [1.2, 1.3])
        np.testing.assert_allclose(g.y, [2.0])
        self.assertEqual(g.data.shape, (2, 1))
        np.testing.assert_allclose(g.dl, [0.1, 0.1])
        self.assertEqual(g.box, PhysBox((g.x[0], g.y[0]), (g.x[-1], g.y[-1])))
        # extent: pixel edges, matplotlib order
        np.testing.assert_allclose(g.extent, (1.15, 1.35, 1.95, 2.05))
        # whole domain, x free in the dict
        g = finest(snap, box=dict(y=(0.0, 1.0)))
        np.testing.assert_allclose([g.x[0], g.x[-1]], [0.0, hier.domain[0]])
        self.assertEqual(len(g.x), 97)
        # coordinates are (first node + i) * dl exactly
        np.testing.assert_array_equal(g.y, np.arange(11) * 0.1)
        with self.assertRaises(ValueError):
            finest(snap, box=dict(x=(100.0, 200.0)))
        with self.assertRaises(TypeError):
            finest(snap, box=((0, 1), (0, 1)))

    def test_extent_is_default_box_and_clips(self):
        hier = make_amr(2, 2, seed=4, extent=dict(x=(1.0, 3.0), y=(2.0, 4.0)))
        snap, _ = snapshot(hier)
        g = finest(snap)
        np.testing.assert_allclose([g.x[0], g.x[-1], g.y[0], g.y[-1]], [1.0, 3.0, 2.0, 4.0])
        g = finest(snap, box=dict(x=(0.0, 2.0)))
        np.testing.assert_allclose([g.x[0], g.x[-1], g.y[0], g.y[-1]], [1.0, 2.0, 2.0, 4.0])
        with self.assertRaises(ValueError):
            finest(snap, box=dict(x=(5.0, 6.0)))


class TestFinestLaziness(unittest.TestCase):
    def expected_reads(self, hier, nlo, nhi, L):
        """file-order ids per level whose closed node box meets that level's window."""
        wins = {L: (nlo, nhi)}
        for l in range(L - 1, -1, -1):
            lo, hi = wins[l + 1]
            N = hier.levels[0].domain_cells * 2**l
            wins[l] = (np.maximum(lo // 2 - 1, 0), np.minimum(-(-hi // 2) + 1, N))
        out = {}
        for geom in hier.levels:
            if geom.level > L:
                continue
            lo, hi = wins[geom.level]
            lower, upper = geom.file_order_boxes()
            touch = ((lower <= hi) & (upper + 1 >= lo)).all(axis=1)
            out[geom.level] = np.nonzero(touch)[0]
        return out

    def test_reads_only_patches_touching_windows(self):
        rng = np.random.default_rng(3)
        for dim, nl, seed, per in CASES:
            hier = make_amr(dim, nl, seed=seed, periodic=per)
            for box in boxes_for(hier, rng):
                for max_level in (None, 0, 1):
                    snap, fb = snapshot(hier, lazy=True)
                    L = hier.finest_level if max_level is None else min(max_level, nl - 1)
                    with self.subTest(dim=dim, nl=nl, box=box, ml=max_level):
                        g = finest(snap, box=box, max_level=max_level)
                        nlo, nhi = node_range(box, hier, L)
                        exp = self.expected_reads(hier, nlo, nhi, L)
                        for geom in hier.levels:
                            got = fb.reads.get((0.0, geom.level))
                            got = (np.nonzero(got)[0] if got is not None
                                   else np.zeros(0, np.int64))
                            if geom.level > L:
                                self.assertEqual(len(got), 0, "level above max_level read")
                            else:
                                np.testing.assert_array_equal(got, exp[geom.level])
                        # one read call per level, each patch read once
                        self.assertEqual(fb.calls, len([v for v in exp.values() if len(v)]))
                        for counts in fb.reads.values():
                            self.assertLessEqual(counts.max(), 1)
                        # same values as from in-memory data
                        ref, _ = snapshot(hier)
                        np.testing.assert_array_equal(
                            g.data, finest(ref, box=box, max_level=max_level).data)

    def test_small_box_reads_few_patches(self):
        hier = make_amr(2, 3, seed=5, domain_cells=(96, 80), patch_cells=8)
        snap, fb = snapshot(hier, lazy=True)
        finest(snap, box=dict(x=(9.0, 9.4), y=(7.0, 7.3)))
        total = sum(int(c.sum()) for c in fb.reads.values())
        self.assertLessEqual(total, 20)
        self.assertLess(total, sum(g.npatch for g in hier.levels) / 20)


class TestFinestApi(unittest.TestCase):
    def setUp(self):
        self.hier = make_amr(2, 2, seed=4)
        self.vec, _ = snapshot(self.hier, {"x": linear, "y": rough, "z": lambda *c: c[1]})
        self.sca, _ = snapshot(self.hier, rough)

    def test_one_component(self):
        g = finest(self.sca)
        self.assertIsInstance(g, FinestGrid)
        self.assertIsInstance(g.data, np.ndarray)
        self.assertIs(g["value"], g.data)
        data, (x, y) = g
        self.assertIs(data, g.data)
        self.assertEqual((len(x), len(y)), data.shape)
        np.testing.assert_array_equal(np.asarray(g), g.data)
        self.assertEqual(np.asarray(g, dtype=np.float32).dtype, np.float32)
        with self.assertRaises(AttributeError):
            _ = g.z
        self.assertTrue(g.data.flags.f_contiguous)
        self.assertEqual(g.time, 0.0)
        self.assertIn("FinestGrid", repr(g))

    def test_components(self):
        g = finest(self.vec)
        self.assertEqual(set(g.data), {"x", "y", "z"})
        ref = finest(self.sca)
        np.testing.assert_array_equal(g["y"], ref.data)
        np.testing.assert_array_equal(g.level_map, ref.level_map)
        with self.assertRaises(TypeError):
            np.asarray(g)
        g2 = finest(self.vec, names=["z"])
        self.assertEqual(set(g2.data), {"z"})
        np.testing.assert_array_equal(g2["z"], g["z"])
        with self.assertRaises(KeyError):
            finest(self.vec, names="w")

    def test_dtype_and_arguments(self):
        g = finest(self.sca, dtype=np.float32)
        self.assertEqual(g.data.dtype, np.float32)
        np.testing.assert_allclose(g.data, finest(self.sca).data, rtol=1e-6, atol=1e-6)
        with self.assertRaises(TypeError):
            finest(self.sca, dtype=np.int32)
        with self.assertRaises(ValueError):
            finest(self.sca, fill="linear")
        with self.assertRaises(ValueError):
            finest(self.sca, max_level=-1)
        self.assertEqual(finest(self.sca, max_level=7).level, 1)

    def test_one_and_three_dims(self):
        g1 = finest(snapshot(make_amr(1, 2, seed=1))[0])
        self.assertEqual(g1.data.ndim, 1)
        with self.assertRaises(AttributeError):
            _ = g1.y
        self.assertEqual(len(g1.extent), 2)
        g3 = finest(snapshot(make_amr(3, 2, seed=1))[0])
        self.assertEqual(len(g3.z), g3.data.shape[2])
        with self.assertRaises(ValueError):
            g3.imshow()

    def test_imshow(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        g = finest(self.sca, box=dict(x=(1.0, 3.0), y=(2.0, 5.0)))
        im = g.imshow(ax)
        np.testing.assert_array_equal(im.get_array(), g.data.T)
        self.assertEqual(tuple(im.get_extent()), g.extent)
        self.assertEqual(im.origin, "lower")
        self.assertEqual(im.get_interpolation(), "nearest")
        gv = finest(self.vec)
        with self.assertRaises(ValueError):
            gv.imshow(ax)
        np.testing.assert_array_equal(gv.imshow(ax, name="y").get_array(), gv["y"].T)
        plt.close(fig)

    def test_missing_level0_data_is_nan(self):
        # drop one level-0 patch: its interior shows NaN and NO_DATA
        from pyphare.pharesee.amr.geometry import HierGeometry, LevelGeometry

        hier = make_amr(2, 1, seed=2)
        g0 = hier.levels[0]
        lower, upper = g0.file_order_boxes()
        big = int(np.argmax(np.prod(upper - lower + 1, axis=1)))
        keep = np.arange(len(lower)) != big
        geom = LevelGeometry.from_file_order(0, lower[keep], upper[keep], g0.dl,
                                             g0.domain_cells, g0.periodic)
        h = HierGeometry([geom])
        snap, _ = snapshot(h)
        g = finest(snap)
        inner = tuple(slice(lower[big, a] + 1, upper[big, a] + 1) for a in range(2))
        self.assertTrue(np.isnan(g.data[inner]).all())
        self.assertTrue((g.level_map[inner] == NO_DATA).all())
        self.assertFalse(np.isnan(g.data[g.level_map != NO_DATA]).any())


class TestFinestMemory(unittest.TestCase):
    def test_peak_memory_large_output(self):
        """>= 50 M output nodes: traced numpy peak <= 1.6 x output (+ the 1-byte level map)."""
        hier = make_amr(2, 2, seed=11, domain_cells=(5000, 2500), patch_cells=250,
                        fraction=(0.2, 0.3))
        snap, _ = snapshot(hier, lambda x, y: x + y, dtype=np.float32)
        gc.collect()
        tracemalloc.start()
        try:
            g = finest(snap)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        n = g.data.size
        self.assertGreaterEqual(n, 50_000_000)
        out = g.data.nbytes
        self.assertLessEqual(peak, 1.6 * out + g.level_map.nbytes,
                             f"peak {peak / 1e6:.0f} MB for output {out / 1e6:.0f} MB")
        del g
        gc.collect()

    def test_budget(self):
        snap, fb = snapshot(make_amr(2, 2, seed=4), lazy=True)
        with mock.patch.dict(os.environ, {"PHARESEE_MEMORY_LIMIT": "1K"}):
            with self.assertRaises(PhareseeMemoryError) as cm:
                finest(snap)
            self.assertIn("force=True", str(cm.exception))
            self.assertEqual(fb.calls, 0, "nothing may be read before the memory check")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                g = finest(snap, force=True)
            self.assertTrue(any(issubclass(x.category, ResourceWarning) for x in w))
        ref, _ = snapshot(make_amr(2, 2, seed=4))
        np.testing.assert_array_equal(g.data, finest(ref).data)


if __name__ == "__main__":
    unittest.main()
