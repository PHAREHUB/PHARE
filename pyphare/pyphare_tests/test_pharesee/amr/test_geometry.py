#!/usr/bin/env python3
"""Tests of pharesee.amr.geometry: store order, views, runs, end layers, faces, coverage."""

import os
import pickle
import unittest
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.boxes import (
    BinIndex,
    box_flat_indices,
    fortran_strides,
    periodic_shifts,
)
from pyphare.pharesee.amr.geometry import (
    GeometryCache,
    HierGeometry,
    LevelGeometry,
    build_faces,
)
from pyphare.pharesee.amr.testing import make_amr, make_large_amr


def layouts():
    """(name, HierGeometry) covering 1D/2D/3D, 2-node patches, singleton shapes, periodicity."""
    out = []
    for dim in (1, 2, 3):
        out.append((f"{dim}d", make_amr(dim, 3, seed=dim)))
        out.append((f"{dim}d-2node-periodic", make_amr(dim, 2, seed=10 + dim, two_node=True,
                                                       periodic=True)))
        out.append((f"{dim}d-singletons", make_amr(dim, 2, seed=20 + dim, singletons=True,
                                                   periodic=[True] + [False] * (dim - 1))))
    return out


def node_index(geom, axis):
    """Global integer node index along axis of every flat node."""
    return np.rint(geom.node_coords(axis) / geom.dl[axis]).astype(np.int64)


class TestStoreOrder(unittest.TestCase):
    def test_round_trip_and_sort(self):
        for name, hier in layouts():
            for g in hier.levels:
                n = g.npatch
                np.testing.assert_array_equal(g.file_index[g.store_index], np.arange(n))
                np.testing.assert_array_equal(g.store_index[g.file_index], np.arange(n))
                flo, fhi = g.file_order_boxes()
                np.testing.assert_array_equal(flo[g.file_index], g.lower)
                np.testing.assert_array_equal(fhi[g.file_index], g.upper)
                # lexicographic by (nx, ny, nz), stable within equal shapes
                keys = [tuple(s) for s in g.shape]
                self.assertEqual(keys, sorted(keys), name)
                for row in g.groups:
                    fi = g.file_index[row["p0"] : row["p1"]]
                    self.assertTrue((np.diff(fi) > 0).all(), name)
                    self.assertTrue((g.shape[row["p0"] : row["p1"]] == row["shape"]).all())
                self.assertEqual(g.offset[-1], g.nnodes)
                self.assertEqual(g.groups["n1"][-1], g.nnodes)
                self.assertFalse(g.lower.flags.writeable)

    def test_from_file_order_rejects_outside(self):
        with self.assertRaises(ValueError):
            LevelGeometry.from_file_order(0, [[0]], [[10]], [0.1], [10], [False])
        with self.assertRaises(ValueError):
            LevelGeometry.from_file_order(0, [[3]], [[2]], [0.1], [10], [False])


class TestViews(unittest.TestCase):
    def test_view_and_stack_share_memory(self):
        for name, hier in layouts():
            g = hier.levels[0]
            flat = np.arange(g.nnodes, dtype=np.float64)
            for gi, row in enumerate(g.groups):
                st = g.stack(flat, gi)
                self.assertTrue(np.shares_memory(st, flat))
                self.assertEqual(st.shape, tuple(row["shape"]) + (row["p1"] - row["p0"],))
                for k, slot in enumerate(range(row["p0"], row["p1"])):
                    v = g.view(flat, slot)
                    self.assertTrue(np.shares_memory(v, flat))
                    np.testing.assert_array_equal(st[..., k], v, name)
            # writes through a view land in the flat buffer, x fastest
            slot = g.npatch - 1
            v = g.view(flat, slot)
            v[(1,) + (0,) * (g.dim - 1)] = -7.0
            self.assertEqual(flat[g.offset[slot] + 1], -7.0)

    def test_coords(self):
        g = make_amr(2, 1, seed=5).levels[0]
        x = g.node_coords(0)
        y = g.node_coords(1)
        for slot in (0, g.npatch // 2, g.npatch - 1):
            vx = g.view(x, slot)
            vy = g.view(y, slot)
            np.testing.assert_allclose(vx[:, 0], g.coords(slot, 0))
            np.testing.assert_allclose(vy[0, :], g.coords(slot, 1))
            self.assertAlmostEqual(g.coords(slot, 0)[0], g.lower[slot, 0] * g.dl[0])


class TestRunsAndEnds(unittest.TestCase):
    def test_runs_centred_difference(self):
        for name, hier in layouts():
            for g in hier.levels:
                c = [g.node_coords(a) for a in range(g.dim)]
                f = sum((a + 1) * c[a] ** 2 for a in range(g.dim)) + c[0] * c[-1]
                for axis in range(g.dim):
                    exact = 2 * (axis + 1) * c[axis] + (c[-1] if axis == 0 else 0)
                    if axis == g.dim - 1:
                        exact = exact + c[0]
                    if g.dim == 1:
                        exact = 2 * c[0] + 2 * c[0]
                    d = np.full(g.nnodes, np.nan)
                    h = g.dl[axis]
                    stride, n0, n1 = g.runs(axis)
                    for s, a, b in zip(stride, n0, n1):
                        d[a + s : b - s] = (f[a + 2 * s : b] - f[a : b - 2 * s]) / (2 * h)
                    lo, hi = g.end_nodes(axis)
                    interior = np.ones(g.nnodes, bool)
                    interior[lo] = interior[hi] = False
                    np.testing.assert_allclose(d[interior], exact[interior], rtol=1e-9,
                                               atol=1e-9, err_msg=f"{name} L{g.level} {axis}")
                    if axis > 0:
                        self.assertLess(len(stride), g.npatch + 1)

    def test_end_layers_brute_force(self):
        for name, hier in layouts():
            g = hier.levels[0]
            flat = np.arange(g.nnodes)
            for axis in range(g.dim):
                lo, hi, stride, count = g.end_layers(axis)
                ref_lo, ref_hi = [], []
                for slot in range(g.npatch):
                    v = g.view(flat, slot)
                    ref_lo.append(np.take(v, 0, axis=axis).ravel(order="F"))
                    ref_hi.append(np.take(v, -1, axis=axis).ravel(order="F"))
                np.testing.assert_array_equal(lo, np.concatenate(ref_lo), name)
                np.testing.assert_array_equal(hi, np.concatenate(ref_hi), name)
                np.testing.assert_array_equal(hi, lo + (count - 1) * stride)
                self.assertTrue((count >= 2).all())


def brute_faces(g, axis):
    """{face flat index: set of acceptable outer global node positions} by direct search."""
    shifts = periodic_shifts(g.domain_cells, g.periodic)
    node_lo = g.lower[:, None, :] + shifts[None, :, :]  # (npatch, nshift, dim)
    node_hi = g.upper[:, None, :] + 1 + shifts[None, :, :]
    e = np.eye(g.dim, dtype=np.int64)[axis]
    flat = np.arange(g.nnodes)
    pos = np.stack([node_index(g, a) for a in range(g.dim)], axis=1)
    expected = {}
    for slot in range(g.npatch):
        v = g.view(flat, slot)
        for side, sign in ((0, -1), (-1, 1)):
            for fidx in np.take(v, side, axis=axis).ravel():
                G = pos[fidx]
                O = G + sign * e
                inside = (node_lo <= G) & (node_hi >= G) & (node_lo <= O) & (node_hi >= O)
                if inside.all(axis=2).any():
                    expected[int(fidx)] = sign
    return expected, pos


def reference_faces(g, axis):
    """The first vectorised build_faces (sorted pairs, two box_flat_indices, np.unique),
    kept verbatim as the reference the optimised one must match bit for bit."""
    n, dim = g.npatch, g.dim
    e = np.zeros(dim, dtype=np.int64)
    e[axis] = 1
    lower, upper, shape, offset = g.lower, g.upper, g.shape, g.offset[:-1]
    strides = fortran_strides(shape)
    edge_lo, edge_hi = lower.copy(), upper + 1
    edge_hi[:, axis] = upper[:, axis]
    index = BinIndex(edge_lo, edge_hi)
    q_lo = np.concatenate([lower, lower])
    q_hi = np.concatenate([upper + 1, upper + 1])
    q_lo[:n, axis] = q_hi[:n, axis] = upper[:, axis] + 1
    q_lo[n:, axis] = q_hi[n:, axis] = lower[:, axis] - 1
    shifts = periodic_shifts(g.domain_cells, g.periodic)
    pairs = index.pairs(q_lo, q_hi, shifts)
    hi_side = pairs.a < n
    p = np.where(hi_side, pairs.a, pairs.a - n)
    q = pairs.b
    s = shifts[pairs.shift]
    size = pairs.hi - pairs.lo + 1
    outer_lo = pairs.lo + np.where(hi_side[:, None], e, 0) - lower[q]
    face_lo = pairs.lo - s + np.where(hi_side[:, None], 0, e) - lower[p]
    run = int(np.argmax(size.sum(axis=0)))
    outer = box_flat_indices(offset[q], shape[q], outer_lo, size, run_axis=run)
    face = box_flat_indices(offset[p], shape[p], face_lo, size, run_axis=run)
    counts = np.prod(size, axis=1)
    sign = np.repeat(np.where(hi_side, 1, -1).astype(np.int8), counts)
    patch = np.repeat(p, counts)
    inner = face - sign.astype(np.int64) * np.repeat(strides[p, axis], counts)
    face, first = np.unique(face, return_index=True)
    inner, outer, sign, patch = inner[first], outer[first], sign[first], patch[first]
    nodes_per_face = np.prod(np.delete(shape, axis, axis=1), axis=1) if dim > 1 else np.ones(n)
    listed_hi = np.bincount(patch[sign > 0], minlength=n)
    listed_lo = np.bincount(patch[sign < 0], minlength=n)
    return dict(face=face, inner=inner, outer=outer, sign=sign, patch=patch,
                one_sided_lo=listed_lo < nodes_per_face, one_sided_hi=listed_hi < nodes_per_face)


class TestFaces(unittest.TestCase):
    def check(self, g, name):
        pos = np.stack([node_index(g, a) for a in range(g.dim)], axis=1)
        period = np.where(g.periodic, g.domain_cells, 0)
        patch_of = np.repeat(np.arange(g.npatch), g.size)
        for axis in range(g.dim):
            t = g.faces(axis)
            expected, _ = brute_faces(g, axis)
            got = dict(zip(t.face.tolist(), t.sign.tolist()))
            self.assertEqual(got, expected, f"{name} L{g.level} axis {axis}")
            e = np.eye(g.dim, dtype=np.int64)[axis]
            # inner: same patch, one step inside
            np.testing.assert_array_equal(patch_of[t.inner], patch_of[t.face])
            np.testing.assert_array_equal(patch_of[t.face], t.patch)
            np.testing.assert_array_equal(pos[t.inner], pos[t.face] - t.sign[:, None] * e)
            # outer: one step outside, modulo the period on periodic axes
            diff = pos[t.outer] - (pos[t.face] + t.sign[:, None] * e)
            ok = np.where(period > 0, diff % np.maximum(period, 1) == 0, diff == 0)
            self.assertTrue(ok.all(), f"{name} outer position axis {axis}")
            self.assertTrue((patch_of[t.outer] != patch_of[t.face]).any() or g.npatch == 1)
            np.testing.assert_array_equal(t.plus, np.where(t.sign > 0, t.outer, t.inner))
            # one-sided flags agree with the listed nodes
            lo, hi = g.end_nodes(axis)
            listed = np.zeros(g.nnodes, bool)
            listed[t.face] = True
            for side_nodes, flags in ((lo, t.one_sided_lo), (hi, t.one_sided_hi)):
                missing = np.bincount(patch_of[side_nodes[~listed[side_nodes]]],
                                      minlength=g.npatch)
                np.testing.assert_array_equal(flags, missing > 0)

    def test_random_layouts(self):
        for name, hier in layouts():
            for g in hier.levels:
                self.check(g, name)

    def test_own_periodic_neighbour(self):
        # one patch spanning a periodic axis, and a 2-node-wide patch beside a wide one
        g = LevelGeometry.from_file_order(0, [[0, 0], [0, 3]], [[9, 2], [9, 3]], [0.1, 0.1],
                                          [10, 4], [True, False])
        self.check(g, "self-neighbour")
        t = g.faces(0)
        self.assertEqual(len(t), 2 * (4 + 2))  # every x-face node of both patches
        self.assertTrue((g.size > 0).all())
        t1 = g.faces(1)
        self.assertTrue(t1.one_sided_lo[g.store_index[0]])
        g1 = LevelGeometry.from_file_order(0, [[0]], [[0]], [1.0], [1], [True])  # 2 nodes
        self.check(g1, "single 2-node periodic patch")
        self.assertEqual(len(g1.faces(0)), 2)

    def assert_same_as_reference(self, g, name):
        for axis in range(g.dim):
            t = build_faces(g, axis)
            for k, v in reference_faces(g, axis).items():
                got = getattr(t, k)
                self.assertEqual(got.dtype, v.dtype, f"{name} axis {axis} {k}")
                np.testing.assert_array_equal(got, v, f"{name} axis {axis} {k}")

    def test_bit_identical_to_reference(self):
        for name, hier in layouts():
            for g in hier.levels:
                self.assert_same_as_reference(g, f"{name} L{g.level}")
        for dim in (1, 2, 3):  # many patches: shifts, runs and duplicates in bulk
            for periodic in (True, False):
                hier = make_large_amr(dim, 3000, nlevels=2, width=5, periodic=periodic, seed=dim)
                for g in hier.levels:
                    self.assert_same_as_reference(g, f"large {dim}D periodic={periodic}")

    def test_inner_derived_equals_stored_reference(self):
        """F4: ``inner`` (and blocks of it) recomputed from face, sign, stride equals the
        reference table bit for bit, on every axis, 1D/2D/3D, with 2-node patches."""
        for dim in (1, 2, 3):
            hier = make_large_amr(dim, 2000, nlevels=2, width=4, periodic=True, seed=dim + 7)
            for g in list(hier.levels) + [make_amr(dim, 1, seed=3, two_node=True).levels[0]]:
                for axis in range(dim):
                    t = build_faces(g, axis)
                    ref = reference_faces(g, axis)["inner"]
                    np.testing.assert_array_equal(t.inner, ref)
                    self.assertEqual(t.inner.dtype, np.int64)
                    m = len(t)
                    for a, b in ((0, m), (0, 1), (m // 3, 2 * m // 3 + 1), (m - 1, m), (5, 5)):
                        np.testing.assert_array_equal(t.inner_rows(a, b), ref[a:b])
        empty = LevelGeometry.from_file_order(0, np.zeros((0, 2)), np.zeros((0, 2)), [1, 1],
                                              [4, 4], [False, False])
        t = build_faces(empty, 1)
        self.assertEqual((len(t.inner), len(t.patch), t.inner.dtype), (0, 0, np.int64))


class TestCover(unittest.TestCase):
    def test_pieces_equal_brute_force_masks(self):
        for name, hier in layouts():
            for coarse in hier.levels:
                t = hier.covered(coarse.level)
                if not hier.has_level(coarse.level + 1):
                    self.assertEqual(len(t), 0)
                    continue
                fine = hier.level(coarse.level + 1)
                self.assertGreater(len(t), 0)
                for slot in range(coarse.npatch):
                    ext = coarse.upper[slot] - coarse.lower[slot] + 1
                    painted = np.zeros(ext, np.int64)
                    for k in np.nonzero(t.slot == slot)[0]:
                        sl = tuple(slice(t.lo[k, a], t.hi[k, a] + 1) for a in range(hier.dim))
                        painted[sl] += 1
                    self.assertLessEqual(painted.max(initial=0), 1, "pieces overlap")
                    # brute: a coarse cell is covered if its 2^d fine cells are in fine patches
                    cells = np.stack(np.meshgrid(*[np.arange(e) for e in ext], indexing="ij"),
                                     -1).reshape(-1, hier.dim) + coarse.lower[slot]
                    covered = np.zeros(len(cells), bool)
                    for flo, fhi in zip(fine.lower, fine.upper):
                        covered |= ((2 * cells >= flo) & (2 * cells + 1 <= fhi)).all(axis=1)
                    np.testing.assert_array_equal(painted.reshape(-1, order="C") > 0,
                                                  covered, f"{name} L{coarse.level} {slot}")
                np.testing.assert_array_equal(np.sort(t.slot), t.slot)
                self.assertTrue((t.fine_slot < fine.npatch).all())


class TestQueriesAndIdentity(unittest.TestCase):
    def test_patches_touching(self):
        g = make_amr(2, 1, seed=8).levels[0]
        rng = np.random.default_rng(0)
        for _ in range(20):
            lo = rng.integers(0, g.domain_cells)
            hi = lo + rng.integers(0, 6, size=2)
            ref = np.nonzero(((g.upper >= lo) & (g.lower <= hi)).all(axis=1))[0]
            np.testing.assert_array_equal(g.patches_touching(lo, hi), ref)

    def test_key_hash_pickle(self):
        hier = make_amr(2, 2, seed=9)
        g = hier.levels[1]
        flo, fhi = g.file_order_boxes()
        same = LevelGeometry.from_file_order(1, flo, fhi, g.dl, g.domain_cells, g.periodic)
        self.assertEqual(same.key(), g.key())
        self.assertEqual(same, g)
        self.assertEqual(hash(same), hash(g))
        perm = np.roll(np.arange(len(flo)), 1)
        other = LevelGeometry.from_file_order(1, flo[perm], fhi[perm], g.dl, g.domain_cells,
                                              g.periodic)
        self.assertNotEqual(other.key(), g.key())
        back = pickle.loads(pickle.dumps(g))
        self.assertEqual(back, g)
        np.testing.assert_array_equal(back.file_index, g.file_index)

    def test_hier(self):
        hier = make_amr(2, 3, seed=1, extent=dict(x=(1, 4)))
        self.assertEqual(hier.extent.lo[0], 1.0)
        self.assertEqual(hier.level_numbers(), (0, 1, 2))
        with self.assertRaises(KeyError):
            hier.level(5)
        with self.assertRaises(ValueError):
            HierGeometry(hier.levels[1:])
        np.testing.assert_allclose(hier.domain, hier.levels[0].domain_cells * 0.2)
        self.assertIs(hier.with_extent(None).covered(0), hier.covered(0))


class TestGeometryCache(unittest.TestCase):
    def boxes(self, seed):
        g = make_amr(2, 1, seed=seed).levels[0]
        lo, hi = g.file_order_boxes()
        return (0, lo, hi, g.dl, g.domain_cells, g.periodic)

    def test_reuse_and_lru(self):
        a, b, c = self.boxes(1), self.boxes(2), self.boxes(3)
        ga = LevelGeometry.from_file_order(*a)
        cache = GeometryCache(max_bytes=int(2.5 * ga.nbytes()))
        g1 = cache.get(*a)
        self.assertIs(cache.get(*a), g1)
        cache.get(*b)
        cache.get(*a)  # a is now the most recent
        cache.get(*c)  # evicts b
        self.assertIn(g1.key(), cache)
        self.assertNotIn(LevelGeometry.from_file_order(*b).key(), cache)
        self.assertLessEqual(cache.nbytes(), cache.max_bytes)

    def test_level_too_big_for_cap_stays_uncached(self):
        # cap fits the geometry and ONE face table, not all of them: faces() returns tables
        # without keeping any, and evicts nothing else (no thrash between axes)
        a, b = self.boxes(6), self.boxes(7)
        probe = LevelGeometry.from_file_order(*a)
        sizes = [build_faces(probe, axis).nbytes for axis in (0, 1)]
        other = LevelGeometry.from_file_order(*b)
        cap = probe.nbytes() + other.nbytes() + int(1.5 * min(sizes))
        cache = GeometryCache(max_bytes=cap)
        gb = cache.get(*b)
        g = cache.get(*a)
        for _ in range(2):
            for axis in (0, 1):
                t = g.faces(axis)
                self.assertGreater(len(t), 0)
                np.testing.assert_array_equal(t.face, build_faces(probe, axis).face)
        self.assertEqual(g._faces, {})
        self.assertIsNot(g.faces(0), g.faces(0))
        self.assertIn(gb.key(), cache)  # nothing evicted to make room
        self.assertLessEqual(cache.nbytes(), cache.max_bytes)
        # a cap holding the geometry and both tables keeps them
        big = GeometryCache(max_bytes=probe.nbytes() + 2 * max(sizes))
        g2 = big.get(*a)
        self.assertIs(g2.faces(0), g2.faces(0))
        self.assertIs(g2.faces(1), g2.faces(1))

    def test_eviction_drops_face_tables(self):
        """F4: an evicted geometry holds no face table, and does not keep new ones."""
        a, b = self.boxes(8), self.boxes(9)
        probe = LevelGeometry.from_file_order(*a)
        sizes = [build_faces(probe, axis).nbytes for axis in (0, 1)]
        # room for g and its two tables, not for another geometry beside them
        cache = GeometryCache(max_bytes=probe.nbytes() + 2 * max(sizes) + 64)
        self.assertLess(2 * max(sizes) - sum(sizes) + 64,
                        LevelGeometry.from_file_order(*b).nbytes())
        g = cache.get(*a)
        for axis in (0, 1):
            g.faces(axis)
        self.assertEqual(set(g._faces), {0, 1})
        self.assertEqual(cache.nbytes(), g.nbytes())
        g2 = cache.get(*b)  # evicts g (least recently used)
        self.assertNotIn(g.key(), cache)
        self.assertEqual(g._faces, {})
        t = g.faces(0)  # still usable, built again, not kept
        np.testing.assert_array_equal(t.face, build_faces(probe, 0).face)
        self.assertEqual(g._faces, {})
        self.assertIsNot(g.faces(0), t)
        self.assertEqual(cache.nbytes(), sum(x.nbytes() for x in (g2,)))
        self.assertLessEqual(cache.nbytes(), cache.max_bytes)
        self.assertIs(cache.add(g), g)  # back in a cache: tables kept again
        self.assertIs(g.faces(1), g.faces(1))
        cache.clear()
        self.assertEqual((g._faces, g2._faces), ({}, {}))
        g.faces(0)
        self.assertEqual(g._faces, {})

    def test_env(self):
        with mock.patch.dict(os.environ, {"PHARESEE_GEOM_CACHE_MB": "3"}):
            self.assertEqual(GeometryCache().max_bytes, 3 * 2**20)
        with mock.patch.dict(os.environ, {"PHARESEE_GEOM_CACHE_MB": "lots"}), \
                self.assertRaises(ValueError):
            GeometryCache()


if __name__ == "__main__":
    unittest.main()
