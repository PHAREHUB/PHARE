#!/usr/bin/env python3
"""Tests of pharesee.amr.boxes: BinIndex, flat indices, PhysBox."""

import math
import unittest

import numpy as np

from pyphare.pharesee.amr.boxes import (
    BinIndex,
    PhysBox,
    box_flat_indices,
    coarsen,
    coarsen_inner,
    intersect_tables,
    iter_box_flat_indices,
    periodic_shifts,
    refine,
)


def random_boxes(rng, n, dim, span, wmax):
    lo = rng.integers(0, span, size=(n, dim))
    hi = lo + rng.integers(0, wmax, size=(n, dim))
    return lo.astype(np.int64), hi.astype(np.int64)


def brute_pairs(alo, ahi, blo, bhi, shifts):
    out = set()
    for k, s in enumerate(shifts):
        for i in range(len(alo)):
            lo = np.maximum(alo[i] + s, blo)
            hi = np.minimum(ahi[i] + s, bhi)
            for j in np.nonzero((lo <= hi).all(axis=1))[0]:
                out.add((i, int(j), k, tuple(lo[j]), tuple(hi[j])))
    return out


def as_set(p):
    return {
        (int(a), int(b), int(s), tuple(lo), tuple(hi))
        for a, b, s, lo, hi in zip(p.a, p.b, p.shift, p.lo, p.hi)
    }


class TestBinIndex(unittest.TestCase):
    def test_pairs_equal_brute_force(self):
        rng = np.random.default_rng(1)
        for dim in (1, 2, 3):
            for trial in range(4):
                blo, bhi = random_boxes(rng, 60, dim, 40, 9)
                alo, ahi = random_boxes(rng, 50, dim, 40, 12)
                if trial == 3:  # one huge query box takes the scan path
                    ahi[0] = alo[0] + 10**4
                idx = BinIndex(blo, bhi, width=2 if trial == 2 else None)
                if trial == 3:
                    idx.big_query_bins = 8
                got = idx.pairs(alo, ahi)
                ref = brute_pairs(alo, ahi, blo, bhi, np.zeros((1, dim), np.int64))
                self.assertEqual(as_set(got), ref, (dim, trial))
                self.assertEqual(len(got.a), len(ref), "a pair was reported twice")

    def test_sparse_and_dense_bins_unsorted(self):
        # dense bin-start table (compact layouts) and binary search (sparse ones) agree with
        # brute force; sort=False gives the same rows in another order
        rng = np.random.default_rng(3)
        for dim in (1, 2, 3):
            n = np.array([40, 30, 20][:dim])
            lo, hi = random_boxes(rng, 50, dim, 30, 8)
            hi = np.minimum(hi, n - 1)
            lo = np.minimum(lo, hi)
            shifts = periodic_shifts(n, [True] * dim)
            ref = brute_pairs(lo - 1, hi + 1, lo, hi, shifts)
            dense = BinIndex(lo, hi)
            self.assertIsNotNone(dense.starts)
            far = 10**7 if dim == 1 else 10**5  # huge bounding box: sparse
            far_lo = np.concatenate([lo, [[far] * dim]])
            far_hi = np.concatenate([hi, [[far] * dim]])
            sparse = BinIndex(far_lo, far_hi, width=1)
            self.assertIsNone(sparse.starts)
            for idx in (dense, sparse):
                got = idx.pairs(lo - 1, hi + 1, shifts, sort=False)
                self.assertEqual(as_set(got), ref)
                self.assertEqual(len(got.a), len(ref))
                srt = idx.pairs(lo - 1, hi + 1, shifts)
                key = np.lexsort((srt.shift, srt.b, srt.a))
                np.testing.assert_array_equal(key, np.arange(len(key)))

    def test_query_and_empty(self):
        lo = np.array([[0, 0], [5, 5], [10, 0]])
        hi = np.array([[4, 4], [9, 9], [12, 3]])
        idx = BinIndex(lo, hi)
        np.testing.assert_array_equal(idx.query([4, 4], [5, 5]), [0, 1])
        np.testing.assert_array_equal(idx.query([13, 13], [20, 20]), [])
        np.testing.assert_array_equal(idx.query([-100, -100], [100, 100]), [0, 1, 2])
        empty = BinIndex(np.zeros((0, 2)), np.zeros((0, 2)))
        self.assertEqual(len(empty.pairs(lo, hi).a), 0)
        self.assertEqual(len(intersect_tables(lo, hi, lo, hi).a), 3)

    def test_rejects_empty_box(self):
        with self.assertRaises(ValueError):
            BinIndex([[3]], [[2]])


class TestFlatIndices(unittest.TestCase):
    def reference(self, offsets, shapes, lo, size):
        out = []
        for o, s, l, z in zip(offsets, shapes, lo, size):
            a = np.arange(np.prod(s)).reshape(s, order="F")
            sl = tuple(slice(l[i], l[i] + z[i]) for i in range(len(s)))
            out.append(o + a[sl].ravel(order="F"))
        return np.concatenate(out) if out else np.zeros(0, np.int64)

    def random_pieces(self, rng, dim, n):
        shapes = rng.integers(2, 7, size=(n, dim))
        lo = (rng.random((n, dim)) * shapes).astype(np.int64)
        size = (rng.random((n, dim)) * (shapes - lo + 1)).astype(np.int64)
        offsets = np.concatenate([[0], np.cumsum(np.prod(shapes, 1))[:-1]])
        return offsets, shapes, lo, size

    def test_equals_reshape_slices(self):
        rng = np.random.default_rng(3)
        for dim in (1, 2, 3):
            pieces = self.random_pieces(rng, dim, 40)
            np.testing.assert_array_equal(box_flat_indices(*pieces), self.reference(*pieces))

    def test_run_axis_same_set_per_piece(self):
        rng = np.random.default_rng(5)
        for dim in (2, 3):
            pieces = self.random_pieces(rng, dim, 30)
            counts = np.prod(pieces[3], axis=1)
            ref = np.split(box_flat_indices(*pieces), np.cumsum(counts)[:-1])
            for run in list(range(1, dim)) + ["auto"]:
                got = np.split(box_flat_indices(*pieces, run_axis=run), np.cumsum(counts)[:-1])
                for r, g in zip(ref, got):
                    np.testing.assert_array_equal(np.sort(r), np.sort(g))
            # run axis 1 of a 2D sub-box: y fastest
            np.testing.assert_array_equal(
                box_flat_indices([0], [[4, 3]], [[1, 1]], [[2, 2]], run_axis=1), [5, 9, 6, 10]
            )

    def test_blocks(self):
        rng = np.random.default_rng(4)
        pieces = self.random_pieces(rng, 3, 50)
        whole = box_flat_indices(*pieces)
        counts = np.prod(pieces[3], axis=1)
        blocks = list(iter_box_flat_indices(*pieces, max_nodes=37))
        self.assertGreaterEqual(len(blocks), max(2, counts.sum() // 37))
        np.testing.assert_array_equal(np.concatenate([b[2] for b in blocks]), whole)
        for i0, i1, idx in blocks:
            self.assertTrue(len(idx) <= 37 or i1 == i0 + 1)
            self.assertEqual(len(idx), counts[i0:i1].sum())
        self.assertEqual(blocks[-1][1], 50)


class TestCoarsenRefine(unittest.TestCase):
    def test_values(self):
        lo, hi = coarsen(np.array([[4], [5]]), np.array([[9], [8]]))
        np.testing.assert_array_equal(lo, [[2], [2]])
        np.testing.assert_array_equal(hi, [[4], [4]])
        lo, hi = coarsen_inner(np.array([4, 5, 5]), np.array([9, 8, 5]))
        np.testing.assert_array_equal(lo, [2, 3, 3])
        np.testing.assert_array_equal(hi, [4, 3, 2])  # last is empty
        lo, hi = refine(np.array([2]), np.array([4]))
        np.testing.assert_array_equal((lo, hi), ([4], [9]))


class TestPhysBox(unittest.TestCase):
    def test_parse(self):
        b = PhysBox.parse(dict(x=(1.0, 2.0)), 2)
        self.assertEqual(b.lo, (1.0, -math.inf))
        self.assertEqual(b.hi, (2.0, math.inf))
        self.assertIsNone(PhysBox.parse(None, 2))
        self.assertIs(PhysBox.parse(b, 2), b)
        self.assertEqual(PhysBox.parse(dict(y=(None, 3)), 2).hi, (math.inf, 3.0))

    def test_rejects_old_syntax(self):
        class Box:  # stands for pyphare.core.box.Box
            lower, upper = (0, 0), (10, 10)

        for bad in (Box(), (0, 10), [(0, 1), (0, 1)], "x"):
            with self.assertRaises(TypeError) as cm:
                PhysBox.parse(bad, 2)
            self.assertIn("box=dict(x=(x0, x1)", str(cm.exception))
        with self.assertRaises(TypeError):
            PhysBox.parse(dict(x=5), 2)
        with self.assertRaises(ValueError):
            PhysBox.parse(dict(z=(0, 1)), 2)
        with self.assertRaises(ValueError):
            PhysBox.parse(dict(x=(2, 1)), 1)

    def test_cells_rule(self):
        dl, n = [0.5], [10]
        cells = lambda x0, x1: tuple(int(v[0]) for v in PhysBox((x0,), (x1,)).cells(dl, n))
        self.assertEqual(cells(1.2, 2.0), (2, 3))
        self.assertEqual(cells(1.0, 1.0), (1, 2))  # plane on a face: both cells
        self.assertEqual(cells(1.2, 1.2), (2, 2))
        self.assertEqual(cells(0.0, 0.0), (0, 0))
        self.assertEqual(cells(5.0, 5.0), (9, 9))
        self.assertEqual(cells(-3.0, 100.0), (0, 9))
        # 0.6 / 0.2 == 2.9999999999999996: must still start at cell 3
        lo, hi = PhysBox((0.6,), (1.0,)).cells([0.2], [10])
        self.assertEqual((int(lo[0]), int(hi[0])), (3, 4))
        with self.assertRaises(ValueError):
            PhysBox((7.0,), (8.0,)).cells(dl, n)

    def test_selects(self):
        box = PhysBox((1.2,), (2.0,))
        lower = np.array([[0], [4], [2]])
        upper = np.array([[1], [7], [3]])
        # cells [2, 3]; patch [4, 7] starts at x = 2.0: touches only, not selected
        np.testing.assert_array_equal(box.selects(lower, upper, [0.5], [10]), [False, False, True])

    def test_intersect(self):
        a = PhysBox.parse(dict(x=(0, 2), y=(0, 2)), 2)
        b = PhysBox.parse(dict(x=(1, 3)), 2)
        self.assertEqual(a.intersect(b), PhysBox((1.0, 0.0), (2.0, 2.0)))
        self.assertIsNone(a.intersect(PhysBox((3.0, 0.0), (4.0, 1.0))))
        self.assertIs(a.intersect(None), a)
        self.assertEqual(hash(a), hash(PhysBox((0, 0), (2, 2))))


if __name__ == "__main__":
    unittest.main()
