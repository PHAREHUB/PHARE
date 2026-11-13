#!/usr/bin/env python3
"""Tests of pharesee.amr.ops.reduce: composite-grid integrals, means, min/max."""

import itertools
import multiprocessing
import unittest
import warnings
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.boxes import PhysBox
from pyphare.pharesee.amr.ops import arith, reduce
from pyphare.pharesee.amr.store import LevelData, Snapshot
from pyphare.pharesee.amr.testing import fill, make_amr, make_large_amr


def snapshot(hier, funcs, dtype=np.float64):
    return Snapshot(hier, tuple(fill(hier, funcs, dtype=dtype)))


def random_snapshot(hier, names=("value",), seed=0, dtype=np.float64):
    rng = np.random.default_rng(seed)
    return Snapshot(hier, tuple(
        LevelData.from_arrays(g, 0.0, {n: rng.standard_normal(g.nnodes).astype(dtype)
                                       for n in names}) for g in hier.levels))


def layouts():
    out = []
    for dim in (1, 2, 3):
        out.append(make_amr(dim, 3, seed=dim))
        out.append(make_amr(dim, 2, seed=10 + dim, two_node=True, singletons=True,
                            periodic=True))
        out.append(make_amr(dim, 1, seed=20 + dim, two_node=True))
    return out


def extents(hier):
    """Boxes cutting cells on both sides, narrower than a cell, cutting a covered piece."""
    dim, L, dl = hier.dim, hier.domain, hier.dl0
    out = [None]
    out.append(PhysBox(tuple(0.13 * L + 0.37 * dl), tuple(0.71 * L - 0.21 * dl)))
    out.append(PhysBox(tuple(0.5 * L + 0.1 * dl), tuple(0.5 * L + 0.45 * dl)))  # < 1 cell
    if len(hier.levels) > 1:
        g = hier.levels[1]
        k = int(np.argmax(g.size))
        mid = (g.lower[k] + g.upper[k] + 1) / 2 * g.dl
        out.append(PhysBox(tuple(mid + 0.3 * g.dl), tuple(np.minimum(mid + 0.4 * L, L))))
    out.append(PhysBox(tuple(-np.ones(dim)), tuple(0.3 * L)))  # partly outside the domain
    return out


def covered_cells(hier, ilvl):
    """Dense bool grid of level ilvl's cells under level ilvl+1 (brute force)."""
    g = hier.level(ilvl)
    grid = np.zeros(tuple(g.domain_cells), dtype=bool)
    if not hier.has_level(ilvl + 1):
        return grid
    f = hier.level(ilvl + 1)
    for lo, hi in zip(f.lower, f.upper):
        clo, chi = -(-lo // 2), (hi + 1) // 2 - 1
        grid[tuple(slice(a, b + 1) for a, b in zip(clo, chi))] = True
    return grid


def ref_integral(snap):
    """Per-patch, per-cell reference: sum over uncovered cells of the cut multilinear rule."""
    hier = snap.hier
    lo, hi = (np.zeros(hier.dim), hier.domain) if hier.extent is None else \
        hier.extent.clipped(hier.domain)
    total = np.zeros(len(snap.names))
    for ld in snap.levels:
        g = ld.geom
        cov = covered_cells(hier, g.level)
        for slot in range(g.npatch):
            cells = [np.arange(g.lower[slot, a], g.upper[slot, a] + 1) for a in range(g.dim)]
            w = []  # per axis: (ncells, 2) corner weights
            for a in range(g.dim):
                s = np.clip(lo[a] / g.dl[a] - cells[a], 0, 1)
                e = np.clip(hi[a] / g.dl[a] - cells[a], 0, 1)
                e = np.maximum(e, s)
                right = (e ** 2 - s ** 2) / 2
                w.append(np.stack([(e - s) - right, right], axis=1) * g.dl[a])
            sel = np.ix_(*cells)
            alive = ~cov[sel]
            for k, n in enumerate(snap.names):
                v = ld.view(slot, n).astype(np.float64)
                for corner in itertools.product((0, 1), repeat=g.dim):
                    wc = np.ones(alive.shape)
                    for a in range(g.dim):
                        shape = [1] * g.dim
                        shape[a] = -1
                        wc = wc * w[a][:, corner[a]].reshape(shape)
                    vals = v[tuple(slice(c, c + len(cells[a])) for a, c in enumerate(corner))]
                    total[k] += (vals * wc)[alive].sum()
    return total


class TestIntegral(unittest.TestCase):
    def test_equals_reference(self):
        for hier in layouts():
            for ext in extents(hier):
                h = hier.with_extent(ext)
                snap = random_snapshot(h, ("x", "y", "z"), seed=hier.dim)
                got = reduce.integral(snap)
                ref = ref_integral(snap)
                np.testing.assert_allclose(got, ref, rtol=1e-12,
                                           atol=1e-12 * np.abs(ref).max(), err_msg=f"{ext}")

    def test_linear_exact_and_measure(self):
        for hier in layouts():
            dim = hier.dim
            c = np.array([0.7, -1.3, 2.1])[:dim]
            snap = snapshot(hier.with_extent(None), lambda *X, c=c, dim=dim: 1.0 + sum(
                c[a] * X[a] for a in range(dim)))
            for ext in extents(hier):
                h = hier.with_extent(ext)
                s = Snapshot(h, snap.levels)
                lo, hi = (np.zeros(dim), h.domain) if ext is None else ext.clipped(h.domain)
                vol = np.prod(hi - lo)
                exact = vol * (1.0 + (c * (lo + hi) / 2).sum())
                self.assertAlmostEqual(reduce.measure(h), vol, delta=1e-14 * vol)
                self.assertAlmostEqual(reduce.integral(s), exact, delta=1e-12 * abs(exact) + 1e-14,
                                       msg=f"{dim}d {ext}")
                ones = snapshot(h, lambda *X: 1.0 + 0 * X[0], dtype=np.float32)
                self.assertAlmostEqual(reduce.integral(ones), vol, delta=1e-12 * vol)
                self.assertAlmostEqual(reduce.mean(s), exact / vol, delta=1e-12 * abs(exact / vol)
                                       + 1e-14)

    def test_coarse_data_under_fine_does_not_contribute(self):
        for hier in layouts():
            if len(hier.levels) < 2:
                continue
            for ext in extents(hier):
                h = hier.with_extent(ext)
                snap = random_snapshot(h, seed=4)
                before = reduce.integral(snap)
                for ld in snap.levels[:-1]:
                    idx = h.covered(ld.level).node_indices(ld.geom, closed=False)
                    buf = ld.buffer("value")
                    buf[idx] = np.where(np.arange(len(idx)) % 2, np.nan, 1e30)
                self.assertEqual(reduce.integral(snap), before)
                self.assertTrue(np.isfinite(reduce.integral(snap)))

    def test_vector_and_blocks(self):
        hier = make_amr(3, 3, seed=9, two_node=True)
        h = hier.with_extent(extents(hier)[3])
        snap = random_snapshot(h, ("x", "y", "z"), seed=5)
        ref = reduce.integral(snap)
        self.assertEqual(ref.shape, (3,))
        with mock.patch.object(arith, "BLOCK_NODES", 41):
            np.testing.assert_allclose(reduce.integral(snap), ref, rtol=1e-13)
            np.testing.assert_array_equal(reduce.maximum(snap), reduce.maximum(snap))
        self.assertIsInstance(reduce.integral(arith.component(snap, "x")), float)

    def test_zero_measure_mean_raises(self):
        hier = make_amr(2, 1, seed=1)
        h = hier.with_extent(PhysBox((1.0, 0.0), (1.0, 5.0)))
        s = random_snapshot(h)
        self.assertEqual(reduce.measure(h), 0.0)
        self.assertEqual(reduce.integral(s), 0.0)
        with self.assertRaises(ValueError):
            reduce.mean(s)
        with self.assertRaises(TypeError):
            reduce.integral(3.0)


class TestMinMax(unittest.TestCase):
    def reference(self, snap, how):
        lo, hi = (np.zeros(snap.hier.dim), snap.hier.domain) if snap.hier.extent is None else \
            snap.hier.extent.clipped(snap.hier.domain)
        vals = []
        for ld in snap.levels:
            g = ld.geom
            keep = np.ones(g.nnodes, dtype=bool)
            for a in range(g.dim):
                x = g.node_coords(a)
                keep &= (x >= lo[a] - 1e-9 * g.dl[a]) & (x <= hi[a] + 1e-9 * g.dl[a])
            keep[snap.hier.covered(ld.level).node_indices(g, closed=True)] = False
            vals.append(ld.buffer("value")[keep])
        v = np.concatenate(vals)
        if len(v) == 0:
            return None
        return v.min() if how == "min" else v.max()

    def test_ignores_covered_nodes_and_respects_extent(self):
        for hier in layouts():
            for ext in extents(hier):
                h = hier.with_extent(ext)
                snap = random_snapshot(h, seed=6)
                for ld in snap.levels[:-1]:  # extreme values on closed covered coarse nodes
                    idx = h.covered(ld.level).node_indices(ld.geom, closed=True)
                    ld.buffer("value")[idx[::2]] = 1e9
                    ld.buffer("value")[idx[1::2]] = -1e9
                if self.reference(snap, "max") is None:  # box between two node planes
                    with self.assertRaises(ValueError):
                        reduce.maximum(snap)
                    continue
                self.assertEqual(reduce.maximum(snap), self.reference(snap, "max"))
                self.assertEqual(reduce.minimum(snap), self.reference(snap, "min"))
                self.assertLess(reduce.maximum(snap), 1e8)
                self.assertGreater(reduce.minimum(snap), -1e8)

    def test_per_level(self):
        hier = make_amr(2, 3, seed=3)
        snap = random_snapshot(hier, ("x", "y", "z"), seed=7, dtype=np.float32)
        for how, func in (("mean", np.mean), ("sum", np.sum), ("min", np.min), ("max", np.max)):
            out = reduce.per_level(snap, how)
            self.assertEqual(sorted(out), [0, 1, 2])
            for ld in snap.levels:
                expect = [func(ld.buffer(c).astype(np.float64)) for c in "xyz"]
                np.testing.assert_allclose(out[ld.level], expect, rtol=1e-12)
        scalar = reduce.per_level(arith.component(snap, "x"), "max")
        self.assertIsInstance(scalar[0], float)
        with self.assertRaises(ValueError):
            reduce.per_level(snap, "median")


class TestGroupPath(unittest.TestCase):
    """F6: whole uncovered patches are summed per shape group with cached 1-D weights; the
    result equals the block path (every patch through per-node weights)."""

    def test_group_path_equals_block_path(self):
        cases = layouts() + [make_large_amr(2, 3000, nlevels=2, width=6, seed=3),
                             make_large_amr(3, 500, nlevels=2, width=4, seed=4, periodic=True)]
        used = 0
        for hier in cases:
            for ext in extents(hier):
                h = hier.with_extent(ext)
                snap = random_snapshot(h, ("x", "y"), seed=hier.dim, dtype=np.float32)
                with mock.patch.object(reduce, "_clean_integral",
                                       wraps=reduce._clean_integral) as spy:
                    got = reduce.integral(snap)
                used += spy.call_count
                with mock.patch.object(reduce, "GROUP_LOOP_MAX", -1):
                    ref = reduce.integral(snap)
                np.testing.assert_allclose(got, ref, rtol=1e-13, atol=1e-13 * np.abs(ref).max(),
                                           err_msg=f"{hier} {ext}")
        self.assertGreater(used, 0)


_FORK_SNAPSHOT = {}


def _child_mean(_):
    return reduce.mean(_FORK_SNAPSHOT["snap"])


@unittest.skipUnless("fork" in multiprocessing.get_all_start_methods(), "no fork on this OS")
class TestForkSafety(unittest.TestCase):
    """A5: a reduction run in a forked child after the parent used numpy (and BLAS) returns.

    A threaded BLAS (OpenBLAS/MKL) whose threads were started before fork() deadlocks in
    the child; reduce.py must not call BLAS at all."""

    def test_mean_in_fork_child(self):
        # big enough (~0.8 M nodes) for a threaded BLAS to use its threads on a dot product
        hier = make_large_amr(2, 10000, nlevels=2, width=8, seed=5)
        snap = random_snapshot(hier, seed=1)
        _FORK_SNAPSHOT["snap"] = snap
        a = np.ones((300, 300))
        _ = a @ a  # start the BLAS thread pool in the parent, as a real session does
        parent = reduce.mean(snap)
        with warnings.catch_warnings():
            # Python >= 3.12 warns that forking a multi-threaded process may deadlock:
            # that is exactly what this test checks does not happen
            warnings.simplefilter("ignore", DeprecationWarning)
            pool = multiprocessing.get_context("fork").Pool(1)
        with pool:
            res = pool.map_async(_child_mean, [0])
            try:
                child = res.get(timeout=20)[0]
            except multiprocessing.TimeoutError:
                pool.terminate()
                self.fail("reduce.mean hung in a fork child (BLAS after fork?)")
        self.assertEqual(child, parent)
        _FORK_SNAPSHOT.clear()


if __name__ == "__main__":
    unittest.main()
