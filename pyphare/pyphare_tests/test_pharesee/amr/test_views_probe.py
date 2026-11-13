#!/usr/bin/env python3
"""Tests of pharesee.amr.views.probe: point probes and line cuts."""

import itertools
import unittest

import numpy as np

from pyphare.pharesee.amr.store import LevelData, Snapshot
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr
from pyphare.pharesee.amr.views.probe import probe


def linear(*c):
    return 1.0 + sum((a + 1.5) * x for a, x in enumerate(c))


def multilinear(*c):
    prod = 1.0
    for x in c:
        prod = prod * x
    return linear(*c) + 0.7 * prod


def rough(*c):
    v = np.sin(37.0 * c[0] + 1.0) * 3.0
    for a, x in enumerate(c[1:]):
        v = v + np.cos((53.0 + 7 * a) * x * (1 + 0.3 * c[0]))
    return v


def snapshot(hier, f=linear, lazy=False):
    src = fill(hier, f)
    if not lazy:
        return Snapshot(hier, tuple(src)), None
    fb = FakeBackend(src)
    return Snapshot(hier, tuple(LevelData.from_backend(fb, g, 0.0) for g in hier.levels)), fb


def brute_value(snap, p):
    """Reference (loops, test only): interpolate in EVERY patch of the finest level whose
    closed node box holds p; returns the list of values (all must agree)."""
    for geom, ld in zip(reversed(snap.hier.levels), reversed(snap.levels)):
        xi = p / geom.dl
        xi = np.where(np.abs(xi - np.round(xi)) < 1e-9, np.round(xi), xi)
        vals = []
        for s in range(geom.npatch):
            lo, hi = geom.lower[s], geom.upper[s] + 1
            if not ((xi >= lo) & (xi <= hi)).all():
                continue
            v = ld.view(s, "value")
            u = xi - lo
            i0 = np.minimum(np.floor(u).astype(int), v.shape - np.ones(geom.dim, int) * 2)
            t = u - i0
            acc = 0.0
            for corner in itertools.product((0, 1), repeat=geom.dim):
                w = np.prod([t[a] if c else 1 - t[a] for a, c in enumerate(corner)])
                acc += w * v[tuple(i0 + np.array(corner))]
            vals.append(acc)
        if vals:
            return geom.level, vals
    raise AssertionError(f"no patch holds {p}")


CASES = [(1, 2, 1, True), (1, 3, 2, False), (2, 1, 3, False), (2, 2, 4, True),
         (2, 3, 5, False), (3, 2, 7, True), (3, 3, 8, False)]


class TestPoints(unittest.TestCase):
    def test_exact_on_multilinear_fields_random_points(self):
        rng = np.random.default_rng(0)
        for dim, nl, seed, per in CASES:
            hier = make_amr(dim, nl, seed=seed, periodic=per, two_node=True)
            snap, _ = snapshot(hier, multilinear)
            P = rng.uniform(0, 1, (5000, dim)) * hier.domain
            with self.subTest(dim=dim, nl=nl):
                v = probe(snap, points=P)
                self.assertEqual(v.shape, (5000,))
                self.assertEqual(v.dtype, np.float64)
                np.testing.assert_allclose(v, multilinear(*P.T), rtol=0, atol=1e-11)

    def node_points(self, hier, rng, n=60):
        """Points on nodes, faces and corners of random patches of every level."""
        pts = []
        for geom in hier.levels:
            slots = rng.integers(0, geom.npatch, n)
            lo, hi = geom.lower[slots], geom.upper[slots] + 1
            # a node of the patch, one of its corners, a point on one face
            node = lo + (rng.random(lo.shape) * (hi - lo + 1)).astype(int)
            corner = np.where(rng.random(lo.shape) < 0.5, lo, hi)
            face = lo + rng.random(lo.shape) * (hi - lo)
            ax = rng.integers(0, geom.dim, n)
            face[np.arange(n), ax] = np.where(rng.random((n, 1)) < 0.5, lo, hi)[np.arange(n), ax]
            pts += [node * geom.dl, corner * geom.dl, face * geom.dl]
        return np.concatenate(pts)

    def test_nodes_faces_coarse_fine_borders_against_brute_force(self):
        rng = np.random.default_rng(1)
        for dim, nl, seed, per in CASES:
            hier = make_amr(dim, nl, seed=seed, periodic=per, two_node=True)
            snap, _ = snapshot(hier, rough)
            P = self.node_points(hier, rng)
            v = probe(snap, points=P)
            for p, got in zip(P, v):
                _, vals = brute_value(snap, p)
                # every patch of the finest level holding p gives the same value...
                self.assertTrue(np.all(np.asarray(vals) == vals[0]), f"{p}: {vals}")
                # ...and so does probe(), bit for bit
                self.assertEqual(got, vals[0], f"{dim}D {p}")

    def test_value_comes_from_finest_level(self):
        """On a coarse-fine border (fine patch face) the value is the fine level's, which
        differs from the coarse value when the coarse data is not the restriction."""
        hier = make_amr(2, 2, seed=4)
        src = fill(hier, linear)
        src[1].buffer("value")[:] += 100.0  # make the fine level distinguishable
        snap = Snapshot(hier, tuple(src))
        g1 = hier.levels[1]
        s = 0
        lo, hi = g1.lower[s], g1.upper[s] + 1
        pts = np.array([lo, [lo[0], (lo[1] + hi[1]) / 2], hi]) * g1.dl
        np.testing.assert_allclose(probe(snap, points=pts), linear(*pts.T) + 100.0)
        # just outside the level-1 boxes: level 0
        lower, upper = g1.lower, g1.upper + 1
        out = []
        for p in np.random.default_rng(0).uniform(0, 1, (2000, 2)) * hier.domain:
            xi = p / g1.dl
            if not ((xi >= lower) & (xi <= upper)).all(axis=1).any():
                out.append(p)
        out = np.array(out)
        np.testing.assert_allclose(probe(snap, points=out), linear(*out.T), atol=1e-12)

    def test_coordinate_arguments(self):
        hier = make_amr(3, 2, seed=7)
        snap, _ = snapshot(hier, multilinear)
        xs = np.linspace(0.1, 3.0, 7)
        v = probe(snap, x=xs, y=1.3, z=np.full(7, 2.0))
        np.testing.assert_allclose(v, multilinear(xs, 1.3, 2.0), atol=1e-12)
        v = probe(snap, x=0.5, y=0.6, z=0.7)
        self.assertEqual(np.shape(v), ())
        self.assertAlmostEqual(float(v), multilinear(0.5, 0.6, 0.7), places=12)
        grid = np.meshgrid(xs, [1.0, 2.0], indexing="ij")
        v = probe(snap, x=grid[0], y=grid[1], z=0.4)
        self.assertEqual(v.shape, (7, 2))
        h1 = make_amr(1, 2, seed=1)
        s1, _ = snapshot(h1)
        np.testing.assert_allclose(probe(s1, x=[1.0, 2.5]), linear(np.array([1.0, 2.5])))
        np.testing.assert_allclose(probe(s1, points=[1.0, 2.5]), linear(np.array([1.0, 2.5])))

    def test_components(self):
        hier = make_amr(2, 2, seed=4)
        snap = Snapshot(hier, tuple(fill(hier, {"x": linear, "y": lambda x, y: 2 * y})))
        P = np.random.default_rng(2).uniform(0, 1, (50, 2)) * hier.domain
        v = probe(snap, points=P)
        self.assertEqual(set(v), {"x", "y"})
        np.testing.assert_allclose(v["y"], 2 * P[:, 1], atol=1e-12)
        c, v = probe(snap, y=1.1)
        np.testing.assert_allclose(v["x"], linear(c, 1.1), atol=1e-12)

    def test_periodic_wrap_and_domain_errors(self):
        hier = make_amr(2, 2, seed=4, periodic=(True, False))
        L = hier.domain
        # periodic in x like real periodic data (nodes x = 0 and x = L hold the same value)
        snap, _ = snapshot(hier, lambda x, y: np.sin(2 * np.pi * x / L[0]) * np.cos(3 * y) + y)
        P = np.array([[0.37, 1.1], [3.3, 0.2], [L[0], 2.0]])
        base = probe(snap, points=P)
        # x + L wraps back to x (to the rounding of x + L - L)
        np.testing.assert_allclose(probe(snap, points=P + [L[0], 0.0]), base, rtol=0, atol=1e-12)
        np.testing.assert_allclose(probe(snap, points=P[:2] - [2 * L[0], 0.0]), base[:2],
                                   rtol=0, atol=1e-12)
        # the domain end x = L is kept as is (same node as x = 0 in periodic data)
        np.testing.assert_allclose(base[2], probe(snap, x=0.0, y=2.0), rtol=0, atol=1e-12)
        with self.assertRaisesRegex(ValueError, r"outside the domain.*along y"):
            probe(snap, points=[[1.0, L[1] + 0.5]])
        with self.assertRaisesRegex(ValueError, "-0.5"):
            probe(snap, x=1.0, y=-0.5)

    def test_extent(self):
        hier = make_amr(2, 2, seed=4, extent=dict(x=(1.0, 3.0), y=(2.0, 4.0)))
        snap, _ = snapshot(hier)
        np.testing.assert_allclose(probe(snap, x=1.0, y=4.0), linear(1.0, 4.0))
        with self.assertRaisesRegex(ValueError, r"\[0\.5, 3\.0\].*extent"):
            probe(snap, points=[[1.5, 2.5], [0.5, 3.0]])

    def test_argument_errors(self):
        s1, _ = snapshot(make_amr(1, 1, seed=0))
        s2, _ = snapshot(make_amr(2, 1, seed=0))
        s3, _ = snapshot(make_amr(3, 1, seed=0))
        with self.assertRaises(TypeError):
            probe(s1)
        with self.assertRaises(TypeError):
            probe(s2)
        with self.assertRaises(TypeError):
            probe(s2, x=1.0, z=1.0)
        with self.assertRaises(TypeError):
            probe(s2, x=1.0, points=[[1.0, 1.0]])
        with self.assertRaises(TypeError):
            probe(s2, y=[1.0, 2.0])
        with self.assertRaisesRegex(NotImplementedError, "finest"):
            probe(s3, z=1.0)
        with self.assertRaises(ValueError):
            probe(s2, points=[[1.0, 2.0, 3.0]])
        with self.assertRaises(ValueError):
            probe(s2, x=np.nan, y=1.0)

    def test_reads_only_needed_patches(self):
        hier = make_amr(2, 3, seed=5, domain_cells=(96, 80), patch_cells=8)
        snap, fb = snapshot(hier, rough, lazy=True)
        rng = np.random.default_rng(3)
        pts = rng.uniform([5.0, 5.0], [6.0, 5.5], (40, 2))
        v = probe(snap, points=pts)
        ref, _ = snapshot(hier, rough)
        np.testing.assert_array_equal(v, probe(ref, points=pts))
        total = 0
        for geom in hier.levels:
            counts = fb.reads.get((0.0, geom.level))
            if counts is None:
                continue
            self.assertLessEqual(counts.max(), 1)
            read = geom.store_index[np.nonzero(counts)[0]]
            total += len(read)
            # every patch read holds at least one point
            for s in read:
                lo, hi = geom.lower[s] * geom.dl, (geom.upper[s] + 1) * geom.dl
                self.assertTrue(((pts >= lo) & (pts <= hi)).all(axis=1).any())
        self.assertLessEqual(fb.calls, len(hier.levels))
        self.assertLess(total, 25)


class TestLines(unittest.TestCase):
    def reference_coords(self, hier, fixed, free):
        """Every finest-grid position whose finest holding level (among patches crossing the
        line) has a node there (loops, test only)."""
        top = hier.levels[-1]
        Ntop = top.domain_cells[free]
        keep = []
        for P in range(Ntop + 1):
            for geom in reversed(hier.levels):
                r = 2 ** (top.level - geom.level)
                held = False
                for s in range(geom.npatch):
                    ok = geom.lower[s, free] * r <= P <= (geom.upper[s, free] + 1) * r
                    for a, v in fixed.items():
                        xi = v / geom.dl[a]
                        ok &= geom.lower[s, a] <= xi + 1e-9 and xi - 1e-9 <= geom.upper[s, a] + 1
                    held |= ok
                if held:
                    if P % r == 0:
                        keep.append(P)
                    break
        return np.array(keep) * top.dl[free]

    def test_line_coords_and_values(self):
        rng = np.random.default_rng(4)
        for dim, nl, seed, per in CASES:
            if dim == 1:
                continue
            hier = make_amr(dim, nl, seed=seed, periodic=per)
            snap, _ = snapshot(hier, multilinear)
            for free in range(dim):
                fixed = {a: float(rng.uniform(0.2, 0.8) * hier.domain[a])
                         for a in range(dim) if a != free}
                if free == 0 and dim == 2:  # a fixed coordinate exactly on a fine node
                    fixed[1] = float(hier.levels[-1].lower[0, 1] * hier.levels[-1].dl[1])
                kw = {"xyz"[a]: v for a, v in fixed.items()}
                with self.subTest(dim=dim, nl=nl, free=free, fixed=fixed):
                    c, v = probe(snap, **kw)
                    self.assertTrue((np.diff(c) > 0).all(), "sorted, unique")
                    np.testing.assert_allclose(c, self.reference_coords(hier, fixed, free),
                                               rtol=0, atol=1e-12)
                    steps = np.unique(np.round(np.diff(c) / hier.levels[-1].dl[free]))
                    allowed = {2 ** (hier.finest_level - g.level) for g in hier.levels}
                    self.assertTrue(set(steps.astype(int)) <= allowed, steps)
                    pts = [c if a == free else np.full_like(c, fixed[a]) for a in range(dim)]
                    np.testing.assert_allclose(v, multilinear(*pts), rtol=0, atol=1e-11)

    def test_line_extent_and_laziness(self):
        hier = make_amr(2, 2, seed=4, extent=dict(x=(1.0, 3.0)), domain_cells=(96, 80))
        snap, fb = snapshot(hier, lazy=True)
        c, v = probe(snap, y=5.03)
        self.assertGreaterEqual(c[0], 1.0)
        self.assertLessEqual(c[-1], 3.0)
        np.testing.assert_allclose(v, linear(c, 5.03), atol=1e-12)
        # only patches crossing the line at x in [1, 3] were read
        for geom in hier.levels:
            counts = fb.reads.get((0.0, geom.level))
            if counts is None:
                continue
            for s in geom.store_index[np.nonzero(counts)[0]]:
                self.assertLessEqual(geom.lower[s, 1] * geom.dl[1], 5.03)
                self.assertGreaterEqual((geom.upper[s, 1] + 1) * geom.dl[1], 5.03)
                self.assertLessEqual(geom.lower[s, 0] * geom.dl[0], 3.0)
                self.assertGreaterEqual((geom.upper[s, 0] + 1) * geom.dl[0], 1.0)


if __name__ == "__main__":
    unittest.main()
