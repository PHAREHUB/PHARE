#!/usr/bin/env python3
"""Tests of pharesee.amr.ops.arith: elementwise operations between AMR snapshots."""

import unittest
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.boxes import PhysBox
from pyphare.pharesee.amr.geometry import HierGeometry
from pyphare.pharesee.amr.ops import arith, reduce
from pyphare.pharesee.amr.store import Snapshot
from pyphare.pharesee.amr.testing import make_amr
from pyphare_tests.test_pharesee.amr._common import (
    assert_all_equal,
    comps,
    random_snapshot,
    snapshot,
    vec,
)

XYZ = ("x", "y", "z")


def layouts():
    out = []
    for dim in (1, 2, 3):
        out.append(make_amr(dim, 3, seed=dim))
        out.append(make_amr(dim, 2, seed=10 + dim, two_node=True, singletons=True,
                            periodic=True))
    return out


def vec_funcs(shift=0.0):
    return vec(lambda *X: 1.0 + X[0] + shift, lambda *X: np.sin(X[-1]) - shift,
               lambda *X: X[0] * X[-1] + 2 * shift)


def cell_grid(geom):
    """Dense bool grid of the level's cells covered by its patches."""
    grid = np.zeros(tuple(geom.domain_cells), dtype=bool)
    for lo, hi in zip(geom.lower, geom.upper):  # test helper: per-patch loop is fine
        grid[tuple(slice(a, b + 1) for a, b in zip(lo, hi))] = True
    return grid


def coarsened(grid):
    out = grid
    for a in range(grid.ndim):
        out = np.logical_or.reduce(
            np.stack([np.take(out, np.arange(k, out.shape[a], 2), axis=a) for k in (0, 1)]),
            axis=0)
    return out


class TestSameLayout(unittest.TestCase):
    def test_binary_equals_numpy_on_buffers(self):
        for hier in layouts():
            a, b = random_snapshot(hier, XYZ, seed=1), random_snapshot(hier, XYZ, seed=2)
            s = random_snapshot(hier, seed=3)
            sv = comps(s, "value")
            for uf in (np.add, np.subtract, np.multiply, np.true_divide, np.arctan2):
                r = arith.binary(uf, a, b)
                self.assertIs(r.hier, a.hier)
                for c in XYZ:
                    assert_all_equal(comps(r, c), map(uf, comps(a, c), comps(b, c)))
            # scalar field broadcast over components, both orders, and Python scalars
            cases = [  # (result, expected from (component of a, value of s))
                (arith.binary(np.multiply, s, a), lambda x, v: v * x),
                (arith.binary(np.subtract, a, s), lambda x, v: x - v),
                (arith.binary(np.power, a, 2), lambda x, v: x**2),
                (arith.binary(np.subtract, 1.5, a), lambda x, v: 1.5 - x),
            ]
            self.assertEqual(cases[0][0].names, XYZ)
            for r, expect in cases:
                for c in XYZ:
                    assert_all_equal(comps(r, c), map(expect, comps(a, c), sv))

    def test_dtype_rules(self):
        hier = make_amr(2, 2, seed=5)
        f32 = random_snapshot(hier, dtype=np.float32)
        f64 = random_snapshot(hier, dtype=np.float64)
        i64 = snapshot(hier, lambda x, y: np.rint(10 * x), dtype=np.int64)
        cases = [
            ((f32, 2.0), np.float32), ((f32, np.float64(2.0)), np.float64),
            ((f32, np.array(2.0)), np.float64), ((f32, f64), np.float64),
            ((f32, f32), np.float32), ((i64, 2), np.int64), ((i64, 2.5), np.float64),
            ((2.0, f32), np.float32),
        ]
        for (x, y), dt in cases:
            r = arith.binary(np.add, x, y)
            self.assertEqual(r.levels[0].dtype, dt, (x if not isinstance(x, Snapshot) else
                                                      x.levels[0].dtype, y))
            self.assertEqual(r.levels[-1].buffer("value").dtype, dt)
        self.assertEqual(arith.binary(np.greater, f32, 0.0).levels[0].dtype, bool)

    def test_type_errors(self):
        hier = make_amr(2, 1, seed=1)
        s = random_snapshot(hier)
        v = random_snapshot(hier, XYZ)
        t = random_snapshot(hier, ("xx", "xy", "xz", "yy", "yz", "zz"))
        with self.assertRaises(TypeError) as ctx:
            arith.binary(np.add, s, np.ones(3))
        self.assertIn("finest()", str(ctx.exception))
        with self.assertRaises(TypeError) as ctx:
            arith.binary(np.add, v, t)
        self.assertIn("('x', 'y', 'z')", str(ctx.exception))
        self.assertIn("('xx'", str(ctx.exception))
        calls = [  # each raises TypeError
            *[(arith.binary, (np.add, s, bad), {}) for bad in ("a", [1.0], None)],
            (arith.binary, (np.add, 1.0, 2.0), {}),
            (arith.unary, (np.sqrt, s), dict(out=None)),
            (arith.dot, (s, v), {}), (arith.cross, (s, v), {}), (arith.modulus, (s,), {}),
        ]
        for func, args, kw in calls:
            with self.subTest(func=func.__name__, args=args), self.assertRaises(TypeError):
                func(*args, **kw)

    def test_different_grids_raise(self):
        a = random_snapshot(make_amr(2, 1, seed=1))
        for other in (make_amr(2, 1, seed=1, dl0=0.1), make_amr(2, 1, seed=1, periodic=True),
                      make_amr(2, 1, seed=1, domain_cells=(40, 40))):
            with self.assertRaises(ValueError) as ctx:
                arith.binary(np.add, a, random_snapshot(other))
            self.assertIn("different grids", str(ctx.exception))

    def test_extent_intersection(self):
        base = make_amr(2, 2, seed=6)
        a = random_snapshot(base.with_extent(PhysBox.parse(dict(x=(1.0, 6.0)), 2)))
        b = random_snapshot(base.with_extent(PhysBox.parse(dict(x=(2.0, 9.0), y=(0, 3)), 2)))
        r = arith.binary(np.add, a, b)
        self.assertEqual(r.hier.extent, PhysBox((2.0, 0.0), (6.0, 3.0)))
        self.assertIs(r.hier.levels[1], base.levels[1])
        self.assertIsNone(arith.binary(np.add, random_snapshot(base), 1.0).hier.extent)
        far = random_snapshot(base.with_extent(PhysBox.parse(dict(x=(7.0, 8.0)), 2)))
        with self.assertRaises(ValueError):
            arith.binary(np.add, a, far)


class TestDifferentLayouts(unittest.TestCase):
    def pairs(self):
        out = []
        for dim in (1, 2, 3):
            out.append((make_amr(dim, 3, seed=dim), make_amr(dim, 3, seed=50 + dim)))
            out.append((make_amr(dim, 3, seed=60 + dim, two_node=True, periodic=True),
                        make_amr(dim, 2, seed=70 + dim, singletons=True, periodic=True)))
        return out

    def test_result_layout_values_and_nesting(self):
        nonempty_fine = 0
        for ha, hb in self.pairs():
            a = snapshot(ha, vec_funcs(), time=1.0)
            b = snapshot(hb, vec_funcs(0.25), time=2.0)
            r = arith.binary(np.subtract, a, b)
            hr = r.hier
            self.assertEqual(r.time, 1.0)
            self.assertLessEqual(len(hr.levels), min(len(ha.levels), len(hb.levels)))
            self.assertEqual(hr.level_numbers(), tuple(range(len(hr.levels))))
            nonempty_fine += len(hr.levels) > 1
            # level 0 covers the domain, every level = intersection of both operands' cover
            self.assertEqual(hr.levels[0].size.sum() > 0, True)
            self.assertTrue(cell_grid(hr.levels[0]).all())
            for g in hr.levels:
                expect = cell_grid(ha.level(g.level)) & cell_grid(hb.level(g.level))
                np.testing.assert_array_equal(cell_grid(g), expect)
                cells = (g.upper - g.lower + 1).prod(axis=1).sum()
                self.assertEqual(cells, expect.sum())  # patches do not overlap
                arith.locate(ha.level(g.level), g)  # inside one patch of each operand
                arith.locate(hb.level(g.level), g)
            for g0, g1 in zip(hr.levels[:-1], hr.levels[1:]):  # properly nested
                self.assertFalse((coarsened(cell_grid(g1)) & ~cell_grid(g0)).any())
            # values: equal the operands at every result node (same coordinates)
            fa, fb = snapshot(hr, vec_funcs()), snapshot(hr, vec_funcs(0.25))
            for c in XYZ:
                assert_all_equal(comps(r, c), map(np.subtract, comps(fa, c), comps(fb, c)))
        self.assertGreater(nonempty_fine, 2)

    def test_integral_of_difference(self):
        for ha, hb in self.pairs():
            dim = ha.dim
            a1 = snapshot(ha, {"value": lambda *X: np.cos(X[0]) + X[-1] ** 2})
            a0 = snapshot(hb, {"value": lambda *X: np.sin(X[-1]) * X[0]})
            r = arith.binary(np.subtract, a1, a0)
            i1 = reduce.integral(arith.restrict(a1, r.hier))
            i0 = reduce.integral(arith.restrict(a0, r.hier))
            self.assertAlmostEqual(reduce.integral(r), i1 - i0,
                                   delta=1e-12 * (abs(i1) + abs(i0)), msg=f"{dim}d")

    def test_level_present_in_one_operand_dropped(self):
        ha = make_amr(2, 3, seed=2)
        hb = HierGeometry(make_amr(2, 3, seed=9).levels[:1])
        r = arith.binary(np.add, random_snapshot(ha), random_snapshot(hb))
        self.assertEqual(r.hier.level_numbers(), (0,))
        # disjoint level-1 patches: level 1 and everything above it are dropped
        g1 = ha.levels[1]
        free = np.argwhere(~cell_grid(g1)[::2, ::2] & ~cell_grid(g1)[1::2, 1::2])[0] * 2
        lone = type(g1).from_file_order(1, free[None, :], free[None, :] + 1, g1.dl,
                                        g1.domain_cells, g1.periodic)
        self.assertIsNone(arith.intersect_levels(g1, lone))
        hd = HierGeometry([ha.levels[0], lone])
        r = arith.binary(np.add, random_snapshot(ha), random_snapshot(hd))
        self.assertEqual(r.hier.level_numbers(), (0,))
        self.assertEqual(len(ha.levels), 3)

    def test_blocks_give_same_result(self):
        ha, hb = make_amr(3, 3, seed=3), make_amr(3, 3, seed=33)
        a = random_snapshot(ha, XYZ, np.float32, seed=1)
        b = random_snapshot(hb, ("value",), np.float32, seed=2)
        ref = arith.binary(np.multiply, a, b)
        with mock.patch.object(arith, "BLOCK_NODES", 29):
            small = arith.binary(np.multiply, a, b)
            restricted = arith.restrict(a, ref.hier)
        for c in XYZ:
            assert_all_equal(comps(small, c), comps(ref, c))
        assert_all_equal(comps(restricted, "y"), comps(arith.restrict(a, ref.hier), "y"))

    def test_restrict_and_align(self):
        ha, hb = make_amr(2, 3, seed=2), make_amr(2, 3, seed=12)
        a, b = random_snapshot(ha, seed=1), random_snapshot(hb, seed=2)
        hier, la, lb = arith.align(a, a)
        self.assertIs(hier, ha)
        self.assertTrue(all(x is y for x, y in zip(la, a.levels)))
        hier, la, lb = arith.align(a, b)
        added = [x.buffer("value") + y.buffer("value") for x, y in zip(la, lb)]
        assert_all_equal(added, comps(arith.binary(np.add, a, b), "value"))
        with self.assertRaises(ValueError):  # b's patches are not inside a's
            arith.restrict(a, hb)
        with self.assertRaises(ValueError):
            arith.restrict(a, make_amr(2, 1, seed=1, dl0=0.3))


class TestHelpers(unittest.TestCase):
    def setUp(self):
        self.hier = make_amr(3, 2, seed=7, two_node=True)
        self.a = random_snapshot(self.hier, XYZ, seed=1)
        self.b = random_snapshot(self.hier, XYZ, seed=2)
        self.s = random_snapshot(self.hier, seed=3)

    def stacked(self, snap, k):
        return np.stack([snap.levels[k].buffer(c) for c in snap.names])

    def test_vector_algebra(self):
        d, c, m = arith.dot(self.a, self.b), arith.cross(self.a, self.b), arith.modulus(self.a)
        for k in range(len(self.hier.levels)):
            A, B = self.stacked(self.a, k), self.stacked(self.b, k)
            np.testing.assert_allclose(d.levels[k].buffer("value"), (A * B).sum(axis=0),
                                       rtol=1e-14)
            np.testing.assert_allclose(self.stacked(c, k), np.cross(A.T, B.T).T, rtol=1e-12,
                                       atol=1e-14)
            np.testing.assert_allclose(m.levels[k].buffer("value"),
                                       np.sqrt((A * A).sum(axis=0)), rtol=1e-14)
        a32 = random_snapshot(self.hier, XYZ, np.float32)
        self.assertEqual(arith.modulus(a32).levels[0].dtype, np.float32)
        with mock.patch.object(arith, "BLOCK_NODES", 17):
            np.testing.assert_array_equal(arith.dot(self.a, self.b).levels[1].buffer("value"),
                                          d.levels[1].buffer("value"))

    def test_component_and_stack_share_buffers(self):
        x = arith.component(self.a, "y")
        self.assertEqual(x.names, ("value",))
        for got, ref in zip(comps(x, "value"), comps(self.a, "y")):
            self.assertTrue(np.shares_memory(got, ref))
        with self.assertRaises(KeyError):
            arith.component(self.a, "w")
        v = arith.stack_components([arith.component(self.a, c) for c in "zyx"], XYZ)
        for got, ref in zip(comps(v, "x"), comps(self.a, "z")):
            self.assertTrue(np.shares_memory(got, ref))
        assert_all_equal(comps(v, "z"), comps(self.a, "x"))
        with self.assertRaises(TypeError):
            arith.stack_components([self.a], ("x",))

    def test_where_and_clip(self):
        cond = arith.binary(np.greater, self.s, 0.0)
        w = arith.where(cond, self.a, 0.0)
        c = arith.clip(self.a, -0.5, 0.5)
        c_lo = arith.clip(self.a, self.s, None)
        sv = comps(self.s, "value")
        for comp in XYZ:
            av = comps(self.a, comp)
            assert_all_equal(comps(w, comp), [np.where(s > 0, a, 0) for a, s in zip(av, sv)])
            assert_all_equal(comps(c, comp), [np.clip(a, -0.5, 0.5) for a in av])
            assert_all_equal(comps(c_lo, comp), map(np.maximum, av, sv))
        with self.assertRaises(ValueError):
            arith.clip(self.a, None, None)


if __name__ == "__main__":
    unittest.main()
