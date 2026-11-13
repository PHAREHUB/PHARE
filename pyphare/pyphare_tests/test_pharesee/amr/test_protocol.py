#!/usr/bin/env python3
"""Tests of pharesee.amr._protocol: operators, ufuncs, numpy functions and what must raise.

Values are compared with numpy applied to the flat per-level buffers of the operands.
"""

import operator
import unittest

import numpy as np

from pyphare.pharesee.amr.field import (
    AMRField,
    ScalarField,
    VectorField,
    make_field,
)
from pyphare.pharesee.amr.ops import arith, deriv
from pyphare.pharesee.amr.store import Snapshot
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr


def funcs(t=0.0):
    return {
        "x": lambda x, y: 2.5 + 0.5 * x - 0.25 * y + t,
        "y": lambda x, y: 2.0 + 0.3 * x * y + t,
        "z": lambda x, y: 0.75 + np.sin(x) * np.cos(y) + t,
    }


def vector(seed=0, dtype=np.float64, t=0.0, hier=None):
    hier = hier or make_amr(2, 3, seed)
    return make_field([Snapshot(hier, tuple(fill(hier, funcs(t), dtype=dtype, time=t)))],
                      quantity="B")


def scalar(f, seed=0, dtype=np.float64, hier=None):
    hier = hier or make_amr(2, 3, seed)
    return make_field([Snapshot(hier, tuple(fill(hier, f, dtype=dtype)))], quantity="n")


def bufs(field, name=None):
    """Concatenation of every level's flat buffer of one component."""
    name = name or field.names[0]
    return np.concatenate([ld.buffer(name) for ld in field.levels_data()])


def fake_multi(times, seed=3):
    hier = make_amr(2, 2, seed)
    data = {float(t): fill(hier, funcs(t), time=t) for t in times}
    backend = FakeBackend(data, quantity="EM_B")
    return AMRField.from_backend(backend, times=list(times)), backend


class LegacyThing:
    """Stands for an object of master's pyphare.pharesee.hierarchy classes."""

    def __add__(self, other):
        return NotImplemented  # what the patched legacy operators do for AMRField operands

    __radd__ = __add__


LegacyThing.__module__ = "pyphare.pharesee.hierarchy.hierarchy"


class TestOperators(unittest.TestCase):
    def setUp(self):
        self.B = vector()
        self.a, self.b = self.B.x, self.B.y  # positive-valued on this layout
        self.assertGreater(bufs(self.a).min(), 0)
        self.assertGreater(bufs(self.b).min(), 0)

    def test_binary_arithmetic(self):
        ops = [operator.add, operator.sub, operator.mul, operator.truediv,
               operator.floordiv, operator.mod, operator.pow]
        A, Bv = bufs(self.a), bufs(self.b)
        for op in ops:
            with self.subTest(op=op.__name__):
                for res, want in ((op(self.a, self.b), op(A, Bv)),
                                  (op(self.a, 2.5), op(A, 2.5)),
                                  (op(2.5, self.a), op(2.5, A)),
                                  (op(self.a, np.float64(1.5)), op(A, 1.5)),
                                  (op(np.array(1.5), self.a), op(1.5, A))):
                    self.assertIsInstance(res, ScalarField)
                    np.testing.assert_allclose(bufs(res), want, rtol=1e-14)

    def test_result_classes(self):
        B = self.B
        self.assertIsInstance(B + B, VectorField)
        self.assertIsInstance(self.a * B, VectorField)
        self.assertIsInstance(B * self.a, VectorField)
        self.assertIsInstance(2 * B, VectorField)
        self.assertIsInstance(-B, VectorField)
        np.testing.assert_allclose(bufs(self.a * B, "z"), bufs(self.a) * bufs(B.z))
        lt = B < 2.0
        self.assertIsInstance(lt, VectorField)
        self.assertEqual(lt.dtype, np.bool_)
        dot = make_field([arith.dot(B.snapshot(), B.snapshot())])
        self.assertIsInstance(dot, ScalarField)
        self.assertIsInstance(make_field([arith.cross(B.snapshot(), B.snapshot())]),
                              VectorField)
        np.testing.assert_allclose(bufs(dot), bufs(B.x)**2 + bufs(B.y)**2 + bufs(B.z)**2)
        with self.assertRaises(TypeError):
            B + make_field([arith.stack_components([self.a.snapshot()] * 6,
                                                   ("xx", "xy", "xz", "yy", "yz", "zz"))])

    def test_comparisons(self):
        A, Bv = bufs(self.a), bufs(self.b)
        for op in (operator.lt, operator.le, operator.gt, operator.ge):
            with self.subTest(op=op.__name__):
                r = op(self.a, self.b)
                self.assertIsInstance(r, ScalarField)
                self.assertEqual(r.dtype, np.bool_)
                np.testing.assert_array_equal(bufs(r), op(A, Bv))
                np.testing.assert_array_equal(bufs(op(2.0, self.a)), op(2.0, A))

    def test_unary_and_bitwise(self):
        A = bufs(self.a)
        np.testing.assert_array_equal(bufs(-self.a), -A)
        np.testing.assert_array_equal(bufs(+self.a), A)
        np.testing.assert_array_equal(bufs(abs(self.a - 3.0)), np.abs(A - 3.0))
        m1, m2 = self.a > 3.0, self.b > 3.0
        M1, M2 = A > 3.0, bufs(self.b) > 3.0
        np.testing.assert_array_equal(bufs(~m1), ~M1)
        np.testing.assert_array_equal(bufs(m1 & m2), M1 & M2)
        np.testing.assert_array_equal(bufs(m1 | m2), M1 | M2)
        np.testing.assert_array_equal(bufs(m1 ^ m2), M1 ^ M2)
        n = scalar(lambda x, y: np.floor(3 * x + y), dtype=np.int64)
        N = bufs(n)
        np.testing.assert_array_equal(bufs(n << 2), N << 2)
        np.testing.assert_array_equal(bufs(n >> 1), N >> 1)
        np.testing.assert_array_equal(bufs(1 << n % 5), 1 << N % 5)

    def test_ufuncs(self):
        A, Bv = bufs(self.a), bufs(self.b)
        for uf in (np.sqrt, np.exp, np.sin, np.log, np.negative, np.isnan):
            np.testing.assert_allclose(bufs(uf(self.a)), uf(A))
        for uf in (np.add, np.maximum, np.arctan2, np.hypot, np.equal):
            np.testing.assert_allclose(bufs(uf(self.a, self.b)), uf(A, Bv))
        eq = np.equal(self.a, self.a)
        self.assertIsInstance(eq, ScalarField)
        self.assertTrue(bufs(eq).all())
        v = np.sqrt(abs(self.B))
        self.assertIsInstance(v, VectorField)
        np.testing.assert_allclose(bufs(v, "y"), np.sqrt(Bv))

    def test_dtype_rules(self):
        B32 = vector(dtype=np.float32)
        self.assertEqual((B32 * 2.0).dtype, np.float32)
        self.assertEqual((2.0 * B32.x).dtype, np.float32)
        self.assertEqual((B32 * np.float32(2.0)).dtype, np.float32)
        self.assertEqual((B32 * np.float64(2.0)).dtype, np.float64)
        self.assertEqual((B32 + self.B).dtype, np.float64)
        self.assertEqual(np.sqrt(B32.x).dtype, np.float32)


class TestRaises(unittest.TestCase):
    def setUp(self):
        self.B = vector(seed=1)
        self.a = self.B.x

    def test_ufunc_misuse(self):
        a = self.a
        with self.assertRaisesRegex(TypeError, "returned value"):
            np.add(a, 1.0, out=np.empty(3))
        with self.assertRaisesRegex(TypeError, "np.where"):
            np.add(a, 1.0, where=True)
        with self.assertRaisesRegex(TypeError, "integral"):
            np.add.reduce(a)
        with self.assertRaisesRegex(TypeError, "finest"):
            np.add.accumulate(a)
        with self.assertRaisesRegex(TypeError, "finest"):
            np.multiply.outer(a, a)
        with self.assertRaisesRegex(TypeError, "outputs"):
            np.modf(a)
        with self.assertRaisesRegex(TypeError, "outputs"):
            np.divmod(a, 2.0)
        with self.assertRaisesRegex(TypeError, "//"):
            divmod(a, 2.0)
        with self.assertRaisesRegex(TypeError, "dot"):
            a @ a
        with self.assertRaisesRegex(TypeError, "keywords"):
            np.add(a, 1.0, dtype=np.float32)

    def test_conversions(self):
        a = self.a
        with self.assertRaisesRegex(TypeError, "np.any"):
            bool(a)
        with self.assertRaisesRegex(TypeError, "np.any"):
            if a > 1.0:
                pass
        with self.assertRaisesRegex(TypeError, "finest"):
            iter(a)
        with self.assertRaisesRegex(TypeError, "finest"):
            list(self.B)
        with self.assertRaisesRegex(TypeError, "nbytes"):
            len(a)
        with self.assertRaisesRegex(TypeError, "finest"):
            np.asarray(a)
        with self.assertRaisesRegex(TypeError, "finest"):
            np.array(self.B)

    def test_array_operands(self):
        for other in (np.ones(3), np.ones((2, 2)), [1.0, 2.0], (1.0,)):
            with self.subTest(other=type(other).__name__):
                with self.assertRaisesRegex(TypeError, "finest"):
                    self.a + other
                with self.assertRaisesRegex(TypeError, "finest"):
                    other * self.a
        with self.assertRaises(TypeError):
            self.a + "text"
        np.testing.assert_allclose(bufs(self.a + np.array(2.0)), bufs(self.a) + 2.0)

    def test_legacy_operands(self):
        old = LegacyThing()
        for call in (lambda: self.a + old, lambda: old + self.a, lambda: self.a * old,
                     lambda: np.add(old, self.a), lambda: np.where(self.a > 1, old, 0.0),
                     lambda: np.isclose(self.a, old)):
            with self.assertRaisesRegex(TypeError, "legacy reader") as cm:
                call()
            self.assertIn("same path", str(cm.exception))

    def test_unsupported_numpy_functions(self):
        a = self.a
        with self.assertRaisesRegex(TypeError, "integral"):
            np.sum(a)
        for f in (np.std, np.var):
            with self.assertRaisesRegex(TypeError, "mean"):
                f(a)
        for f in (np.median, np.cumsum, np.reshape, np.concatenate, np.sort, np.shape):
            with self.subTest(f=f.__name__):
                args = ([a, a],) if f is np.concatenate else (a, -1) if f is np.reshape else (a,)
                with self.assertRaisesRegex(TypeError, "finest"):
                    f(*args)
        for f in (np.mean, np.min, np.max, np.any, np.all):
            with self.assertRaisesRegex(TypeError, "axis"):
                f(a, axis=0)
        with self.assertRaisesRegex(TypeError, "np.where"):
            np.where(a > 1)
        with self.assertRaisesRegex(TypeError, "scalar"):
            np.gradient(self.B)
        with self.assertRaisesRegex(TypeError, "spacing"):
            np.gradient(a, 0.1)
        with self.assertRaisesRegex(TypeError, "keywords"):
            np.mean(a, keepdims=True)


class TestIdentityEquality(unittest.TestCase):
    def test_eq_is_identity_without_reading(self):
        hier = make_amr(2, 2, 4)
        backend = FakeBackend(fill(hier, funcs()), quantity="EM_B")
        B = AMRField.from_backend(backend)
        B2 = AMRField.from_backend(backend)
        self.assertFalse(B == B2)
        self.assertTrue(B != B2)
        self.assertTrue(B == B)
        self.assertFalse(B.x == B.x)  # two component objects
        self.assertEqual(len({B, B2, B}), 2)
        self.assertEqual({B: 1}[B], 1)
        self.assertEqual(hash(B), hash(B))
        self.assertIn(B, [B2, B])
        self.assertIn(B, [None, "B", B])  # non-numeric objects compare by identity
        self.assertFalse(B == None)  # noqa: E711 - the comparison itself is tested
        self.assertTrue(B != "B")
        self.assertEqual(backend.calls, 0)
        self.assertIs(type(B).__hash__, object.__hash__)

    def test_eq_with_number_or_array_raises(self):
        dens = scalar(lambda x, y: 1.0 + x * y)
        for other in (1.0, 2, np.float64(1.0), np.float32(1.0), np.array(1.0), True,
                      np.ones(3), [1.0, 2.0], (1.0,)):
            with self.subTest(other=repr(other)):
                with self.assertRaisesRegex(TypeError, r"np\.equal.*np\.isclose"):
                    dens == other  # noqa: B015 - the comparison must raise
                with self.assertRaisesRegex(TypeError, "np.not_equal"):
                    dens != other  # noqa: B015
        # the reviewer's trap: locating the maximum with == must not silently give 0
        with self.assertRaisesRegex(TypeError, "np.equal"):
            np.where(dens == dens.max(), 1, 0)

    def test_elementwise_equality_functions_return_bool_fields(self):
        N = scalar(lambda x, y: np.where(x > 2.0, 1.0, 0.0))
        v = bufs(N)
        for func in (np.equal, np.not_equal):
            got = func(N, 1.0)
            self.assertIsInstance(got, ScalarField)
            self.assertEqual(got.dtype, np.bool_)
            np.testing.assert_array_equal(bufs(got), func(v, 1.0))
        close = np.isclose(N, 1.0 + 1e-12)
        self.assertEqual(close.dtype, np.bool_)
        np.testing.assert_array_equal(bufs(close), np.isclose(v, 1.0 + 1e-12))
        # where the maximum is: the elementwise idiom
        at_max = np.where(np.equal(N, N.max()), 1, 0)
        np.testing.assert_array_equal(bufs(at_max), np.where(v == v.max(), 1, 0))

    def test_numpy_scalar_on_the_left_is_elementwise(self):
        """Documented: numpy turns np.float64(x) == f into np.equal(x, f) before asking the
        field, so it gives a field of bool (the same as np.equal), not a TypeError."""
        N = scalar(lambda x, y: np.where(x > 2.0, 1.0, 0.0))
        got = np.float64(1.0) == N
        self.assertIsInstance(got, ScalarField)
        np.testing.assert_array_equal(bufs(got), bufs(np.equal(1.0, N)))


class TestNumpyFunctions(unittest.TestCase):
    def setUp(self):
        self.B = vector(seed=2)
        self.a, self.b = self.B.x, self.B.z

    def test_reductions(self):
        try:
            from pyphare.pharesee.amr.ops import reduce
        except ImportError:
            self.skipTest("ops/reduce.py not written yet")
        a, snap = self.a, self.a.snapshot()
        self.assertEqual(np.mean(a), reduce.mean(snap))
        self.assertEqual(np.mean(a), a.mean())
        for f in (np.min, np.amin, np.nanmin):
            self.assertEqual(f(a), reduce.minimum(snap))
        for f in (np.max, np.amax, np.nanmax):
            self.assertEqual(f(a), reduce.maximum(snap))
        self.assertEqual(np.nanmean(a), a.mean())
        np.testing.assert_array_equal(np.mean(self.B), reduce.mean(self.B.snapshot()))
        holed = np.where(a > 2.0, np.nan, a)
        with self.assertRaisesRegex(ValueError, "NaN"):
            np.nanmean(holed)

    def test_gradient_copy_clip(self):
        g = np.gradient(self.a)
        self.assertIsInstance(g, VectorField)
        ref = deriv.grad(self.a.snapshot())
        for k, ld in enumerate(g.levels_data()):
            for c in "xyz":
                np.testing.assert_array_equal(ld.buffer(c), ref.levels[k].buffer(c))
        c = np.copy(self.a)
        self.assertIsInstance(c, ScalarField)
        c.levels_data()[0].buffer("value")[0] = 1e9
        self.assertNotEqual(bufs(self.a)[0], 1e9)
        A = bufs(self.a)
        np.testing.assert_array_equal(bufs(np.clip(self.a, 1.0, 2.0)), np.clip(A, 1.0, 2.0))
        np.testing.assert_array_equal(bufs(np.clip(self.a, None, 2.0)), np.clip(A, None, 2.0))
        np.testing.assert_array_equal(bufs(np.clip(self.a, min=1.2)), np.clip(A, 1.2, None))
        np.testing.assert_array_equal(bufs(np.clip(self.a, self.b, 3.0)),
                                      np.clip(A, bufs(self.b), 3.0))

    def test_where_isclose_any_all(self):
        A, Bz = bufs(self.a), bufs(self.b)
        w = np.where(self.a > 2.0, self.a, 0.0)
        self.assertIsInstance(w, ScalarField)
        np.testing.assert_array_equal(bufs(w), np.where(A > 2.0, A, 0.0))
        wv = np.where(self.a > 2.0, self.B, -1.0)
        self.assertIsInstance(wv, VectorField)
        np.testing.assert_array_equal(bufs(wv, "z"), np.where(A > 2.0, Bz, -1.0))
        near = np.isclose(self.a, self.a + 1e-9)
        self.assertIsInstance(near, ScalarField)
        self.assertTrue(bufs(near).all())
        np.testing.assert_array_equal(bufs(np.isclose(self.a, self.b, atol=0.1)),
                                      np.isclose(A, Bz, atol=0.1))
        self.assertIs(np.allclose(self.a, self.a + 1e-9), True)
        self.assertIs(np.allclose(self.a, self.a + 1.0), False)
        self.assertIs(np.allclose(self.B, self.B * (1 + 1e-12)), True)
        self.assertIs(np.any(self.a > A.max() - 1e-9), True)
        self.assertIs(np.any(self.a > A.max()), False)
        self.assertIs(np.all(self.a > A.min() - 1.0), True)
        self.assertIs(np.all(self.a > A.min()), False)


class TestMultiTime(unittest.TestCase):
    def test_pairs_equal_time_lists(self):
        B, backend = fake_multi([0.5, 0.0])
        self.assertEqual(backend.calls, 0)
        s = B + B
        np.testing.assert_array_equal(s.times(), [0.5, 0.0])
        for t in (0.0, 0.5):
            np.testing.assert_allclose(bufs(s.at(t), "y"), 2 * bufs(B.at(t), "y"))
        u = -B.x
        np.testing.assert_array_equal(u.times(), [0.5, 0.0])
        np.testing.assert_allclose(bufs(u.at(0.5)), -bufs(B.x.at(0.5)))
        np.testing.assert_allclose(bufs((B * 2.0).at(0.0), "x"), 2 * bufs(B.at(0.0), "x"))

    def test_different_lists_raise(self):
        B, _ = fake_multi([0.0, 0.5])
        C, _ = fake_multi([0.5, 0.0])
        D, _ = fake_multi([0.0, 0.5, 1.0])
        for other in (C, D, D.at(1.0)):
            with self.assertRaisesRegex(ValueError, "at"):
                B + other
        with self.assertRaisesRegex(ValueError, "0.5"):
            np.mean(B.x)
        with self.assertRaisesRegex(ValueError, "0.5"):
            B.x.max()

    def test_single_times_combine(self):
        B, _ = fake_multi([0.0, 1.0])
        d = B.at(1.0) - B.at(0.0)
        np.testing.assert_array_equal(d.times(), [1.0])
        np.testing.assert_allclose(bufs(d, "x"), 1.0, rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
