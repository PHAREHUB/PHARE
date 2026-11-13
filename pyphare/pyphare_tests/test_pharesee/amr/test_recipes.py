"""Unit tests of ``pharesee.amr.recipes`` (derived quantities served by Run), on synthetic data.

Analytic fields on random valid AMR layouts (``testing.make_amr``) behind an in-memory
backend that counts reads (``testing.FakeBackend``): no file, no build needed.
"""

import os
import unittest
from types import SimpleNamespace

import numpy as np

from pyphare.pharesee.amr import recipes
from pyphare.pharesee.amr.field import AMRField, ScalarField, VectorField
from pyphare.pharesee.amr.geometry import GeometryCache
from pyphare.pharesee.amr.ops.deriv import one_sided_mask
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr

os.environ.setdefault("PHARESEE_CACHE", "0")


def _backend(dim, funcs, names, nlevels=3, seed=4, periodic=True, **kw):
    hier = make_amr(dim, nlevels, seed, periodic=periodic, **kw)
    return FakeBackend(fill(hier, funcs, names=names))


def _vector_backend(dim, fx, fy, fz, **kw):
    return _backend(dim, {"x": fx, "y": fy, "z": fz}, None, **kw)


class Opener:
    """The ``opener`` Run gives to recipes, over fake backends; records every call."""

    def __init__(self, backends):
        self.backends = backends
        self.calls = []
        self.cache = GeometryCache()

    def __call__(self, stem, times, extent, **kw):
        self.calls.append((stem, times, extent, dict(kw)))
        kw.pop("derived", None)  # fake backends compute nothing
        kw.pop("names", None)
        return AMRField.from_backend(self.backends[stem], times=times, extent=extent,
                                     cache=self.cache, **kw)


def _nodes(field):
    """{(level, lower, upper): {component: patch array}} of a single-time field."""
    out = {}
    snap = field.snapshot()
    for ld in snap.levels:
        g = ld.geom
        for slot in range(g.npatch):
            key = (g.level, tuple(g.lower[slot]), tuple(g.upper[slot]))
            out[key] = {n: np.array(g.view(ld.buffer(n), slot)) for n in ld.names}
    return out


class CurrentRecipeTest(unittest.TestCase):
    def test_curl_of_linear_field_is_exact_and_labelled(self):
        # B = (0, a x, -b y): J = curl B = (-b, 0, a) at every node, including patch faces
        a, b = 1.5, 0.25
        be = _vector_backend(2, lambda x, y: 0 * x, lambda x, y: a * x, lambda x, y: -b * y,
                             periodic=False)
        op = Opener({"EM_B": be})
        J = recipes.evaluate(recipes.RECIPES["J"], op, 0.0, None, be.meta().dl0,
                             {"EM_B": "vtkhdf"})
        self.assertIsInstance(J, VectorField)
        self.assertEqual((J.quantity, J.stencil), ("J", "node"))
        for comps in _nodes(J).values():
            np.testing.assert_allclose(comps["x"], -b, rtol=1e-12, atol=1e-12)
            np.testing.assert_allclose(comps["y"], 0.0, atol=1e-12)
            np.testing.assert_allclose(comps["z"], a, rtol=1e-12)


class DivergenceRecipeTest(unittest.TestCase):
    def test_divergence_of_linear_field(self):
        # div (2x, 3y, z) in 2D = 5 everywhere (z does not vary along a missing axis)
        be = _vector_backend(2, lambda x, y: 2 * x, lambda x, y: 3 * y, lambda x, y: x,
                             periodic=False)
        d = recipes.evaluate(recipes.RECIPES["divB"], Opener({"EM_B": be}), 0.0, None,
                             be.meta().dl0, {"EM_B": "vtkhdf"}, mask_faces=False)
        self.assertIsInstance(d, ScalarField)
        self.assertEqual((d.quantity, d.stencil), ("divB", "node"))
        for comps in _nodes(d).values():
            np.testing.assert_allclose(comps["value"], 5.0, rtol=1e-12)

    def test_mask_faces_default_depends_on_the_format(self):
        """vtkhdf: one-sided nodes are NaN by default (the file cannot give divB there);
        h5: values kept (Yee divB from the file's ghosts)."""
        self.assertTrue(recipes.default_mask_faces("vtkhdf"))
        self.assertFalse(recipes.default_mask_faces("h5"))
        be = _vector_backend(2, lambda x, y: 2 * x, lambda x, y: 3 * y, lambda x, y: x,
                             periodic=False, nlevels=2)
        d = recipes.evaluate(recipes.RECIPES["divB"], Opener({"EM_B": be}), 0.0, None,
                             be.meta().dl0, {"EM_B": "vtkhdf"})
        masked = 0
        for ld in d.snapshot().levels:
            mask = one_sided_mask(ld.geom)
            v = ld.buffer("value")
            np.testing.assert_array_equal(np.isnan(v), mask)
            np.testing.assert_allclose(v[~mask], 5.0, rtol=1e-12)
            masked += int(mask.sum())
        self.assertGreater(masked, 0)
        # what the reductions do on the masked field (documented in Run.GetDivB)
        self.assertTrue(np.isnan(d.mean()))
        self.assertAlmostEqual(float(d.max()), 5.0, places=10)  # min/max ignore the masked nodes
        self.assertAlmostEqual(float(d.min()), 5.0, places=10)
        self.assertAlmostEqual(float(d.global_max("value")), 5.0, places=10)
        op = Opener({"EM_B": _backend(2, lambda x, y: x + y, None)})
        dh = recipes.evaluate(recipes.RECIPES["divB"], op, 0.0, None, be.meta().dl0,
                              {"EM_B": "h5"})
        self.assertEqual(op.calls[0][3]["derived"], "divB_yee")  # h5: no mask, Yee values
        self.assertEqual((dh.quantity, dh.stencil), ("divB", "yee"))


class PressureRecipeTest(unittest.TestCase):
    def _sim(self, Te):
        return lambda: SimpleNamespace(electrons=SimpleNamespace(
            closure=SimpleNamespace(Te=Te)))

    def test_pe_is_ni_times_te(self):
        be = _backend(2, lambda x, y: 1 + x * y, None)
        Pe = recipes.evaluate(recipes.RECIPES["Pe"], Opener({"ions_charge_density": be}),
                              0.0, None, be.meta().dl0, {"ions_charge_density": "h5"},
                              self._sim(0.12))
        self.assertEqual(Pe.quantity, "Pe")
        src = _nodes(AMRField.from_backend(be))
        for key, comps in _nodes(Pe).items():
            np.testing.assert_array_equal(comps["value"], src[key]["value"] * 0.12)


class MagneticFluxTest(unittest.TestCase):
    def test_flux_of_uniform_by_and_linear_bx(self):
        # Az = -int_0^x By(x', 0) dx' + int_0^y Bx dy' = -c x + y^2 / 2  (trapezoids exact)
        c = 0.7
        be = _vector_backend(2, lambda x, y: y, lambda x, y: c + 0 * x, lambda x, y: 0 * x,
                             nlevels=3, seed=2)
        B = AMRField.from_backend(be)
        Az, (xn, yn) = recipes.magnetic_flux(B)
        X, Y = np.meshgrid(xn, yn, indexing="ij")
        np.testing.assert_allclose(Az, -c * X + Y**2 / 2, atol=1e-12)
        np.testing.assert_allclose(np.diff(xn), be.meta().dl0[0])  # level-0 grid
        # only level 0 was read (finest(max_level=0))
        self.assertEqual({lvl for (_, lvl, _) in be.requests}, {0})

    def test_flux_on_given_grid_interpolates_linearly(self):
        be = _vector_backend(2, lambda x, y: y, lambda x, y: 0.5 + 0 * x, lambda x, y: 0 * x)
        B = AMRField.from_backend(be)
        xn = np.linspace(0, 9.0, 31)
        yn = np.linspace(0, 7.5, 17)
        Az, (gx, _) = recipes.magnetic_flux(B, interp="bilinear", xn=xn, yn=yn)
        X, Y = np.meshgrid(xn, yn, indexing="ij")
        np.testing.assert_array_equal(gx, xn)
        np.testing.assert_allclose(Az, -0.5 * X + Y**2 / 2, atol=1e-12)

    def test_flux_needs_2d(self):
        be = _vector_backend(1, lambda x: x, lambda x: x, lambda x: x)
        with self.assertRaisesRegex(ValueError, "2D"):
            recipes.magnetic_flux(AMRField.from_backend(be))


class RegistryTest(unittest.TestCase):
    def test_registry_names_inputs_and_halos(self):
        R = recipes.RECIPES
        self.assertEqual(set(R), {"J", "J_yee", "divB", "Pe", "Az"})
        self.assertEqual((R["J"].inputs, R["J"].halo), (("EM_B",), 1))
        self.assertEqual((R["divB"].inputs, R["divB"].halo), (("EM_B",), 1))
        self.assertEqual(R["J_yee"].halo, 0)
        self.assertEqual(R["Pe"].inputs, ("ions_charge_density",))

    def test_multi_time_curl_keeps_caller_order(self):
        hier = make_amr(2, 2, 3, periodic=False)  # t*y is not periodic
        data = {t: fill(hier, {"x": lambda x, y, t=t: t * y, "y": lambda x, y: 0 * x,
                               "z": lambda x, y: 0 * x}, time=t) for t in (0.0, 1.0, 2.0)}
        be = FakeBackend(data)
        J = recipes.evaluate(recipes.RECIPES["J"], Opener({"EM_B": be}), [2.0, 0.0], None,
                             be.meta().dl0, {"EM_B": "vtkhdf"})
        np.testing.assert_array_equal(J.times(), [2.0, 0.0])
        # Jz = -dBx/dy = -t
        for t in (2.0, 0.0):
            for ld in J.snapshot(t).levels:
                np.testing.assert_allclose(ld.buffer("z"), -t, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
