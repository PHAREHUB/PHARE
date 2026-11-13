"""Acceptance checks of the 10-19 pharesee API review issues #1446-#1458 (DESIGN 14).

One test per issue. Each docstring states what master does wrong, the test asserts what the
new reader does. Synthetic fields (``testing.make_amr`` + ``FakeBackend``, which counts
reads) unless a real file is needed (then skipped without ``build/sample_data``).
"""

import os
import unittest

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from pyphare.core.operators import cross, dot, grad, modulus, sqrt
from pyphare.pharesee.amr.boxes import PhysBox
from pyphare.pharesee.amr.field import AMRField, ScalarField, VectorField
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr
from pyphare.pharesee.amr.views import plot as plotmod
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES

os.environ.setdefault("PHARESEE_CACHE", "0")
plotmod.install()


def _field(f, dim=2, nlevels=3, seed=3, names=None, periodic=False, **kw):
    hier = make_amr(dim, nlevels, seed, periodic=periodic, **kw)
    be = FakeBackend(fill(hier, f, names=names))
    return AMRField.from_backend(be), be


def _vector(fx, fy, fz, **kw):
    return _field({"x": fx, "y": fy, "z": fz}, **kw)


def _all_values(field, name="value"):
    return np.concatenate([ld.buffer(name) for ld in field.snapshot().levels])


class AcceptanceTest(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_1446_grad_divides_by_the_mesh_size(self):
        """#1446: master's grad (core/operators.py:54-79) calls np.gradient(dataset) without
        spacing, so grad(a x + b y) = (a dx, b dy) and only in 2D. New: (a, b, 0) at every
        node of every level, in 1D, 2D and 3D, also through np.gradient."""
        a, b = 3.0, -2.0
        n, _ = _field(lambda x, y: a * x + b * y, dl0=(0.2, 0.5))
        for g in (grad(n), np.gradient(n)):
            self.assertIsInstance(g, VectorField)
            np.testing.assert_allclose(_all_values(g, "x"), a, rtol=1e-12)
            np.testing.assert_allclose(_all_values(g, "y"), b, rtol=1e-12)
            np.testing.assert_array_equal(_all_values(g, "z"), 0.0)
        n3, _ = _field(lambda x, y, z: a * z, dim=3, nlevels=2)
        np.testing.assert_allclose(_all_values(grad(n3), "z"), a, rtol=1e-12)

    def test_1447_operands_on_different_layouts_are_aligned_by_position_in_space(self):
        """#1447: master pairs patches by their index in the patch list (compute_hier_from,
        hierarchy_utils.py:173-209), so A - B mixes unrelated patches when the layouts differ
        (another time, another refinement). New: the result lives on the intersections of
        A's and B's patches; the same function on two layouts subtracts to 0 everywhere."""
        def f(x, y):
            return np.sin(x) * np.cos(y)

        A, _ = _field(f, seed=1)
        B, _ = _field(f, seed=2)
        self.assertNotEqual(A.snapshot().levels[1].geom.key(), B.snapshot().levels[1].geom.key())
        d = A - B
        for ld in d.snapshot().levels:
            np.testing.assert_allclose(ld.buffer("value"), 0.0, atol=1e-15)
        self.assertAlmostEqual(d.integral(), 0.0, places=12)
        # the fine levels of the result are inside both operands' fine levels
        for ga, gd in zip(A.snapshot().hier.levels[1:], d.snapshot().hier.levels[1:]):
            self.assertLessEqual(int(gd.size.sum()), int(ga.size.sum()) * 4)

    def test_1448_out_argument_raises_instead_of_crashing(self):
        """#1448: np.sqrt(h, out=...) on a master hierarchy reached C code with a wrong
        buffer (segfault). New: out= (and where=) raise TypeError naming what to do."""
        n, _ = _field(lambda x, y: 1.0 + x * x)
        with self.assertRaisesRegex(TypeError, "out="):
            np.sqrt(n, out=n)
        with self.assertRaisesRegex(TypeError, "out="):
            np.add(n, n, out=np.empty(3))
        np.testing.assert_allclose(_all_values(np.sqrt(n)), np.sqrt(_all_values(n)))
        np.testing.assert_allclose(_all_values(sqrt(n)), np.sqrt(_all_values(n)))

    def test_1449_operations_on_file_data(self):
        """#1449: master's moments are h5py Datasets (fromh5.py:110,124): Vi / Ni crashed
        on file data. New: file-backed fields combine like in-memory ones (vector / scalar
        broadcasts over components) and per-patch views hold numpy arrays."""
        V, _ = _vector(lambda x, y: x, lambda x, y: 2 * y, lambda x, y: x * y, seed=4)
        N, _ = _field(lambda x, y: 1.0 + x * y, seed=4)
        self.assertTrue(all(ld.backend is not None for ld in V.snapshot().levels))
        q = V / N
        self.assertIsInstance(q, VectorField)
        for c in "xyz":
            np.testing.assert_allclose(_all_values(q, c), _all_values(V, c) / _all_values(N))
        if os.path.isdir(os.path.join(SAMPLES, "2d_tagged", "vtk")):
            from pyphare.pharesee.run import Run

            r = Run(os.path.join(SAMPLES, "2d_tagged", "h5"), legacy=False)
            ratio = r.GetVi(0.1) / r.GetNi(0.1)
            self.assertIsInstance(ratio, VectorField)
            ds = r.GetNi(0.1).levels()[0].patches[0].patch_datas["value"].dataset
            self.assertIsInstance(ds, np.ndarray)

    def test_1450_colour_range_is_that_of_the_shown_data(self):
        """#1450: master's plot2d takes vmin/vmax from global_min/max over every level and
        ghost node of the whole hierarchy (hierarchy.py:476-477, 506-507): a zoom on a quiet
        region is one flat colour. New: min/max of the data on screen."""
        n, _ = _field(lambda x, y: np.where(x < 2.0, 100.0 * y, np.sin(y)), seed=6)
        xlim, ylim = (5.0, 8.0), (2.0, 6.0)
        _, ax = n.plot(xlim=xlim, ylim=ylim)
        g = n.finest(box=dict(x=xlim, y=ylim))
        self.assertEqual(ax.get_images()[0].get_clim(),
                         (float(np.nanmin(g.data)), float(np.nanmax(g.data))))
        self.assertLess(ax.get_images()[0].get_clim()[1], 1.01)  # not 100 * y

    def test_1451_float_boxes_are_not_truncated(self):
        """#1451: selection_box is turned into a pyphare Box, which casts its bounds to int
        (core/box.py:15-16): Box((10, 5), (18, 12.5)).upper == [18, 12], so a patch starting
        between y=12 and 12.5 is lost. New: box=dict(...) keeps float bounds; the selected
        patches are exactly those whose physical extent overlaps the box."""
        from pyphare.core.box import Box

        self.assertEqual(list(Box((10, 5), (18, 12.5)).upper), [18, 12])  # the old trap
        hier = make_amr(2, 2, 7, domain_cells=(80, 80), dl0=0.2)
        g0 = hier.levels[0]
        starts = np.unique(g0.lower[:, 1] * g0.dl[1])
        y_edge = float(next(s for s in starts if s % 1.0 > 0.05))  # e.g. 3.2
        box = PhysBox((1.1, 0.3), (8.7, y_edge + 0.1))
        be = FakeBackend(fill(hier, lambda x, y: x))
        f = AMRField.from_backend(be, extent=box.to_dict())
        self.assertEqual(f.extent, box)
        for g_full, g_sel in zip(hier.levels, f.snapshot().hier.levels):
            lo, hi = g_full.lower * g_full.dl, (g_full.upper + 1) * g_full.dl
            want = ((hi > np.array(box.lo)) & (lo < np.array(box.hi))).all(axis=1)
            got = {tuple(v) for v in g_sel.lower.tolist()}
            self.assertEqual(got, {tuple(v) for v in g_full.lower[want].tolist()})
        rows = f.snapshot().hier.levels[0].lower[:, 1] * g0.dl[1]
        self.assertIn(y_edge, rows.tolist())  # the patch Box truncation would lose

    def test_1452_each_patch_is_read_once(self):
        """#1452: master's GetB(all_primal=True) reads the Yee datasets, then GetJ re-reads B
        through GetB(all_primal=False), and every merged/compute pass re-reads h5py datasets
        (2-4 reads per patch). New: one read per patch, kept for every later use."""
        B, be = _vector(lambda x, y: x, lambda x, y: y, lambda x, y: x * y, seed=9)
        from pyphare.pharesee.amr.ops import deriv

        J = B._like([deriv.curl(B.snapshot())])
        _ = dot(B, B), cross(B, J), modulus(B), B.x + B.y, B.mean(), B.finest()
        for (_, level), counts in be.reads.items():
            np.testing.assert_array_equal(counts, 1, f"level {level}")
        self.assertEqual(sum(c.size for c in be.reads.values()),
                         sum(g.npatch for g in B.snapshot().hier.levels))

    def test_1453_small_plot_bugs(self):
        """#1453: master's plot silently does nothing on a 3D field (hierarchy.py:543-547),
        plot1d returns None and slices [ng:-ng], empty when ng == 0. New: 3D raises asking
        for a cut; 1D returns (fig, ax) with every node value drawn."""
        f3, _ = _field(lambda x, y, z: x, dim=3, nlevels=2)
        with self.assertRaisesRegex(ValueError, "cut"):
            f3.plot()
        f1, _ = _field(lambda x: x * x, dim=1, nlevels=1)
        fig, ax = f1.plot()
        self.assertIsNotNone(fig)
        y = ax.lines[0].get_ydata()
        self.assertEqual(int(np.isfinite(y).sum()), f1.snapshot().hier.levels[0].nnodes)

    def test_1454_functions_and_reductions_in_3d(self):
        """#1454: hierarchy_from(func=...) (fromfunc.py:33-104) and master's merged views
        stop at 2D. New: 3D fields are first class: an analytic 3D field integrates and
        samples exactly (multilinear), on every level."""
        def f(x, y, z):
            return 1.0 + 2.0 * x - y + 0.5 * z

        n, _ = _field(f, dim=3, nlevels=2, seed=2)
        dom = n.domain
        exact = np.prod(dom) * (1.0 + dom[0] - dom[1] / 2 + dom[2] / 4)
        self.assertAlmostEqual(n.integral() / exact, 1.0, places=12)
        rng = np.random.default_rng(0)
        pts = rng.uniform(0, 1, (200, 3)) * dom
        np.testing.assert_allclose(n(points=pts), f(*pts.T), rtol=1e-12)

    def test_1455_numpy_semantics(self):
        """#1455: on master np.mean(B.x) returns a dict of per-patch means, np.sum adds
        ghost nodes and duplicated shared nodes, == compares by hand-written rules. New:
        np.mean is the area-weighted mean (a float), np.sum raises pointing to integral(),
        == is identity, elementwise equality is np.isclose (a field of bool)."""
        n, _ = _field(lambda x, y: 2.0 + 0 * x, seed=5)
        self.assertIsInstance(np.mean(n), float)
        self.assertAlmostEqual(np.mean(n), 2.0, places=12)
        with self.assertRaisesRegex(TypeError, "integral"):
            np.sum(n)
        m = n.copy()
        self.assertFalse(n == m)
        self.assertTrue(n == n)
        close = np.isclose(n, m)
        self.assertIsInstance(close, ScalarField)
        self.assertTrue(np.all(close))
        with self.assertRaises(TypeError):
            bool(n)

    def test_1456_region_loading_reads_only_the_region(self):
        """#1456: master reads every patch of every level before any selection (h5py
        datasets of the whole file, or whole-hierarchy merges). New: box= reads only the
        patches meeting the box, at every level."""
        n, be = _field(lambda x, y: x * y, seed=12, domain_cells=(64, 64), patch_cells=4)
        box = PhysBox((2.0, 3.0), (4.0, 5.0))
        n.load(box=box.to_dict())
        total = read = 0
        for g in n.snapshot().hier.levels:
            counts = be.reads.get((0.0, g.level), np.zeros(g.npatch, np.int64))
            was_read = counts[g.file_index] > 0
            selected = np.zeros(g.npatch, bool)
            selected[g.patches_in(box)] = True
            np.testing.assert_array_equal(was_read, selected)
            total += g.npatch
            read += int(was_read.sum())
        self.assertLess(read, total / 4)

    def test_1457_sampling_at_points_uses_the_finest_level(self):
        """#1457: master samples points through merged=True interpolators built from every
        patch (nearest/linear on scattered finest data, 2D at most). New: f(x=xs, y=ys)
        interpolates multilinearly in the finest patch holding each point, reading only
        those patches; exact on linear fields; line cuts return (coords, values)."""
        def lin(x, y):
            return 0.5 + 3.0 * x - 2.0 * y

        n, be = _field(lin, seed=4)
        rng = np.random.default_rng(1)
        xs, ys = rng.uniform(0, n.domain[0], 100), rng.uniform(0, n.domain[1], 100)
        np.testing.assert_allclose(n(x=xs, y=ys), lin(xs, ys), rtol=1e-12)
        coords, values = n(y=1.3)
        np.testing.assert_allclose(values, lin(coords, 1.3), rtol=1e-12)
        self.assertTrue(np.all(np.diff(coords) > 0))
        # a point inside a level-2 patch is answered by level 2
        g2 = n.snapshot().hier.levels[2]
        p = (g2.lower[0] + 0.5) * g2.dl
        before = be.reads.get((0.0, 2), np.zeros(g2.npatch)).sum()
        n(points=p[None, :])
        self.assertGreaterEqual(be.reads[(0.0, 2)].sum(), max(before, 1))

    def test_1458_plot_draws_only_visible_patches(self):
        """#1458: master's plot2d draws one pcolormesh and one Rectangle per patch of every
        level, visible or not (hierarchy.py:498-529). New: one image, and per level one
        LineCollection holding only the outlines of patches inside the view."""
        n, _ = _field(lambda x, y: x, seed=8, domain_cells=(64, 64), patch_cells=4)
        xlim, ylim = (2.0, 4.0), (5.0, 7.0)
        _, ax = n.plot(xlim=xlim, ylim=ylim, plot_patches=True)
        self.assertEqual(len(ax.get_images()), 1)
        lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
        hier = n.snapshot().hier
        self.assertEqual(len(lcs), len(hier.levels))
        drawn = sum(len(lc.get_segments()) for lc in lcs)
        self.assertLess(drawn, sum(g.npatch for g in hier.levels) / 4)
        self.assertGreater(drawn, 0)


if __name__ == "__main__":
    unittest.main()
