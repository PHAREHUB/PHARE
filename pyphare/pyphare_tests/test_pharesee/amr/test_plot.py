"""``field.plot()`` on the new fields (views/plot.py): master's keywords, one image, colour
limits of what is shown, reads limited to the zoom, culled patch outlines.

Synthetic fields on random AMR layouts behind ``testing.FakeBackend`` (which counts the
patches read); one smoke test on the real sample run (``build/sample_data``, an error when absent).
matplotlib runs with the Agg backend.
"""

import os
import unittest

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, QuadMesh
from matplotlib.colors import LogNorm

from pyphare.pharesee.amr.field import AMRField
from pyphare.pharesee.amr.testing import make_amr
from pyphare_tests.test_pharesee.amr._common import (
    PlotTestCase,
    fake_backend,
    fake_field,
    tmp_dir,
    vec,
)
from pyphare_tests.test_pharesee.amr.sample_data import require, run_dirs

os.environ.setdefault("PHARESEE_CACHE", "0")


def _scalar(dim=2, nlevels=3, seed=5, f=None, **kw):
    if f is None:
        f = (lambda x, y: np.sin(x) * np.cos(0.7 * y) + 0.1 * x) if dim == 2 else np.sin
    return fake_field(make_amr(dim, nlevels, seed, periodic=True, **kw), f)


def _image(ax):
    images = ax.get_images()
    assert len(images) == 1, images
    return images[0]


class PlotBasicsTest(PlotTestCase):
    def test_returns_fig_ax_and_draws_one_image(self):
        f, _ = _scalar()
        fig, ax = f.plot()
        self.assertIs(ax.figure, fig)
        _image(ax)
        self.assertFalse([c for c in ax.collections if isinstance(c, QuadMesh)])

    def test_given_axes_is_used(self):
        f, _ = _scalar()
        fig0, ax0 = plt.subplots()
        fig, ax = f.plot(ax=ax0, title="t", xlabel="X", ylabel="Y", aspect="auto")
        self.assertIs(ax, ax0)
        self.assertIs(fig, fig0)
        self.assertEqual((ax.get_title(), ax.get_xlabel(), ax.get_ylabel()), ("t", "X", "Y"))

    def test_unknown_keyword_is_a_type_error(self):
        f, _ = _scalar()
        with self.assertRaisesRegex(TypeError, "colour"):
            f.plot(colour="red")

    def test_vector_needs_qty_and_multi_time_needs_time(self):
        comps = vec(lambda x, y: x, lambda x, y: y, lambda x, y: x * y)
        B, _ = fake_field(make_amr(2, 2, 1, periodic=True), lambda t: comps, times=(0.0, 1.0))
        with self.assertRaisesRegex(ValueError, "holds 2 times"):
            B.plot(qty="x")
        with self.assertRaisesRegex(ValueError, "specify a quantity in x y z"):
            B.plot(time=1.0)
        _, ax = B.plot(time=1.0, qty="Bz")
        g = B.at(1.0)["z"].finest()
        self.assertEqual(_image(ax).get_clim(), (float(np.nanmin(g.data)),
                                                 float(np.nanmax(g.data))))

    def test_component_accepts_parent_name_and_derived_title_shows_stencil(self):
        from pyphare.pharesee.amr.ops import deriv

        comps = vec(lambda x, y: np.sin(y), lambda x, y: np.cos(x), lambda x, y: x * y)
        be = fake_backend(make_amr(2, 2, 1, periodic=True), comps)
        B = AMRField.from_backend(be, quantity="EM_B")
        for qty in ("x", "value", "EM_B_x", None):
            with self.subTest(qty=qty):
                _, ax = B.x.plot(qty=qty)
                self.assertEqual(ax.get_title(), "")  # stored quantity: no default title
                plt.close("all")
        J = B._like([deriv.curl(B.snapshot())], quantity="J", stencil="node")
        _, ax = J.plot(qty="z")
        self.assertEqual(ax.get_title(), "J_z (stencil: node)")
        self.assertIn("stencil node", repr(J))
        _, ax = J.plot(qty="z", title="mine")
        self.assertEqual(ax.get_title(), "mine")
        _, ax = J.z.plot(qty="z")
        self.assertEqual(ax.get_title(), "J_z (stencil: node)")


class ColourLimitsTest(PlotTestCase):
    def setUp(self):
        # a smooth field with two isolated spikes: percentiles clip both sides
        def f(x, y):
            v = np.sin(x) * np.cos(0.7 * y)
            v = np.where((np.abs(x - 2.0) < 0.15) & (np.abs(y - 2.0) < 0.15), 40.0, v)
            return np.where((np.abs(x - 6.0) < 0.15) & (np.abs(y - 4.0) < 0.15), -30.0, v)

        self.f, _ = _scalar(f=f)
        self.shown = self.f.finest().data

    def test_default_is_min_max_of_shown_data(self):
        _, ax = self.f.plot()
        im = _image(ax)
        self.assertEqual(im.get_clim(), (float(np.nanmin(self.shown)),
                                         float(np.nanmax(self.shown))))
        self.assertEqual(im.colorbar.extend, "neither")

    def test_robust_uses_percentiles_and_extends_clipped_sides(self):
        _, ax = self.f.plot(robust=True)
        im = _image(ax)
        lo, hi = np.percentile(self.shown[np.isfinite(self.shown)], (1.0, 99.7))
        np.testing.assert_allclose(im.get_clim(), (lo, hi))
        self.assertEqual(im.colorbar.extend, "both")

    def test_symmetric(self):
        _, ax = self.f.plot(symmetric=True)
        m = max(abs(np.nanmin(self.shown)), abs(np.nanmax(self.shown)))
        self.assertEqual(_image(ax).get_clim(), (-m, m))

    def test_explicit_limits_win(self):
        _, ax = self.f.plot(vmin=-0.5, vmax=0.5, robust=True, symmetric=True)
        im = _image(ax)
        self.assertEqual(im.get_clim(), (-0.5, 0.5))
        self.assertEqual(im.colorbar.extend, "both")
        _, ax = self.f.plot(vmax=0.5)
        self.assertEqual(_image(ax).get_clim(), (float(np.nanmin(self.shown)), 0.5))
        self.assertEqual(_image(ax).colorbar.extend, "max")

    def test_norm(self):
        f, _ = _scalar(f=lambda x, y: 1.0 + x * y)
        norm = LogNorm(vmin=1.0, vmax=10.0)
        _, ax = f.plot(norm=norm)
        self.assertIs(_image(ax).norm, norm)
        with self.assertRaisesRegex(ValueError, "not both"):
            f.plot(norm=LogNorm(), vmin=1.0)


class ZoomLevelsPatchesTest(PlotTestCase):
    def _field(self):
        return _scalar(nlevels=3, seed=8, patch_cells=4, domain_cells=(64, 64))

    def test_zoom_reads_only_patches_under_it(self):
        f, be = self._field()
        xlim, ylim = (3.0, 4.6), (7.1, 8.3)
        f.plot(xlim=xlim, ylim=ylim)
        snap = f.snapshot()
        dl0 = snap.hier.levels[0].dl
        L = snap.hier.finest_level
        total = read = 0
        for g in snap.hier.levels:
            counts = be.reads.get((0.0, g.level), np.zeros(g.npatch, np.int64))
            was_read = counts[g.file_index] > 0  # store order
            lo, hi = g.lower * g.dl, (g.upper + 1) * g.dl
            # finest() windows: the zoom plus at most one coarse node per level of margin
            m = (L + 1) * dl0
            near = ((hi > np.array([xlim[0], ylim[0]]) - m)
                    & (lo < np.array([xlim[1], ylim[1]]) + m)).all(axis=1)
            under = ((hi > np.array([xlim[0], ylim[0]]))
                     & (lo < np.array([xlim[1], ylim[1]]))).all(axis=1)
            self.assertFalse((was_read & ~near).any(), f"level {g.level}")
            self.assertTrue(was_read[under].all(), f"level {g.level}")
            total += g.npatch
            read += int(was_read.sum())
        self.assertLess(read, total / 4)

    def test_levels_shows_only_that_level(self):
        f, _ = self._field()
        _, ax = f.plot(levels=[1])
        g = f.finest(max_level=1)
        arr = _image(ax).get_array()  # data.T
        mask = np.ma.getmaskarray(arr).T
        np.testing.assert_array_equal(mask, g.level_map != 1)
        self.assertTrue((g.level_map == 1).any() and (g.level_map == 0).any())
        shown = np.asarray(g.data)[g.level_map == 1]
        self.assertEqual(_image(ax).get_clim(), (float(shown.min()), float(shown.max())))

    def test_patch_outlines_one_collection_per_level_culled_to_view(self):
        f, _ = self._field()
        xlim, ylim = (2.0, 7.0), (3.0, 9.0)
        _, ax = f.plot(xlim=xlim, ylim=ylim, plot_patches=True, patchcolors=["r", "g", "b"])
        lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
        snap = f.snapshot()
        self.assertEqual(len(lcs), len(snap.hier.levels))
        g = f.finest(box=dict(x=xlim, y=ylim))
        view_lo, view_hi = np.array([g.x[0], g.y[0]]), np.array([g.x[-1], g.y[-1]])
        for lc, geom in zip(lcs, snap.hier.levels):
            n = 0
            for slot in range(geom.npatch):  # reference: one patch at a time
                lo = geom.lower[slot] * geom.dl
                hi = (geom.upper[slot] + 1) * geom.dl
                n += bool(((hi > view_lo) & (lo < view_hi)).all())
            self.assertEqual(len(lc.get_segments()), n, f"level {geom.level}")
            self.assertLess(n, geom.npatch)
        colors = [tuple(lc.get_colors()[0][:3]) for lc in lcs]
        self.assertEqual(colors, [(1.0, 0.0, 0.0), (0.0, 0.5, 0.0), (0.0, 0.0, 1.0)])


class Plot1DTest(PlotTestCase):
    def test_one_line_per_level_with_node_values(self):
        f, _ = _scalar(dim=1, nlevels=2, f=lambda x: x * x)
        _, ax = f.plot(levels=(0, 1), marker="o", color="r")
        self.assertEqual(len(ax.lines), 2)
        for line, ld in zip(ax.lines, f.snapshot().levels):
            x, y = line.get_xdata(), line.get_ydata()
            ok = np.isfinite(y)
            np.testing.assert_allclose(y[ok], x[ok] ** 2)
            self.assertEqual(int((~ok).sum()), ld.geom.npatch - 1)  # breaks between patches
            self.assertEqual(line.get_label(), f"L{ld.geom.level}")
        f.plot()  # default levels (0,) as master
        self.assertEqual(len(plt.gca().lines), 1)

    def test_1d_limits_and_save(self):
        f, _ = _scalar(dim=1, nlevels=2, f=np.cos)
        path = str(tmp_dir(self) / "p.png")
        _, ax = f.plot(xlim=(1.0, 3.0), ylim=(-2, 2), legend=True, filename=path)
        self.assertTrue(os.path.isfile(path))
        self.assertEqual(ax.get_xlim(), (1.0, 3.0))
        self.assertIsNotNone(ax.get_legend())


class RealDataPlotTest(PlotTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        require()

    def test_bz_of_the_sample_run(self):
        from pyphare.pharesee.run import Run

        B = Run(run_dirs("2d_tagged")[1]).GetB(0.1)
        _, ax = B.plot(qty="z", plot_patches=True)
        g = B.z.finest()
        self.assertEqual(_image(ax).get_clim(), (float(np.nanmin(g.data)),
                                                 float(np.nanmax(g.data))))
        self.assertEqual(_image(ax).get_array().shape, g.data.T.shape)


if __name__ == "__main__":
    unittest.main()
