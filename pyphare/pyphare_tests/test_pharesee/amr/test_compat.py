#!/usr/bin/env python3
"""Tests of pharesee.amr.compat: master's per-patch hierarchy API on the flat fields.

Real-data tests compare against master's reader (``Run(path, legacy=True).GetB(t)`` ...) on the sample
runs of ``build/sample_data`` and are skipped when those files are absent. The field used
is ``amr.field.AMRField`` when it exists, else ``StandInField`` below (only the members
``compat.MasterAPI`` may use, see FIELD_CONTRACT).
"""

import gc
import os
import sys
import types
import unittest
import weakref

import numpy as np

from pyphare.pharesee.amr import compat
from pyphare.pharesee.amr.compat import LevelView, MasterAPI, PatchView
from pyphare.pharesee.amr.errors import TimeNotFoundError
from pyphare.pharesee.amr.geometry import HierGeometry, LevelGeometry
from pyphare.pharesee.amr.io import open_backend
from pyphare.pharesee.amr.store import LevelData, Snapshot
from pyphare.pharesee.amr.testing import FakeBackend, fill, make_amr, make_large_amr
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES, SAMPLES_MPI3

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS_1D_2D = ("1d", "2d_static", "2d_tagged")


def _have(path):
    return os.path.isdir(path)


class StandInField(MasterAPI):
    """Minimal field: the FIELD_CONTRACT members MasterAPI uses, over Snapshots."""

    def __init__(self, snaps, meta=None):
        self._snaps = {float(s.time): s for s in snaps}
        s0 = snaps[0]
        self.meta = meta
        self.names = tuple(s0.names)
        self.quantity = None if meta is None else meta.quantity
        self.interp_order = 1 if meta is None else meta.interp_order
        self.sim = None
        self.ndim = s0.hier.dim
        self.dl = np.asarray(s0.hier.dl0, dtype=np.float64)
        self.domain = s0.hier.domain

    @classmethod
    def from_backend(cls, backend, times=None):
        meta = backend.meta()
        if times is None:
            times = backend.times()
        snaps = []
        for t in np.atleast_1d(times):
            index = backend.index(float(t))
            geoms = [LevelGeometry.for_level(li.level, li.lower, li.upper, meta.dl0,
                                             meta.domain_cells0, meta.periodic) for li in index]
            t_file = float(backend.times()[np.argmin(np.abs(backend.times() - float(t)))])
            lds = tuple(LevelData.from_backend(backend, g, t_file) for g in geoms)
            snaps.append(Snapshot(HierGeometry(geoms), lds))
        return cls(snaps, meta)

    def times(self):
        return np.array(sorted(self._snaps), dtype=np.float64)

    def snapshot(self, time=None):
        if time is None:
            if len(self._snaps) != 1:
                raise ValueError(f"several times {list(self.times())}")
            return next(iter(self._snaps.values()))
        t = np.asarray(sorted(self._snaps))
        i = int(np.argmin(np.abs(t - float(time))))
        if abs(t[i] - float(time)) > max(2 * np.spacing(abs(t[i])), 5e-11):
            raise TimeNotFoundError(f"no time {time}")
        return self._snaps[float(t[i])]


def _real_field_module():
    try:
        from pyphare.pharesee.amr import field
    except ImportError:
        return None
    return field if hasattr(field, "AMRField") else None


def open_field(path, times=None):
    """A field over the file ``path``: the real AMRField if available, else the stand-in."""
    field = _real_field_module()
    backend = open_backend(path)
    if field is not None:
        return field.AMRField.from_backend(backend, times=times)
    return StandInField.from_backend(backend, times)


def _key(box):
    return tuple(box.lower.tolist()) + tuple(box.upper.tolist())


def f32_ulp(ref):
    """One float32 ulp at each value of ``ref`` (float64 array)."""
    return np.spacing(np.abs(ref).astype(np.float32)).astype(np.float64)


def _interior(a, ng):
    return a[tuple(slice(ng, -ng if ng else None) for _ in range(a.ndim))]


def _ring_mask(shape, ng):
    m = np.ones(shape, bool)
    m[tuple(slice(ng, -ng if ng else None) for _ in shape)] = False
    return m


@unittest.skipUnless(_have(SAMPLES), "needs build/sample_data")
class TestAgainstMaster(unittest.TestCase):
    """(a) same shapes, coordinates and values as master's all_primal objects."""

    @classmethod
    def setUpClass(cls):
        from pyphare.pharesee.run import Run

        cls.Run = Run

    def _compare(self, run, getter, fname, components, moment):
        path = os.path.join(SAMPLES, run, "h5")
        backend = open_backend(os.path.join(path, fname))
        t = float(backend.times()[-1])
        master = getattr(self.Run(path, legacy=True), getter)(t)
        new = open_field(os.path.join(path, fname), times=[t])
        self.assertEqual(new.levelNbrs(), master.levelNbrs())
        self.assertEqual(new.quantities(), master.quantities())
        npatch = 0
        for ilvl, mlvl in master.levels().items():
            nlvl = new.levels()[ilvl]
            self.assertIsInstance(nlvl, LevelView)
            self.assertEqual(len(nlvl.patches), len(mlvl.patches))
            np.testing.assert_array_equal(nlvl.cell_width, mlvl.cell_width)
            by_box = {_key(p.box): p for p in nlvl.patches}
            for mp in mlvl.patches:
                npatch += 1
                p = by_box[_key(mp.box)]
                self.assertEqual(p.box, mp.box)
                np.testing.assert_allclose(p.origin, mp.origin, rtol=0, atol=1e-12)
                np.testing.assert_array_equal(p.dl, mp.dl)
                self.assertEqual(list(p.patch_datas.keys()), list(mp.patch_datas.keys()))
                self.assertEqual(list(p.patch_datas), components)
                self.assertEqual(p.id, mp.id)
                for k in ("lower", "upper", "nbrCells", "origin", "mpi_rank"):
                    np.testing.assert_allclose(p.attrs[k], mp.attrs[k], rtol=0, atol=1e-12)
                # census: layout attributes
                self.assertEqual(p.layout.box, mp.layout.box)
                np.testing.assert_array_equal(p.layout.dl, mp.layout.dl)
                self.assertEqual(p.layout.nbrGhosts(), mp.layout.nbrGhosts())
                self.assertEqual(p.layout.centering["X"]["Bx"], mp.layout.centering["X"]["Bx"])
                self.assertEqual(p.layout.interp_order, mp.layout.interp_order)
                for c in components:
                    self._compare_pdata(p.patch_datas[c], mp.patch_datas[c], moment)
        self.assertGreater(npatch, 0)
        return npatch

    def _compare_pdata(self, pd, mpd, moment):
        ref = np.asarray(mpd.dataset[:], dtype=np.float64)
        got = pd.dataset
        ng = int(pd.ghosts_nbr[0])
        self.assertEqual(got.shape, ref.shape)
        np.testing.assert_array_equal(pd.ghosts_nbr, mpd.ghosts_nbr)
        np.testing.assert_array_equal(pd.size, mpd.size)
        self.assertEqual(pd.ghost_box, mpd.ghost_box)
        self.assertEqual(pd.centerings, ["primal"] * pd.ndim)
        self.assertEqual(pd.field_name, mpd.field_name)
        for axis in "xyz"[: pd.ndim]:
            np.testing.assert_allclose(getattr(pd, axis), getattr(mpd, axis), rtol=0,
                                       atol=1e-12)
        inner, inner_ref = _interior(got, ng), _interior(ref, ng)
        self.assertGreater(inner.size, 0)
        self.assertFalse(np.isnan(inner).any())
        diff = np.abs(inner.astype(np.float64) - inner_ref)
        self.assertTrue((diff <= f32_ulp(inner_ref)).all(),
                        f"{pd.field_name}: max diff {diff.max()} > 1 float32 ulp")
        ring = _ring_mask(got.shape, ng)
        self.assertTrue(np.isnan(got[ring]).all())
        if moment:
            # documented difference: master's moments carry real ghost values
            self.assertTrue(np.isfinite(ref[ring]).all())
        else:
            np.testing.assert_array_equal(np.isnan(got), np.isnan(ref))
        # master's own box selection gives the same interior on both objects
        sel, sel_ref = pd[pd.box], np.asarray(mpd[mpd.box], dtype=np.float64)
        self.assertEqual(sel.shape, sel_ref.shape)
        self.assertTrue((np.abs(sel - sel_ref) <= f32_ulp(sel_ref)).all())

    def test_B_and_E(self):
        for run in RUNS_1D_2D:
            for getter, fname in (("GetB", "EM_B.h5"), ("GetE", "EM_E.h5")):
                with self.subTest(run=run, getter=getter):
                    self._compare(run, getter, fname, ["x", "y", "z"], moment=False)

    def test_Ni(self):
        for run in RUNS_1D_2D:
            with self.subTest(run=run):
                self._compare(run, "GetNi", "ions_charge_density.h5", ["value"], moment=True)

    def test_global_min_max(self):
        path = os.path.join(SAMPLES, "2d_tagged", "h5")
        t = float(open_backend(os.path.join(path, "EM_B.h5")).times()[-1])
        master = self.Run(path, legacy=True).GetB(t)
        new = open_field(os.path.join(path, "EM_B.h5"), times=[t])
        for c in "xyz":
            ref_min, ref_max = master.global_min(c), master.global_max(c)
            self.assertLessEqual(abs(float(new.global_min(c)) - ref_min),
                                 float(f32_ulp(np.array([ref_min]))[0]))
            self.assertLessEqual(abs(float(new.global_max(c)) - ref_max),
                                 float(f32_ulp(np.array([ref_max]))[0]))
        with self.assertRaises(KeyError):
            new.global_min("Bx")

    def test_hierarchy_attributes(self):
        path = os.path.join(SAMPLES, "2d_tagged", "h5")
        t = float(open_backend(os.path.join(path, "EM_B.h5")).times()[-1])
        master = self.Run(path, legacy=True).GetB(t)
        new = open_field(os.path.join(path, "EM_B.h5"), times=[t])
        self.assertEqual(new.domain_box, master.domain_box)
        self.assertEqual(new.ndim, master.ndim)
        self.assertEqual(new.levelNbr(), master.levelNbr())
        self.assertEqual(new.finest_level(), master.finest_level())
        self.assertEqual(new.nbrPatches(), master.nbrPatches())
        self.assertEqual(new.is_homogeneous(), master.is_homogeneous())
        for ilvl in master.levelNbrs():
            self.assertEqual(new.level_domain_box(ilvl), master.level_domain_box(ilvl))
            self.assertEqual(new.level(ilvl).level_range(), master.level(ilvl).level_range())
        self.assertEqual(new.data_files, {})
        self.assertFalse(new.has_file("EM_B.h5"))
        self.assertEqual(new.refinement_ratio, 2)


def _public(cls):
    return {n for n in dir(cls) if not n.startswith("_")}


class TestApiSurface(unittest.TestCase):
    """Every public name of master's objects exists on the views (minus the listed ones)."""

    def test_surface(self):
        from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy
        from pyphare.pharesee.hierarchy.patch import Patch
        from pyphare.pharesee.hierarchy.patchlevel import PatchLevel

        # FIELD's (field.py / plot layer), particles only, or instance attributes of master
        # resolve_selection_boxes: a helper of master's constructor
        not_compat = {"times", "nbytes", "plot", "dist_plot", "sim", "resolve_selection_boxes"}
        missing = _public(PatchHierarchy) - _public(MasterAPI) - not_compat
        self.assertEqual(missing, set())
        self.assertEqual(_public(PatchLevel) - _public(LevelView) - {"patches",
                                                                       "level_number"}, set())
        self.assertEqual(_public(Patch) - _public(PatchView), set())
        for name in ("time_hier", "selection_box", "domain_box", "data_files"):
            self.assertTrue(hasattr(MasterAPI, name))


@unittest.skipUnless(_have(SAMPLES), "needs build/sample_data")
class TestVtkhdfVsH5(unittest.TestCase):
    """(b) the views are the same on both formats: boxes by index, shapes, coords, values."""

    def _compare(self, root, fname, values=True):
        h5 = open_field(os.path.join(root, "h5", fname + ".h5"))
        vtk = open_field(os.path.join(root, "vtk", fname + ".vtkhdf"))
        np.testing.assert_allclose(h5.times(), vtk.times(), rtol=0, atol=1e-6)
        npatch = 0
        for t_h5, t_vtk in zip(h5.times(), vtk.times()):
            lh, lv = h5.levels(t_h5), vtk.levels(t_vtk)
            self.assertEqual(list(lh), list(lv))
            for ilvl in lh:
                ph, pv = lh[ilvl].patches, lv[ilvl].patches
                self.assertEqual(len(ph), len(pv))
                for i in range(len(ph)):
                    npatch += 1
                    self.assertEqual(ph[i].box, pv[i].box, f"t={t_h5} L{ilvl} patch {i}")
                    np.testing.assert_allclose(ph[i].origin, pv[i].origin, atol=1e-6)
                    for c in h5.names:
                        a, b = ph[i].patch_datas[c], pv[i].patch_datas[c]
                        self.assertEqual(a.dataset.shape, b.dataset.shape)
                        np.testing.assert_array_equal(a.ghosts_nbr, b.ghosts_nbr)
                        for axis in "xyz"[: a.ndim]:
                            np.testing.assert_allclose(getattr(a, axis), getattr(b, axis),
                                                       rtol=1e-6, atol=1e-6)
                        np.testing.assert_array_equal(np.isnan(a.dataset), np.isnan(b.dataset))
                        if not values:
                            continue
                        # float32 rounding of the patch's largest value; the 1e-12 floor
                        # absorbs components that are round-off noise (~1e-23 at t=0)
                        scale = float(np.nanmax(np.abs(a.dataset)))
                        np.testing.assert_allclose(a.dataset, b.dataset, rtol=0,
                                                   atol=2.4e-7 * scale + 1e-12)
        self.assertGreater(npatch, 0)

    def test_runs(self):
        for run in RUNS_1D_2D + ("3d",):
            for fname in ("EM_B", "ions_charge_density"):
                with self.subTest(run=run, q=fname):
                    self._compare(os.path.join(SAMPLES, run), fname)

    @unittest.skipUnless(_have(SAMPLES_MPI3), "needs build/sample_data_mpi3")
    def test_three_ranks(self):
        # the two formats come from two runs on 3 ranks, which are not bitwise identical
        # (issue #1526): compare boxes, shapes and coordinates, not values
        self._compare(os.path.join(SAMPLES_MPI3, "2d_tagged"), "EM_B", values=False)
        h5 = open_field(os.path.join(SAMPLES_MPI3, "2d_tagged", "h5", "EM_B.h5"))
        t = h5.times()[-1]
        ranks = {p.id.split("#")[0] for lvl in h5.levels(t).values() for p in lvl.patches}
        self.assertGreater(len(ranks), 1)
        ids = [p.id for p in h5.level(0, t).patches]
        keys = [tuple(int(v) for v in s[1:].split("#")) for s in ids]
        self.assertEqual(keys, sorted(keys), "h5 file order = (rank, localid)")


@unittest.skipUnless(_have(SAMPLES), "needs build/sample_data")
class TestMasterHelpers(unittest.TestCase):
    """(c) master's hierarchy helpers that only walk levels/patches work on the views."""

    @classmethod
    def setUpClass(cls):
        from pyphare.pharesee.run import Run

        cls.path = os.path.join(SAMPLES, "2d_tagged", "h5")
        cls.t = float(open_backend(os.path.join(cls.path, "EM_B.h5")).times()[-1])
        cls.master = Run(cls.path, legacy=True).GetB(cls.t)
        cls.new = open_field(os.path.join(cls.path, "EM_B.h5"), times=[cls.t])

    def test_is_root_lvl(self):
        from pyphare.pharesee.hierarchy.hierarchy_utils import is_root_lvl

        self.assertTrue(is_root_lvl(self.new.level(0)))
        self.assertFalse(is_root_lvl(self.new.level(1)))

    def test_flat_finest_field(self):
        from pyphare.pharesee.hierarchy.hierarchy_utils import flat_finest_field

        for run in RUNS_1D_2D:
            path = os.path.join(SAMPLES, run, "h5")
            from pyphare.pharesee.run import Run

            t = float(open_backend(os.path.join(path, "EM_B.h5")).times()[-1])
            master = Run(path, legacy=True).GetB(t)
            new = open_field(os.path.join(path, "EM_B.h5"), times=[t])
            d_ref, c_ref = flat_finest_field(master, "y", time=t)
            d_new, c_new = flat_finest_field(new, "y", time=t)
            # same points (patch order may differ: compare sorted by coordinates)
            def order(c, d):
                c2 = np.asarray(c).reshape(len(d), -1)
                return np.lexsort(c2.T[::-1])

            o_ref, o_new = order(c_ref, d_ref), order(c_new, d_new)
            np.testing.assert_allclose(np.asarray(c_new)[o_new], np.asarray(c_ref)[o_ref],
                                       atol=1e-12)
            a, b = d_new[o_new].astype(np.float64), np.asarray(d_ref)[o_ref]
            np.testing.assert_array_equal(np.isnan(a), np.isnan(b))
            ok = ~np.isnan(b)
            self.assertTrue((np.abs(a[ok] - b[ok]) <= f32_ulp(b[ok])).all())

    def test_hierarchy_overlaps(self):
        from pyphare.pharesee.geometry import hierarchy_overlaps

        ref = hierarchy_overlaps(self.master, self.t)
        got = hierarchy_overlaps(self.new, self.t)
        self.assertEqual(sorted(ref), sorted(got))
        for ilvl in ref:
            self.assertEqual(len(got[ilvl]), len(ref[ilvl]))
            self.assertEqual(sorted(_key(o["box"]) + (o["name"],) for o in got[ilvl]),
                             sorted(_key(o["box"]) + (o["name"],) for o in ref[ilvl]))

    def test_level_ghost_boxes(self):
        """geometry.level_ghost_boxes deep-copies patches to shift them: master objects."""
        from pyphare.pharesee.geometry import level_ghost_boxes
        from pyphare.pharesee.run import Run

        # fresh objects: other master helpers mutate the hierarchies they are given
        t = 0.1  # this dump has level ghost boxes on level 2 (the last one has none)
        master = Run(self.path, legacy=True).GetB(t)
        new = open_field(os.path.join(self.path, "EM_B.h5"), times=[t])
        ref = level_ghost_boxes(master, "x", time=t)
        got = level_ghost_boxes(new, "x", time=t)
        self.assertEqual(sorted(ref), sorted(got))

        def boxes(res, ilvl):
            return sorted(_key(b) for e in res[ilvl]["x"] for b in e["boxes"])

        for ilvl in ref:
            self.assertEqual(boxes(got, ilvl), boxes(ref, ilvl))
        self.assertGreater(sum(len(boxes(ref, ilvl)) for ilvl in ref), 0)

    def test_compute_hier_from(self):
        """compute_hier_from builds a master hierarchy from the views (patch by index)."""
        from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

        def double(patch_datas, **kwargs):
            pd = patch_datas["x"]
            return ({"name": "x2", "data": 2 * np.asarray(pd.dataset),
                     "centering": pd.centerings, "ghosts_nbr": pd.ghosts_nbr},)

        h = compute_hier_from(double, self.new)
        p_new = self.new.level(1).patches[0]
        p_h = h.level(1, h.times()[0]).patches[0]
        self.assertEqual(p_h.box, p_new.box)
        np.testing.assert_array_equal(p_h.patch_datas["x2"].dataset,
                                      2 * p_new.patch_datas["x"].dataset)

    def test_patch_call_and_plot2d(self):
        p = self.new.level(0).patches[0]
        pd = p.patch_datas["z"]
        y = float(pd.y[pd.ghosts_nbr[1] + 3])
        cut = p("z", y=y)
        ng = int(pd.ghosts_nbr[0])
        np.testing.assert_array_equal(cut, pd.dataset[ng:-ng, ng + 3])
        import matplotlib

        matplotlib.use("Agg")
        fig, ax = self.new.plot2d(qty="z")
        self.assertEqual(len(ax.collections), self.new.nbrPatches())
        import matplotlib.pyplot as plt

        plt.close(fig)


class TestIsinstanceAndOperators(unittest.TestCase):
    """(d) virtual subclasses and mixing legacy with new objects."""

    def test_registration_code(self):
        """hierarchy._register_amr_fields registers whatever amr.field defines."""
        from pyphare.pharesee import amr, hierarchy
        from pyphare.pharesee.hierarchy import PatchHierarchy, ScalarField, VectorField
        from pyphare.pharesee.hierarchy.tensorfield import AnyTensorField

        name = "pyphare.pharesee.amr.field"
        saved_mod = sys.modules.get(name)
        saved_attr = getattr(amr, "field", None)
        fake = types.ModuleType(name)

        class AMRField:
            pass

        class FakeScalar(AMRField):
            pass

        class FakeVector(AMRField):
            pass

        class FakeTensor(AMRField):
            pass

        fake.AMRField, fake.ScalarField = AMRField, FakeScalar
        fake.VectorField, fake.TensorField = FakeVector, FakeTensor
        sys.modules[name] = fake
        amr.field = fake
        try:
            hierarchy._register_amr_fields()
        finally:
            if saved_mod is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = saved_mod
            if saved_attr is None:
                del amr.field
            else:
                amr.field = saved_attr
        v, s = FakeVector(), FakeScalar()
        self.assertIsInstance(v, PatchHierarchy)
        self.assertIsInstance(v, VectorField)
        self.assertIsInstance(v, AnyTensorField)
        self.assertNotIsInstance(v, ScalarField)
        self.assertIsInstance(s, ScalarField)
        self.assertIsInstance(FakeTensor(), AnyTensorField)
        self.assertNotIsInstance(FakeTensor(), VectorField)

    def test_real_field_registered(self):
        field = _real_field_module()
        if field is None:
            self.skipTest("amr.field not written yet")
        from pyphare.pharesee.hierarchy import PatchHierarchy, ScalarField, VectorField

        self.assertTrue(issubclass(field.AMRField, PatchHierarchy))
        self.assertTrue(issubclass(field.VectorField, VectorField))
        self.assertTrue(issubclass(field.ScalarField, ScalarField))

    @unittest.skipUnless(_have(SAMPLES), "needs build/sample_data")
    def test_legacy_operators_give_way(self):
        from pyphare.pharesee.hierarchy import ScalarField, VectorField
        from pyphare.pharesee.hierarchy.tensorfield import is_amr_field
        from pyphare.pharesee.run import Run

        path = os.path.join(SAMPLES, "1d", "h5")
        t = float(open_backend(os.path.join(path, "EM_B.h5")).times()[0])
        run = Run(path, legacy=True)
        B, N = run.GetB(t), run.GetNi(t)
        newB = open_field(os.path.join(path, "EM_B.h5"), times=[t])
        newN = open_field(os.path.join(path, "ions_charge_density.h5"), times=[t])
        self.assertTrue(is_amr_field(newB))
        self.assertFalse(is_amr_field(B))
        field = _real_field_module()
        if field is not None:
            self.assertNotIsInstance(B, field.AMRField)
            self.assertNotIsInstance(N, field.AMRField)
        for op in ("__add__", "__sub__", "__mul__", "__truediv__", "__radd__", "__rmul__",
                   "__rtruediv__"):
            self.assertIs(getattr(ScalarField, op)(N, newN), NotImplemented, op)
        for op in ("__add__", "__sub__", "__mul__", "__truediv__", "__rmul__"):
            self.assertIs(getattr(VectorField, op)(B, newB), NotImplemented, op)
        # legacy + new: legacy gives way, the new object refuses (TypeError, not an assert)
        for a, b in ((N, newN), (B, newB), (B, newN)):
            with self.assertRaises(TypeError):
                a + b
            with self.assertRaises(TypeError):
                a * b
        if field is None:
            return  # new + legacy needs field.py's operators
        for a, b in ((newN, N), (newB, B)):
            with self.assertRaises(TypeError):
                a + b


def _fake_field(npatch, dim=2, width=8):
    hier = make_large_amr(dim, npatch, nlevels=1, width=width)
    data = fill(hier, lambda *c: sum((i + 1) * x for i, x in enumerate(c)),
                dtype=np.float32)
    backend = FakeBackend(data)
    field = _real_field_module()
    if field is not None:
        return field.AMRField.from_backend(backend), backend, data
    return StandInField.from_backend(backend), backend, data


class TestLaziness(unittest.TestCase):
    """(e) levels() builds no patch object and reads nothing; one patch costs one read."""

    def test_levels_build_nothing(self):
        field, backend, data = _fake_field(19000)
        n0 = len(field.level(0).patches)
        self.assertGreaterEqual(n0, 19000)
        created = PatchView.created
        fdv = compat.field_data_view_class()
        created_fd = fdv.created
        levels = field.levels()
        lvl = levels[0]
        self.assertEqual(len(lvl.patches), n0)
        self.assertEqual(PatchView.created, created)
        self.assertEqual(backend.calls, 0)
        p = lvl.patches[5]
        self.assertEqual(PatchView.created, created + 1)
        pd = p.patch_datas["value"]
        self.assertEqual(fdv.created, created_fd + 1)
        self.assertEqual(backend.calls, 0, "building the FieldData view reads nothing")
        values = pd.dataset
        self.assertEqual(PatchView.created, created + 1)
        # the first dataset of a level reads the whole level in ONE call (it fits)
        self.assertEqual(backend.calls, 1)
        np.testing.assert_array_equal(backend.requests[0][2], np.arange(n0))
        # the value read is file-order patch 5 of the source data
        src = data[0]
        ref = src.geom.view(src.buffer("value"), src.geom.store_index[5])
        ng = int(pd.ghosts_nbr[0])
        np.testing.assert_array_equal(values[ng:-ng, ng:-ng], ref)
        self.assertIs(lvl.patches[5], p, "a live view is reused")
        self.assertIs(field.levels()[0], lvl, "levels are cached per (field, time)")
        # slicing builds only the requested views
        created = PatchView.created
        few = lvl.patches[10:13]
        self.assertEqual(len(few), 3)
        self.assertEqual(PatchView.created, created + 3)
        self.assertEqual(backend.calls, 1)
        # every other dataset of the level is already in memory: no more read
        shapes = [q.patch_datas["value"].dataset.shape for q in lvl.patches]
        self.assertEqual(len(shapes), n0)
        self.assertEqual(backend.calls, 1)
        # structure-only accesses never read
        field2, backend2, _ = _fake_field(300)
        for q in field2.level(0).patches:
            q.box, q.origin, q.layout, list(q.patch_datas), q.patch_datas["value"].x
            q.patch_datas["value"].ghost_box, q.patch_datas["value"].size
        self.assertEqual(backend2.calls, 0)

    def test_level_that_does_not_fit_reads_one_patch(self):
        field, backend, _ = _fake_field(400)

        class TinyBudget:  # a level never fits; single patches always do
            def available(self):
                return 1024

            def check(self, needed, what, force=False):
                return "ok"

        for ld in field.levels_data():
            ld.budget = TinyBudget()
        lvl = field.level(0)
        for i in (7, 9):
            self.assertEqual(lvl.patches[i].patch_datas["value"].dataset.ndim, 2)
        self.assertEqual(backend.calls, 2)
        self.assertEqual([list(r[2]) for r in backend.requests], [[7], [9]])

    def test_views_die_with_field(self):
        field, _, _ = _fake_field(200)
        lvl_ref = weakref.ref(field.level(0))
        p_ref = weakref.ref(field.level(0).patches[3])
        gc.collect()
        self.assertIsNone(p_ref(), "nobody holds the patch view: it is gone")
        self.assertIsNotNone(lvl_ref(), "the field caches its level views")
        f_ref = weakref.ref(field)
        del field
        gc.collect()
        self.assertIsNone(f_ref())
        self.assertIsNone(lvl_ref(), "level views die with the field")

    def test_held_view_survives_field(self):
        field, _, _ = _fake_field(50)
        p = field.level(0).patches[0]
        del field
        gc.collect()
        self.assertEqual(p.patch_datas["value"].dataset.ndim, 2)


class TestDataset(unittest.TestCase):
    """(f) dataset is a read-only copy, independent from the store; edge layouts."""

    def _memory_field(self, dim, two_node=False, interp_order=1):
        hier = make_amr(dim, nlevels=2, seed=3, two_node=two_node)
        lds = fill(hier, lambda *c: 1.0 + sum(c), dtype=np.float64)
        field = StandInField([Snapshot(hier, tuple(lds))])
        field.interp_order = interp_order
        return field, lds

    def test_read_only_and_independent(self):
        field, lds = self._memory_field(2)
        pd = field.level(0).patches[0].patch_datas["value"]
        d = pd.dataset
        with self.assertRaises(ValueError):
            d[0, 0] = 1.0
        with self.assertRaises(AttributeError):
            pd.dataset = np.zeros(3)
        with self.assertRaises(ValueError):
            pd[pd.box] = 0.0
        buf = lds[0].buffer("value")
        self.assertFalse(np.shares_memory(d, buf))
        before = d.copy()
        buf[:] = -5.0  # the in-memory store is writeable: the copy does not follow
        np.testing.assert_array_equal(pd.dataset, before)
        import copy

        c = copy.deepcopy(pd)
        c.dataset[...] = 0.0  # writeable master FieldData
        np.testing.assert_array_equal(pd.dataset, before)

    def test_shapes_and_coords_1d_2d_3d(self):
        for dim in (1, 2, 3):
            for interp, ng in ((1, 2), (2, 4), (3, 4), (None, 2)):
                with self.subTest(dim=dim, interp=interp):
                    field, _ = self._memory_field(dim, two_node=True, interp_order=interp)
                    for lvl in field.levels().values():
                        for p in lvl.patches:
                            pd = p.patch_datas["value"]
                            shape = tuple(p.box.shape + 1 + 2 * ng)
                            self.assertEqual(pd.dataset.shape, shape)
                            self.assertEqual(tuple(pd.size), shape)
                            self.assertEqual(int(pd.ghosts_nbr[0]), ng)
                            inner = _interior(pd.dataset, ng)
                            self.assertEqual(inner.shape, tuple(p.box.shape + 1))
                            self.assertFalse(np.isnan(inner).any())
                            self.assertTrue(np.isnan(pd.dataset[_ring_mask(shape, ng)]).all())
                            # f = 1 + sum of coordinates, coordinates include the ring
                            grids = pd.meshgrid()
                            ref = 1.0 + sum(grids)
                            np.testing.assert_allclose(_interior(pd.dataset, ng),
                                                       _interior(ref, ng), atol=1e-12)
                            x = pd.x
                            self.assertEqual(len(x), shape[0])
                            np.testing.assert_allclose(x[ng], p.origin[0])
            # two-node patches exist in this layout
        field, _ = self._memory_field(1, two_node=True)
        self.assertIn(1, [int(p.box.shape[0]) for p in field.level(0).patches])

    def test_patch_api(self):
        field, _ = self._memory_field(2)
        lvl = field.level(1)
        self.assertEqual(lvl.level_number, 1)
        p = lvl[0]
        self.assertIs(p, lvl.patches[0])
        self.assertIs(lvl.patches[-1], lvl.patches[len(lvl.patches) - 1])
        with self.assertRaises(IndexError):
            lvl.patches[len(lvl.patches)]
        with self.assertRaises(IndexError):
            lvl["a"]
        self.assertEqual(sum(1 for _ in lvl), len(lvl.patches))
        self.assertEqual(p.id, "mem#0")  # in memory: no rank, not parseable as one
        self.assertEqual(set(p.attrs), {"lower", "upper", "nbrCells", "origin"})
        np.testing.assert_array_equal(p.attrs["nbrCells"], p.box.shape)
        self.assertIs(p["value"], p.patch_datas["value"])
        with self.assertRaises(KeyError):
            p["x"]
        self.assertIs(p.layout, p.patch_datas["value"].layout)
        q = p.copy()
        np.testing.assert_array_equal(q.patch_datas["value"].dataset,
                                      p.patch_datas["value"].dataset)


class TestMasterKeysAndViews(unittest.TestCase):
    """Component keys (A7), NaN-aware compare (A3), pickling the views (A10)."""

    def _vector(self):
        hier = make_amr(2, nlevels=2, seed=4)
        lds = fill(hier, {"x": lambda x, y: x, "y": lambda x, y: y, "z": lambda x, y: x + y})
        field = _real_field_module()
        if field is None:
            self.skipTest("needs amr.field")
        return field.make_field([Snapshot(hier, tuple(lds))], quantity="EM_B")

    def test_component_answers_parent_key(self):
        B = self._vector()
        Bx = B.x
        p = Bx.level(0).patches[0]
        self.assertEqual(list(p.patch_datas.keys()), ["x"])  # master: sub-hierarchy key
        self.assertEqual(Bx.quantities(), ["x"])
        self.assertIs(p.patch_datas["x"], p.patch_datas["value"])
        self.assertIn("x", p.patch_datas)
        self.assertIn("value", p.patch_datas)
        self.assertNotIn("y", p.patch_datas)
        self.assertEqual(p.patch_datas["x"].field_name, "x")
        np.testing.assert_array_equal(p["x"].dataset,
                                      B.level(0).patches[0].patch_datas["x"].dataset)
        self.assertEqual(Bx.global_max("x"), Bx.global_max("value"))
        self.assertEqual(Bx.global_min("x"), B.global_min("x"))
        with self.assertRaises(KeyError):
            Bx.global_max("y")
        with self.assertRaises(KeyError):
            p.patch_datas["y"]
        with self.assertRaises(KeyError):  # master: B.x["x"] is a KeyError too
            Bx["x"]
        self.assertEqual(B["x"].quantities(), ["x"])
        self.assertEqual(B["Bx"].quantities(), ["x"])
        # results of operations are new quantities: "value" only
        self.assertEqual((Bx * 2).quantities(), ["value"])
        self.assertEqual(B.quantities(), ["x", "y", "z"])

    def test_compare_treats_nan_as_equal(self):
        from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare

        B = self._vector()
        a = B.level(0).patches[0].patch_datas["x"]
        b = B.copy().level(0).patches[0].patch_datas["x"]
        self.assertTrue(np.isnan(a.dataset).any())  # the NaN ghost ring
        self.assertTrue(a.compare(b))
        self.assertTrue(a == b)
        self.assertTrue(bool(hierarchy_compare(B, B.copy())))
        self.assertFalse(bool(hierarchy_compare(B, B * 2.0)))
        c = (B * 2.0).level(0).patches[0].patch_datas["x"]
        self.assertFalse(a.compare(c))

    def test_pickling_views_says_pickle_the_field(self):
        import copy
        import pickle

        B = self._vector()
        lvl = B.level(0)
        for view in (lvl, lvl.patches[0]):
            with self.assertRaisesRegex(TypeError, "Pickle the field itself"):
                pickle.dumps(view)
        self.assertEqual(copy.copy(lvl).level_number, 0)
        self.assertGreater(len(pickle.dumps(B)), 0)


@unittest.skipUnless(_have(SAMPLES_MPI3), "needs build/sample_data_mpi3")
class TestRealIdsAndCompare(unittest.TestCase):
    def test_vtkhdf_patch_id_is_not_a_rank(self):
        from pyphare.pharesee.run import Run

        root = os.path.join(SAMPLES_MPI3, "2d_tagged")
        Ni = Run(os.path.join(root, "vtk")).GetNi(0.0)
        ids = [p.id for lvl in Ni.levels().values() for p in lvl.patches]
        self.assertEqual(ids[:3], ["vtk#0", "vtk#1", "vtk#2"])
        # master's rank-from-id idiom (run/utils.py:418) must fail, not answer rank 0
        with self.assertRaises(ValueError):
            int(ids[0].split("#")[0][1:])
        with self.assertRaises(ValueError):
            int(ids[0].strip("p").split("#")[0])
        # phareh5 keeps the file's group names (real ranks)
        h5 = Run(os.path.join(root, "h5"), legacy=False).GetNi(0.0)
        ranks = {int(p.id.split("#")[0][1:]) for lvl in h5.levels().values()
                 for p in lvl.patches}
        self.assertEqual(ranks, {0, 1, 2})

    def test_hierarchy_compare_same_data_is_true(self):
        from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare
        from pyphare.pharesee.run import Run

        root = os.path.join(SAMPLES_MPI3, "2d_tagged")
        for r in (Run(os.path.join(root, "h5"), legacy=False), Run(os.path.join(root, "vtk"))):
            with self.subTest(r=repr(r)):
                self.assertTrue(bool(hierarchy_compare(r.GetNi(0.2), r.GetNi(0.2))))
                self.assertTrue(bool(hierarchy_compare(r.GetB(0.2), r.GetB(0.2))))
                self.assertFalse(bool(hierarchy_compare(r.GetNi(0.2), r.GetNi(0.0))))

    def test_component_keys_as_master_on_real_files(self):
        from pyphare.pharesee.run import Run

        root = os.path.join(SAMPLES_MPI3, "2d_tagged")
        legacy = Run(os.path.join(root, "h5"), legacy=True).GetB(0.2)
        ref = legacy.x.level(0).patches[0].patch_datas
        for r in (Run(os.path.join(root, "h5"), legacy=False), Run(os.path.join(root, "vtk"))):
            B = r.GetB(0.2)
            pd = B.x.level(0).patches[0].patch_datas
            self.assertEqual(list(pd.keys()), list(ref.keys()))
            self.assertEqual(pd["x"].dataset.shape, ref["x"].dataset.shape)
            self.assertEqual(B.x.quantities(), legacy.x.quantities())
            self.assertAlmostEqual(float(B.x.global_max("x")), legacy.x.global_max("x"),
                                   places=5)


class TestTimes(unittest.TestCase):
    """(g) time_hier keys are floats, master's strings work as keys, ambiguity raises."""

    def _two_times(self):
        hier = make_amr(2, nlevels=2, seed=1)
        snaps = [Snapshot(hier, tuple(fill(hier, lambda x, y, k=k: x + k, time=t)))
                 for k, t in enumerate((0.05, 0.1))]
        return StandInField(snaps)

    def test_time_hier(self):
        from pyphare.pharesee.hierarchy.hierarchy import format_timestamp

        f = self._two_times()
        th = f.time_hier
        self.assertEqual(list(th.keys()), [0.05, 0.1])
        self.assertTrue(all(isinstance(k, float) for k in th))
        self.assertEqual(len(th), 2)
        lv = th[format_timestamp(0.1)]
        self.assertIs(lv[0], f.levels(0.1)[0])
        self.assertIs(th[0.1][1], f.level(1, "0.1000000000"))
        self.assertIn(format_timestamp(0.05), th)
        self.assertNotIn(0.07, th)
        self.assertNotIn("garbage", th)
        self.assertTrue(f.has_time("0.0500000000"))
        self.assertTrue(f.has_time(0.1))
        self.assertFalse(f.has_time(0.2))
        d0 = f.levels(0.05)[0].patches[0].patch_datas["value"].dataset
        d1 = f.levels(0.1)[0].patches[0].patch_datas["value"].dataset
        np.testing.assert_allclose(d1 - d0, np.where(np.isnan(d0), np.nan, 1.0))

    def test_ambiguity_raises(self):
        f = self._two_times()
        for call in (f.levels, f.levelNbr, f.levelNbrs, f.finest_level,
                     lambda: f.level(0), lambda: f.global_min("value"), f._default_time):
            with self.assertRaises(ValueError) as cm:
                call()
            self.assertIn("0.05", str(cm.exception))
        with self.assertRaises(TimeNotFoundError):
            f.levels(0.3)
        self.assertEqual(f.nbrPatches(), 2 * sum(len(lvl.patches)
                                                 for lvl in f.levels(0.05).values()))


if __name__ == "__main__":
    unittest.main()
