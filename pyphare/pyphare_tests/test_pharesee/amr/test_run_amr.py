"""``Run`` routing to the new reader, on the real sample runs.

Data: ``build/sample_data/{1d,2d_static,2d_tagged,3d}/{h5,vtk}`` (the same runs written in
phareh5 and vtkhdf), written by ``gen_sample_data.py`` (a ctest fixture); location: ``sample_data.py``
(``$PHARESEE_SAMPLE_DATA``). A missing run is an error (``sample_data.require``); the tests
of the 3-rank run are in test_mpi3.py.

Tolerances (measured on these files):
  * stored quantities h5 (new reader) vs vtkhdf: the C++ writer and the reader do the same
    2-point Yee-to-node averages in float32; worst 1.19e-7 x max|level| -> bound 2.5e-7.
  * J (node stencil) and divB between formats, at nodes with same-level neighbours: the
    float32 differences of B (1 ulp) divided by dl; bound 1e-5 x max|J| for J and
    1e-5 x max|B| / min(dl) for divB.
  * GetJ(stencil="yee") and the 2D per-cell GetDivB on the new h5 path vs master: equal,
    master's float64 J rounded to the file's dtype (same operations in the same order).
    Node B and the other averaged quantities: _common.assert_matches_master (each reader vs
    the exact mean of the cells, within its own rounding bound; also on float64 copies).
    GetPe: equal, master's float64 Ni * Te rounded to the file's dtype.
  * magnetic flux new vs master: 5 % of max|Az|. Both integrate the finest data on the
    finest level's nodes; master interpolates (nearest) its Yee-staggered merged data, the
    new reader paints node-averaged data with finest() (coarse-fine gaps filled
    differently); measured 1.1-2.2 %.
"""

import os
import tempfile
import unittest
import warnings

import numpy as np

from pyphare.pharesee import hierarchy as legacy
from pyphare.pharesee.amr import field as amr
from pyphare.pharesee.amr.boxes import PhysBox
from pyphare.pharesee.amr.errors import TimeNotFoundError
from pyphare.pharesee.amr.ops.deriv import one_sided_mask
from pyphare.pharesee.run import Run
from pyphare.pharesee.run.run import LegacyRun
from pyphare_tests.test_pharesee.amr._common import (
    assert_all_equal,
    check_against_master,
    float64_copy,
    tmp_dir,
)
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES, require

os.environ.setdefault("PHARESEE_CACHE", "0")  # never write the index cache under ~/.cache

CASES = ("1d", "2d_static", "2d_tagged", "3d")
POPS = ("protons", "alpha")


def _dirs(case, root=SAMPLES):
    return os.path.join(root, case, "h5"), os.path.join(root, case, "vtk")


def _registered():
    """True when the amr classes are registered as master's (hierarchy/__init__.py)."""
    return issubclass(amr.VectorField, legacy.VectorField)


class LegacyEnv:
    """Remove PHARESEE_LEGACY for the duration of a test (the default must be tested)."""

    def setUp(self):
        self._env = os.environ.pop("PHARESEE_LEGACY", None)

    def tearDown(self):
        if self._env is not None:
            os.environ["PHARESEE_LEGACY"] = self._env


def _level_scale(ld):
    return max(float(np.nanmax(np.abs(ld.buffer(n)))) for n in ld.names) or 1.0


def _patches_by_box(field, time=None):
    """{(level, lower, upper): {component: patch array [x, y, z]}} of a new-path field."""
    out = {}
    for ld in field.snapshot(time).levels:
        g = ld.geom
        for slot in range(g.npatch):
            key = (g.level, tuple(int(v) for v in g.lower[slot]),
                   tuple(int(v) for v in g.upper[slot]))
            out[key] = {n: np.asarray(g.view(ld.buffer(n), slot)) for n in ld.names}
    return out


def _hier_arrays(hier, time):
    """{(level, lower, upper): {quantity: array with ghosts}} of a master hierarchy."""
    out = {}
    for il, lvl in hier.levels(time).items():
        for p in lvl.patches:
            key = (il, tuple(int(v) for v in p.box.lower), tuple(int(v) for v in p.box.upper))
            out[key] = {q: np.asarray(pd.dataset[:]) for q, pd in p.patch_datas.items()}
    return out


def _legacy_by_box(hier, comps, time):
    """Same for a master hierarchy: primal nodes without the ghost ring, float64."""
    out = {}
    for il, lvl in hier.levels(time).items():
        for p in lvl.patches:
            entry = {}
            for c in comps:
                pd = p.patch_datas[c]
                ng = int(pd.ghosts_nbr[0])
                sl = tuple(slice(ng, -ng) for _ in range(pd.ndim))
                entry[c] = np.asarray(pd.dataset[:], dtype=np.float64)[sl]
            out[(il, tuple(int(v) for v in p.box.lower), tuple(int(v) for v in p.box.upper))] = (
                entry
            )
    return out


def _assert_rounded_master(tc, new, ref, msg):
    """new (_patches_by_box) == master's float64 values (_legacy_by_box) rounded to new's
    dtype, at every node of every patch."""
    tc.assertEqual(new.keys(), ref.keys())
    for key in new:
        for c in ref[key]:
            np.testing.assert_array_equal(new[key][c], ref[key][c].astype(new[key][c].dtype),
                                          err_msg=f"{msg} {key} {c}")


# ============================================================================ routing


class RoutingTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require()

    def test_reader_chosen_once_and_shown_in_repr(self):
        h5, vtk = _dirs("2d_tagged")
        self.assertEqual(Run(vtk).reader, "new")
        self.assertIn("vtkhdf: new reader", repr(Run(vtk)))
        r = Run(h5)
        self.assertEqual(r.reader, "new")
        self.assertIn("h5: new reader", repr(r))
        r = Run(h5, legacy=True)
        self.assertEqual(r.reader, "legacy")
        self.assertIn("legacy reader", repr(r))
        os.environ["PHARESEE_LEGACY"] = "1"
        try:
            self.assertEqual(Run(h5).reader, "legacy")
            self.assertEqual(Run(h5, legacy=False).reader, "new")
        finally:
            os.environ.pop("PHARESEE_LEGACY")

    def test_legacy_h5_path_is_masters_code(self):
        h5, _ = _dirs("2d_tagged")
        os.environ["PHARESEE_LEGACY"] = "1"
        try:
            runs = (Run(h5, legacy=True), Run(h5))
        finally:
            os.environ.pop("PHARESEE_LEGACY")
        for r in runs:
            self._check_masters_code(r, h5)

    def _check_masters_code(self, r, h5):
        B, Ni = r.GetB(0.1), r.GetNi(0.1)
        self.assertIs(type(B), legacy.VectorField)
        self.assertIs(type(Ni), legacy.ScalarField)
        self.assertEqual(r._backends, {})  # the new reader never opened a file
        lr = LegacyRun(h5)
        np.testing.assert_array_equal(r.times("B"), lr.times("B"))
        self.assertEqual(r.all_times().keys(), lr.all_times().keys())
        np.testing.assert_array_equal(r.GetDl(1), lr.GetDl(1))

    def test_both_formats_need_format(self):
        h5, vtk = _dirs("2d_tagged")
        with tempfile.TemporaryDirectory() as tmp:
            os.symlink(os.path.join(h5, "EM_B.h5"), os.path.join(tmp, "EM_B.h5"))
            os.symlink(os.path.join(vtk, "EM_B.vtkhdf"), os.path.join(tmp, "EM_B.vtkhdf"))
            os.symlink(os.path.join(vtk, "EM_E.vtkhdf"), os.path.join(tmp, "EM_E.vtkhdf"))
            r = Run(tmp)
            self.assertIsInstance(r.GetE(0.1), amr.VectorField)  # unambiguous
            with self.assertRaisesRegex(ValueError, "format="):
                r.GetB(0.1)
            self.assertEqual(Run(tmp, format="vtkhdf").GetB(0.1).format, "vtkhdf")
            self.assertEqual(Run(tmp, format="h5").GetB(0.1).format, "phareh5")
            self.assertIs(type(Run(tmp, format="h5", legacy=True).GetB(0.1)),
                          legacy.VectorField)
            with self.assertRaisesRegex(ValueError, "cannot read .vtkhdf"):
                Run(tmp, format="vtkhdf", legacy=True)

    def test_legacy_true_reads_h5_when_both_formats_exist(self):
        """Master reads h5 only: legacy=True (argument or PHARESEE_LEGACY=1) must not raise
        on a stem present in both formats. Without legacy it still asks for format=."""
        h5, vtk = _dirs("2d_tagged")
        with tempfile.TemporaryDirectory() as tmp:
            os.symlink(os.path.join(h5, "EM_B.h5"), os.path.join(tmp, "EM_B.h5"))
            os.symlink(os.path.join(vtk, "EM_B.vtkhdf"), os.path.join(tmp, "EM_B.vtkhdf"))
            r = Run(tmp, legacy=True)
            self.assertIs(type(r.GetB(0.1)), legacy.VectorField)
            np.testing.assert_array_equal(r.times("B"), LegacyRun(tmp).times("B"))
            os.environ["PHARESEE_LEGACY"] = "1"
            try:
                self.assertIs(type(Run(tmp).GetB(0.1)), legacy.VectorField)
            finally:
                os.environ.pop("PHARESEE_LEGACY")
            with self.assertRaisesRegex(ValueError, "format="):
                Run(tmp).GetB(0.1)

    def test_missing_quantity_lists_directory(self):
        _, vtk = _dirs("2d_tagged")
        with self.assertRaisesRegex(FileNotFoundError, "EM_B"):
            Run(vtk).GetMHDrho(0.1)

    def test_empty_vtkhdf_file_raises_its_format_error(self):
        _, vtk = _dirs("2d_tagged")
        from pyphare.pharesee.amr.errors import FormatError
        from pyphare.pharesee.amr.io import open_backend
        from pyphare.pharesee.amr.io.files import find_quantity_file

        path = find_quantity_file(vtk, "ions_momentum_tensor")
        with self.assertRaisesRegex(FormatError, "does not write momentum_tensor"):
            open_backend(path)


# ============================================================================ getters


class GettersTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require(cases=CASES)

    def _check_classes(self, r, t, divb_master=False):
        """divb_master: GetDivB is master's per-cell object (2D phareh5), not an AMRField."""
        scalar = [r.GetNi(t), r.GetMassDensity(t), r.GetN(t, "protons"), r.GetN(t, "alpha")]
        if divb_master:
            self.assertIsInstance(r.GetDivB(t), legacy.ScalarField)
            self.assertNotIsInstance(r.GetDivB(t), amr.ScalarField)
        else:
            scalar.append(r.GetDivB(t))
        vector = [r.GetB(t), r.GetE(t), r.GetVi(t), r.GetChargeFlux(t),
                  r.GetFlux(t, "protons"), r.GetJ(t)]
        for f in scalar:
            self.assertIsInstance(f, amr.ScalarField)
        for f in vector:
            self.assertIsInstance(f, amr.VectorField)
        if _registered():
            for f in scalar:
                self.assertIsInstance(f, legacy.ScalarField)
                self.assertIsInstance(f, legacy.PatchHierarchy)
            for f in vector:
                self.assertIsInstance(f, legacy.VectorField)
        for f in scalar + vector:
            self.assertEqual(f.times().tolist(), [t])

    def test_every_getter_on_vtkhdf(self):
        for case in CASES:
            h5, vtk = _dirs(case)
            rv, rl = Run(vtk), Run(h5, legacy=True)
            t = float(rv.times("B")[-1])
            with self.subTest(case=case):
                self._check_classes(rv, t)
                np.testing.assert_array_equal(rv.times("B"), rl.times("B"))
                np.testing.assert_array_equal(rv.times("EM_B"), rl.times("B"))
                np.testing.assert_array_equal(rv.GetDl(), rl.GetDl())
                np.testing.assert_array_equal(rv.GetDl(2), rl.GetDl(2))
                np.testing.assert_array_equal(rv.GetDomainSize(), rl.GetDomainSize())
                for pop in POPS:
                    self.assertEqual(rv.GetMass(pop), rl.GetMass(pop))
                    self.assertEqual(rv.GetCharge(pop), rl.GetCharge(pop))
                vt, lt = rv.all_times(), rl.all_times()
                written = set(lt) - {"tags", "ions_momentum_tensor"} - {
                    f"ions_pop_{p}_momentum_tensor" for p in POPS}
                self.assertEqual(set(vt), written)
                for k in written:
                    np.testing.assert_array_equal(vt[k], np.sort(lt[k]))
                if case.startswith("2d"):
                    Az, (xn, yn) = rv.GetMagneticFlux(t)
                    self.assertEqual(Az.shape, (len(xn), len(yn)))

    def test_every_getter_on_new_h5_path(self):
        for case in CASES:
            h5, _ = _dirs(case)
            r = Run(h5, legacy=False)
            with self.subTest(case=case):
                self._check_classes(r, float(r.times("B")[0]),
                                    divb_master=case.startswith("2d"))
                self.assertIsInstance(r.GetPe(0.0), amr.ScalarField)
                self.assertIsInstance(r.GetJ(0.0, stencil="yee"), amr.VectorField)
                self.assertEqual(r.GetJ(0.0).stencil, "yee")  # h5 default: master's Yee J
                self.assertEqual(r.GetJ(0.0, stencil="node").stencil, "node")

    def test_stored_quantities_agree_between_formats(self):
        getters = [("GetB", ()), ("GetE", ()), ("GetNi", ()), ("GetVi", ()),
                   ("GetMassDensity", ()), ("GetChargeFlux", ()), ("GetN", ("protons",)),
                   ("GetN", ("alpha",)), ("GetFlux", ("protons",)), ("GetFlux", ("alpha",))]
        for case in CASES:
            h5, vtk = _dirs(case)
            rh, rv = Run(h5, legacy=False), Run(vtk)
            for t in rh.times("B"):
                for g, args in getters:
                    a, b = getattr(rh, g)(t, *args), getattr(rv, g)(t, *args)
                    for la, lb in zip(a.levels_data(), b.levels_data(), strict=True):
                        # patch i is the same box in both formats (file order)
                        self.assertEqual(la.geom.key(), lb.geom.key())
                        scale = _level_scale(lb)
                        for n in la.names:
                            d = np.abs(la.buffer(n).astype(np.float64) - lb.buffer(n))
                            self.assertLessEqual(float(d.max()), 2.5e-7 * scale,
                                                 f"{case} t={t} {g}{args} L{la.level} {n}")

    def test_current_and_divergence_agree_at_nodes_with_neighbours(self):
        for case in CASES:
            h5, vtk = _dirs(case)
            rh, rv = Run(h5, legacy=False), Run(vtk)
            for t in rh.times("B"):
                Bscale = max(_level_scale(ld) for ld in rh.GetB(t).levels_data())
                for g in ("GetJ", "GetDivB"):
                    a = (rh.GetJ(t, stencil="node") if g == "GetJ"
                         else rh.GetDivB(t, stencil="node", mask_faces=False))
                    b = getattr(rv, g)(t)
                    for la, lb in zip(a.levels_data(), b.levels_data(), strict=True):
                        inner = ~one_sided_mask(la.geom)
                        if g == "GetJ":
                            bound = 1e-5 * _level_scale(lb)
                        else:
                            bound = 1e-5 * Bscale / float(np.min(la.geom.dl))
                        for n in la.names:
                            d = np.abs(la.buffer(n).astype(np.float64) - lb.buffer(n))[inner]
                            if d.size:
                                self.assertLessEqual(float(d.max()), bound,
                                                     f"{case} t={t} {g} L{la.level} {n}")

    def test_h5_defaults_are_masters_current_and_divergence(self):
        """h5: GetJ() is master's Yee J and GetDivB() master's per-cell divB (same shapes,
        ghosts included, the same values: the same operations in the same order, so exactly
        equal, J rounded to the file's dtype); explicit options still choose the node
        versions."""
        for case in ("2d_static", "2d_tagged"):
            h5, _ = _dirs(case)
            rl, rn = Run(h5, legacy=True), Run(h5, legacy=False)
            for t in rn.times("B"):
                ref, new = rl.GetDivB(t), rn.GetDivB(t)
                ref_levels, new_levels = ref.levels(), new.levels()
                self.assertEqual(ref_levels.keys(), new_levels.keys())
                for il in ref_levels:
                    for pa, pb in zip(ref_levels[il].patches, new_levels[il].patches,
                                      strict=True):
                        a, b = pa.patch_datas["value"].dataset, pb.patch_datas["value"].dataset
                        self.assertEqual(np.shape(a), np.shape(b))
                        np.testing.assert_array_equal(a[:], b[:])  # values are ~roundoff
                J = rn.GetJ(t)
                self.assertEqual(J.stencil, "yee")
                self.assertEqual(_patches_by_box(J).keys(),
                                 _legacy_by_box(rl.GetJ(t), "xyz", t).keys())
                _assert_rounded_master(self, _patches_by_box(J),
                                       _legacy_by_box(rl.GetJ(t), "xyz", t), f"{case} t={t}")
                nodes = rn.GetDivB(t, stencil="node")  # explicit: the node divB
                self.assertIsInstance(nodes, amr.ScalarField)
                self.assertEqual(nodes.stencil, "node")

    def test_yee_current_is_masters_rounded_to_the_file_dtype(self):
        """Same differences, divisions and averages as master, in the same order: master's
        float64 J rounded to the file's dtype, exactly; also on a float64 copy (CI's dtype)."""
        f64 = os.path.join(tmp_dir(self), "2d_static")
        float64_copy(_dirs("2d_static")[0], f64)
        for case, h5 in (("2d_static", _dirs("2d_static")[0]),
                         ("2d_tagged", _dirs("2d_tagged")[0]), ("2d_static float64", f64)):
            rl, rn = Run(h5, legacy=True), Run(h5, legacy=False)
            for t in rn.times("B"):
                _assert_rounded_master(self, _patches_by_box(rn.GetJ(t, stencil="yee")),
                                       _legacy_by_box(rl.GetJ(t), "xyz", t), f"{case} t={t}")

    def test_electron_pressure_is_masters_rounded_to_the_file_dtype(self):
        """Ni * Te in float64, rounded once: master's Pe rounded to the file's dtype, exactly;
        also on a float64 copy (CI's dtype)."""
        f64 = os.path.join(tmp_dir(self), "2d_static")
        float64_copy(_dirs("2d_static")[0], f64)
        cases = [(c, _dirs(c)[0]) for c in ("1d", "2d_static", "2d_tagged")]
        for case, h5 in cases + [("2d_static float64", f64)]:
            rl, rn = Run(h5, legacy=True), Run(h5, legacy=False)
            for t in rn.times("Ni"):
                new = _patches_by_box(rn.GetPe(t))
                ref = _legacy_by_box(rl.GetPe(t), ("value",), f"{t:.10f}")
                self.assertEqual(new.keys(), ref.keys())
                _assert_rounded_master(self, new, ref, f"{case} t={t}")

    def test_new_h5_matches_master_primal_b(self):
        """Every dump; both readers vs the exact node means (_common.assert_matches_master)."""
        for case in ("1d", "2d_static", "2d_tagged"):
            h5, _ = _dirs(case)
            rl, rn = Run(h5, legacy=True), Run(h5, legacy=False)
            for t in rn.times("B"):
                self.assertGreater(check_against_master(rl, rn, "GetB", (), t, case), 0)

    def test_multi_time_keeps_callers_order(self):
        _, vtk = _dirs("2d_tagged")
        r = Run(vtk)
        both = r.GetB([0.1, 0.0])
        np.testing.assert_array_equal(both.times(), [0.1, 0.0])
        for t in (0.0, 0.1):
            one = r.GetB(t)
            for la, lb in zip(both.levels_data(t), one.levels_data()):
                np.testing.assert_array_equal(la.buffer("y"), lb.buffer("y"))
        J = r.GetJ([0.1, 0.0])
        np.testing.assert_array_equal(J.times(), [0.1, 0.0])
        with self.assertRaises(ValueError):
            J.mean()  # needs one time
        self.assertEqual(J.at(0.0).levels_data()[0].geom.key(),
                         r.GetJ(0.0).levels_data()[0].geom.key())

    def test_divb_mask_faces(self):
        h5, vtk = _dirs("2d_tagged")
        for r in (Run(vtk), Run(h5, legacy=False)):
            d = r.GetDivB(0.1, stencil="node", mask_faces=True)
            for ld in d.levels_data():
                np.testing.assert_array_equal(np.isnan(ld.buffer("value")),
                                              one_sided_mask(ld.geom))

    def test_divb_stencil_is_an_explicit_choice(self):
        """GetDivB's result depends on stencil= only: "yee" = master's per-cell object (2D
        h5 only, the 2D h5 default), "node" = the new reader's node divB (default elsewhere);
        box=/mask_faces= never switch from one to the other silently."""
        box = dict(x=(5.0, 15.0), y=(5.0, 15.0))
        h5, vtk = _dirs("2d_tagged")
        rh, rl = Run(h5, legacy=False), Run(h5, legacy=True)
        t = 0.1
        for d in (rh.GetDivB(t), rh.GetDivB(t, stencil="yee"),
                  rh.GetDivB(t, stencil="yee", selection_box=(0, 0, 10, 10))):
            self.assertIsInstance(d, legacy.ScalarField)
            self.assertNotIsInstance(d, amr.ScalarField)
        ref = rl.GetDivB(t).level(0).patches[0].patch_datas["value"].dataset[:]
        np.testing.assert_array_equal(
            rh.GetDivB(t, stencil="yee").level(0).patches[0].patch_datas["value"].dataset[:],
            ref)
        for kw in (dict(), dict(stencil="yee")):
            with self.assertRaisesRegex(TypeError, "selection_box.*stencil='node'"):
                rh.GetDivB(t, box=box, **kw)
            with self.assertRaisesRegex(TypeError, "stencil='node'"):
                rh.GetDivB(t, mask_faces=False, **kw)
        for kw in (dict(), dict(box=box), dict(mask_faces=False)):
            self.assertEqual(rh.GetDivB(t, stencil="node", **kw).stencil, "node")
        with self.assertRaisesRegex(ValueError, "'node' or 'yee'"):
            rh.GetDivB(t, stencil="cell")
        # vtkhdf (any dimension) and h5 in 1D/3D: node divB only
        for case in CASES:
            h5c, vtkc = _dirs(case)
            runs = [Run(vtkc)] + ([] if case.startswith("2d") else [Run(h5c, legacy=False)])
            for r in runs:
                with self.subTest(case=case, run=repr(r)):
                    self.assertEqual(r.GetDivB(0.0).stencil, "node")
                    self.assertEqual(r.GetDivB(0.0, stencil="node").stencil, "node")
                    with self.assertRaisesRegex(ValueError, "stencil='node'"):
                        r.GetDivB(0.0, stencil="yee")

    def test_all_primal_false_on_new_h5_is_masters_yee_object(self):
        h5, _ = _dirs("2d_tagged")
        r = Run(h5, legacy=False)
        raw = r.GetB(0.1, all_primal=False)
        self.assertEqual(type(raw).__module__.split(".")[2], "hierarchy")
        self.assertIn("Bx", raw.quantities())
        with self.assertRaisesRegex(TypeError, "legacy reader"):
            raw.Bx + r.GetB(0.1).x
        Jraw = r.GetJ(0.1, all_primal=False)
        self.assertIn("Jx", Jraw.quantities())


# ============================================================================ errors


class RaisesTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require()

    def setUp(self):
        super().setUp()
        h5, vtk = _dirs("2d_tagged")
        self.vtk, self.h5 = Run(vtk), Run(h5, legacy=False)

    def test_merged_points_to_finest_and_probes(self):
        for getter in ("GetB", "GetNi", "GetJ", "GetDivB"):
            with self.assertRaisesRegex(ValueError, "finest"):
                getattr(self.vtk, getter)(0.1, merged=True)

    def test_merged_on_h5_is_masters_code(self):
        """merged=True on h5: master's result (or master's exception), not an error of ours."""
        h5, _ = _dirs("2d_tagged")
        master = LegacyRun(h5)
        for getter in ("GetB", "GetE", "GetNi", "GetVi", "GetJ", "GetDivB", "GetPe"):
            with self.subTest(getter):
                outcome = []
                for r in (self.h5, master):
                    try:
                        x = getattr(r, getter)(0.1, merged=True)
                        outcome.append((type(x), sorted(x) if isinstance(x, dict) else None))
                    except (AttributeError, KeyError) as e:  # master's own failures (dict)
                        outcome.append((type(e), str(e)))
                self.assertEqual(outcome[0], outcome[1])

    def test_unknown_keyword_named(self):
        for r in (self.vtk, self.h5):
            with self.assertRaisesRegex(TypeError, "colour.*accepts: time"):
                r.GetB(0.1, colour="red")
            with self.assertRaisesRegex(TypeError, "stencl"):
                r.GetJ(0.1, stencl="yee")

    def test_interp_only_default(self):
        self.vtk.GetB(0.1, interp="nearest")
        with self.assertRaisesRegex(TypeError, "merged=True"):
            self.vtk.GetB(0.1, interp="linear")

    def test_all_primal_false_on_vtkhdf(self):
        with self.assertRaisesRegex(ValueError, "no Yee-centred data"):
            self.vtk.GetB(0.1, all_primal=False)

    def test_yee_current_on_vtkhdf(self):
        with self.assertRaisesRegex(ValueError, "stencil='node'"):
            self.vtk.GetJ(0.1, stencil="yee")

    def test_pe_on_vtkhdf(self):
        with self.assertRaisesRegex(ValueError, "no simulation parameters"):
            self.vtk.GetPe(0.1)

    def test_unavailable_getters(self):
        calls = {
            "GetPressure": lambda r: r.GetPressure(0.1, "protons"),
            "GetPi": lambda r: r.GetPi(0.1),
            "GetTags": lambda r: r.GetTags(0.1),
            "GetRanks": lambda r: r.GetRanks(0.1),
            "GetParticles": lambda r: r.GetParticles(0.1, "protons"),
            "GetParticleCount": lambda r: r.GetParticleCount(0.1),
        }
        for name, call in calls.items():
            with self.subTest(name), self.assertRaisesRegex(ValueError,
                                                            "vtkhdf writer does not write"):
                call(self.vtk)
        # h5: master's reader serves them, with master's inputs (not the new reader's
        # node-centred fields): the same hierarchy as Run(path, legacy=True)
        h5, _ = _dirs("2d_tagged")
        master = Run(h5, legacy=True)
        for name in ("GetPressure", "GetPi", "GetTags", "GetRanks"):
            with self.subTest(name):
                got, ref = calls[name](self.h5), calls[name](master)
                self.assertIs(type(got), type(ref))
                self.assertEqual(_hier_arrays(got, 0.1).keys(), _hier_arrays(ref, 0.1).keys())
                for key, arrays in _hier_arrays(got, 0.1).items():
                    for qty, a in arrays.items():
                        np.testing.assert_array_equal(a, _hier_arrays(ref, 0.1)[key][qty],
                                                      err_msg=f"{name} {key} {qty}")
        for name in ("GetParticles", "GetParticleCount"):  # not written by these runs
            with self.subTest(name), self.assertRaises(FileNotFoundError):
                calls[name](self.h5)

    def test_time_not_found_lists_nearest(self):
        with self.assertRaisesRegex(TimeNotFoundError, "nearest"):
            self.vtk.GetB(0.123)
        with self.assertRaises(KeyError):
            self.h5.GetNi(0.123)

    def test_missing_time_is_caught_as_masters_time_not_found_error(self):
        """A post_advance hook catches pyphare.pharesee.hierarchy.TimeNotFoundError (master's
        Simulator docs): every reader raises that class for a time without a dump."""
        h5, _ = _dirs("2d_tagged")
        runs = {"vtkhdf": self.vtk, "h5": self.h5, "h5 legacy": Run(h5, legacy=True)}
        for name, r in runs.items():
            for getter in ("GetB", "GetNi"):
                with self.subTest(name, getter=getter), \
                        self.assertRaises(legacy.TimeNotFoundError):
                    getattr(r, getter)(0.123)


# ============================================================================ configuration


class ConfigurationTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require()

    """The legacy path never reads the new reader's configuration (A8)."""

    def setUp(self):
        super().setUp()
        self._saved = {k: os.environ.get(k) for k in ("PHARESEE_GEOM_CACHE_MB",)}

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        os.environ.pop("PHARESEE_LEGACY", None)
        super().tearDown()

    def test_bad_new_reader_setting_does_not_break_the_legacy_path(self):
        h5, vtk = _dirs("2d_tagged")
        os.environ["PHARESEE_GEOM_CACHE_MB"] = "xyz"
        r = Run(h5, legacy=True)
        self.assertIs(type(r.GetB(0.1)), legacy.VectorField)
        self.assertIsNone(r._geom_cache)
        # the new reader does read it, and says which variable is wrong
        with self.assertRaisesRegex(ValueError, "PHARESEE_GEOM_CACHE_MB"):
            Run(vtk).GetB(0.1)

    def test_legacy_argument_is_a_bool_or_none(self):
        h5, _ = _dirs("2d_tagged")
        for bad in ("0", "False", 0, 1, np.bool_(False)):
            with self.subTest(bad=bad), self.assertRaisesRegex(TypeError, "True, False or None"):
                Run(h5, legacy=bad)
        self.assertEqual(Run(h5, legacy=False).reader, "new")
        self.assertEqual(Run(h5, legacy=True).reader, "legacy")

    def test_environment_variable_parsed_strictly(self):
        h5, _ = _dirs("2d_tagged")
        cases = {"0": "new", "no": "new", "False": "new", " NO ": "new", "1": "legacy",
                 "yes": "legacy", "TRUE": "legacy", "": "new"}
        for value, reader in cases.items():
            os.environ["PHARESEE_LEGACY"] = value
            with self.subTest(value=value):
                self.assertEqual(Run(h5).reader, reader)
        for value in ("2", "abc", "off"):
            os.environ["PHARESEE_LEGACY"] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError,
                                                                   "PHARESEE_LEGACY"):
                Run(h5)

    def test_all_times_skips_a_file_that_cannot_be_opened(self):
        _, vtk = _dirs("2d_tagged")
        with tempfile.TemporaryDirectory() as tmp:
            for name in ("EM_B.vtkhdf", "ions_charge_density.vtkhdf"):
                os.symlink(os.path.join(vtk, name), os.path.join(tmp, name))
            bad = os.path.join(tmp, "EM_E.vtkhdf")
            with open(os.path.join(vtk, "EM_E.vtkhdf"), "rb") as f:
                head = f.read(4096)
            with open(bad, "wb") as f:
                f.write(head)  # a truncated file
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                times = Run(tmp).all_times()
            self.assertEqual(set(times), {"B", "Ni"})
            msgs = [str(x.message) for x in w if issubclass(x.category, RuntimeWarning)]
            self.assertEqual(len(msgs), 1, msgs)
            self.assertIn(bad, msgs[0])
            self.assertIn("Error", msgs[0])  # the error's type and text
            np.testing.assert_array_equal(times["B"], Run(vtk).times("B"))

    def test_all_times_keeps_h5_files_the_new_reader_cannot_parse(self):
        """tags.h5 (no dataset) and a particle-like file are listed with master's times."""
        import h5py

        h5, _ = _dirs("2d_tagged")
        with tempfile.TemporaryDirectory() as tmp:
            for name in ("EM_B.h5", "tags.h5"):
                os.symlink(os.path.join(h5, name), os.path.join(tmp, name))
            with h5py.File(os.path.join(tmp, "ions_pop_protons_domain.h5"), "w") as f:
                f.create_group("t/0.0000000000/pl0/p0#0")
                f.create_group("t/0.1000000000/pl0/p0#0")
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                times = Run(tmp).all_times()
            master = LegacyRun(tmp).all_times()
            self.assertEqual(set(times), set(master))
            self.assertEqual(set(times), {"B", "tags", "ions_pop_protons_domain"})
            for key, t in master.items():
                np.testing.assert_array_equal(times[key], np.sort(t))

    def test_empty_time_list_raises_clearly(self):
        h5, vtk = _dirs("2d_tagged")
        for r in (Run(vtk), Run(h5, legacy=False)):
            for getter, empty in (("GetB", []), ("GetJ", []), ("GetNi", ()),
                                  ("GetB", np.array([]))):
                with self.assertRaisesRegex(ValueError, "no time given"):
                    getattr(r, getter)(empty)


# ============================================================================ reuse (F3)


class _CountingRun:
    """Counts every backend read of a Run (read_into and read_derived), all files."""

    def __init__(self, run):
        self.run = run
        self.calls = 0
        orig = run._backend

        def backend(path):
            b = orig(path)
            if not getattr(b, "_counted", False):
                for name in ("read_into", "read_derived"):
                    f = getattr(b, name, None)
                    if f is not None:
                        setattr(b, name, self._wrap(f))
                b._counted = True
            return b

        run._backend = backend

    def _wrap(self, f):
        def wrapped(*args, **kw):
            self.calls += 1
            return f(*args, **kw)

        return wrapped


class ReuseTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require()

    BOX = PhysBox((8.0, 6.0), (22.0, 18.0)).to_dict()

    def _runs(self):
        h5, vtk = _dirs("2d_tagged")
        return (Run(vtk), Run(h5, legacy=False))

    def test_second_getter_call_shares_the_data(self):
        import gc
        import weakref

        for r in self._runs():
            with self.subTest(r=repr(r)):
                c = _CountingRun(r)
                B = r.GetB(0.1).load()
                self.assertGreater(c.calls, 0)
                c.calls = 0
                B2 = r.GetB(0.1)
                self.assertIs(B2.levels_data()[0], B.levels_data()[0])
                B2.load()
                B2.x.mean()
                self.assertEqual(c.calls, 0, "second GetB(t) reads nothing")
                self.assertIsNot(B2, B)  # a new object sharing the data
                # derived quantities built on the held B read nothing
                J = r.GetJ(0.1, stencil="node")  # (h5's default Yee J reads the file's Yee B)
                J.x.max()
                self.assertEqual(c.calls, 0, "GetJ after B.load() reads nothing")
                ref = Run(r.path, legacy=False).GetJ(0.1, stencil="node")
                for la, lb in zip(J.levels_data(), ref.levels_data(), strict=True):
                    np.testing.assert_array_equal(la.buffer("z"), lb.buffer("z"))
                # a boxed request is cut from the loaded whole-domain field in memory
                part = r.GetB(0.1, box=self.BOX)
                fresh = Run(r.path, legacy=False).GetB(0.1, box=self.BOX)
                self.assertEqual(part.extent, fresh.extent)
                for la, lb in zip(part.levels_data(), fresh.levels_data(), strict=True):
                    self.assertEqual(la.geom.key(), lb.geom.key())
                    np.testing.assert_array_equal(la.buffer("y"), lb.buffer("y"))
                Jb = r.GetJ(0.1, box=self.BOX, stencil="node")
                Jb.load()
                self.assertEqual(c.calls, 0, "boxed reads come from memory")
                # no strong reference: dropping the fields frees them, then reads again
                wr = weakref.ref(B)
                del B, B2, J, part, Jb
                gc.collect()
                self.assertIsNone(wr())
                r.GetB(0.1).load()
                self.assertGreater(c.calls, 0)

    def test_derived_from_held_density_and_flux(self):
        _, rh = self._runs()
        c = _CountingRun(rh)
        Ni = rh.GetNi(0.1).load()
        c.calls = 0
        Pe = rh.GetPe(0.1)
        Pe.mean()
        self.assertEqual(c.calls, 0)
        rv = self._runs()[0]
        cv = _CountingRun(rv)
        B = rv.GetB(0.1).load()
        cv.calls = 0
        rv.GetMagneticFlux(0.1)
        self.assertEqual(cv.calls, 0)
        del Ni, B

    def test_boxed_request_reads_file_when_whole_is_not_loaded(self):
        r = self._runs()[0]
        c = _CountingRun(r)
        B = r.GetB(0.1)  # alive but nothing loaded
        part = r.GetB(0.1, box=self.BOX)
        part.load()
        self.assertGreater(c.calls, 0)
        self.assertLess(part.nbytes(), B.nbytes())
        self.assertEqual(B.loaded_nbytes, 0)  # the boxed read did not fill B
        # a whole field held only weakly by derived results is not kept alive
        part2 = r.GetB(0.1, box=self.BOX)
        self.assertIsNot(part2.levels_data()[0], part.levels_data()[0])


class DivergenceFormatsTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require(cases=CASES)

    def test_vtkhdf_masks_one_sided_nodes_by_default_and_agrees_with_h5(self):
        for case in CASES:
            h5, vtk = _dirs(case)
            rh, rv = Run(h5, legacy=False), Run(vtk)
            for t in rh.times("B"):
                with self.subTest(case=case, t=t):
                    dv, dh = rv.GetDivB(t), rh.GetDivB(t, stencil="node", mask_faces=False)
                    Bscale = max(_level_scale(ld) for ld in rh.GetB(t).levels_data())
                    worst_v = worst_h = 0.0
                    for lv, lh in zip(dv.levels_data(), dh.levels_data(), strict=True):
                        mask = one_sided_mask(lv.geom)
                        v = lv.buffer("value").astype(np.float64)
                        h = lh.buffer("value").astype(np.float64)
                        np.testing.assert_array_equal(np.isnan(v), mask)
                        self.assertFalse(np.isnan(h).any())  # h5: values kept
                        if (~mask).any():
                            worst_v = max(worst_v, float(np.abs(v[~mask]).max()))
                            worst_h = max(worst_h, float(np.abs(h[~mask]).max()))
                    dl = float(np.min(rh.GetDl()))
                    self.assertLessEqual(abs(worst_v - worst_h), 1e-5 * Bscale / dl)
                    # mask_faces=False gives values everywhere on vtkhdf too
                    full = rv.GetDivB(t, mask_faces=False)
                    self.assertFalse(any(np.isnan(ld.buffer("value")).any()
                                         for ld in full.levels_data()))


# ============================================================================ flux function


class MagneticFluxTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require(cases=("2d_static", "2d_tagged"))

    def test_flux_new_vs_master(self):
        for case in ("2d_static", "2d_tagged"):
            h5, vtk = _dirs(case)
            rl, rn, rv = Run(h5, legacy=True), Run(h5, legacy=False), Run(vtk)
            t = float(rn.times("B")[-1])
            Al, (xl, yl) = rl.GetMagneticFlux(t)
            An, (xn, yn) = rn.GetMagneticFlux(t)
            Av, _ = rv.GetMagneticFlux(t)
            # master's grid: the finest level's nodes
            np.testing.assert_allclose(xl, xn)
            np.testing.assert_allclose(yl, yn)
            self.assertEqual(Al.shape, An.shape)
            scale = np.abs(Al).max()
            print(case, "flux max|new-master|/max|Az| =", np.abs(Al - An).max() / scale)
            self.assertLess(np.abs(Al - An).max(), 0.05 * scale, case)
            self.assertLess(np.abs(Av - An).max(), 1e-6 * scale, case)

    def test_flux_interp(self):
        """interp samples given xn/yn grids ("nearest" or "linear"/"bilinear", others raise
        as on master); on the default grid (the finest nodes) there is nothing to sample and
        the result is the same for every interp (measured against master's "bilinear":
        0.6 % of max|Az| on 2d_static, 2e-8 on 2d_tagged; master's Delaunay interpolator
        takes ~15 s here, so it is not run)."""
        h5, _ = _dirs("2d_static")
        rn = Run(h5, legacy=False)
        t = float(rn.times("B")[-1])
        A0, (xn, yn) = rn.GetMagneticFlux(t)
        for interp in ("nearest", "linear", "bilinear"):
            np.testing.assert_array_equal(rn.GetMagneticFlux(t, interp=interp)[0], A0)
        xs, ys = xn[1:-1:7] + 0.3 * (xn[1] - xn[0]), yn[1:-1:5]
        near, lin = (rn.GetMagneticFlux(t, interp=i, xn=xs, yn=ys)[0]
                     for i in ("nearest", "bilinear"))
        self.assertEqual(near.shape, (len(xs), len(ys)))
        self.assertGreater(np.abs(near - lin).max(), 0.0)
        for kw in (dict(), dict(xn=xs, yn=ys)):
            with self.assertRaisesRegex(ValueError, "'nearest' or 'bilinear'"):
                rn.GetMagneticFlux(t, interp="cubic", **kw)

    def test_reconnection_rate_uses_masters_algorithm(self):
        _, vtk = _dirs("2d_tagged")
        r = Run(vtk)
        times = r.times("B")
        tc, rates, flux, traj = r.GetReconnectionRate(times)
        np.testing.assert_allclose(tc, (times[:-1] + times[1:]) / 2)
        self.assertEqual(traj.shape, (len(times), 2))
        for i, t in enumerate(times):
            Az, (xn, yn) = r.GetMagneticFlux(t)
            x, y, idx = r.FindPrimaryXPoint(Az, xn, yn)
            self.assertEqual(flux[i], Az[idx])
            self.assertEqual(list(traj[i]), [x, y])
        np.testing.assert_allclose(rates, np.diff(flux) / np.diff(times))


def _copy_h5(src, dst, ntimes=None):
    """Copy of the phareh5 file ``src`` at ``dst`` with only its first ``ntimes`` times."""
    import h5py

    with h5py.File(src, "r") as fs, h5py.File(dst, "w") as fd:
        fd.attrs.update(fs.attrs)
        for key in fs:
            if key != "t":
                fs.copy(fs[key], fd, name=key)
        fd.create_group("t")
        for k in list(fs["t"])[:ntimes]:
            fs.copy(fs["t"][k], fd["t"], name=k)


def _append_h5_time(src, dst, index):
    """Append time number ``index`` of ``src`` to the phareh5 file ``dst`` (a run writing)."""
    import h5py

    with h5py.File(src, "r") as fs, h5py.File(dst, "a") as fd:
        k = list(fs["t"])[index]
        fs.copy(fs["t"][k], fd["t"], name=k)


class TestRoutingEdges(LegacyEnv, unittest.TestCase):
    """Calls at the seam between the new reader and master's (review threads of #1595)."""

    def setUp(self):
        super().setUp()
        require(cases=("2d_static",))
        self.h5, self.vtk = _dirs("2d_static")

    def test_box_on_masters_path_raises(self):
        # master's reader has no box=: it would read the whole domain without saying so
        r = Run(self.h5)
        t = float(r.times("B")[0])
        box = dict(x=(0.0, 10.0), y=(0.0, 10.0))
        with self.assertRaisesRegex(TypeError, "selection_box"):
            r.GetB(t, merged=True, box=box)
        with self.assertRaisesRegex(TypeError, "selection_box"):
            Run(self.h5, legacy=True).GetB(t, box=box)

    def test_arguments_handed_to_master_are_bound_by_name(self):
        # zip() alone would drop an extra value without a word
        from pyphare.pharesee.run.run import _master_arguments

        self.assertEqual(_master_arguments("GetTags", (0.5, True)),
                         {"time": 0.5, "merged": True})
        with self.assertRaisesRegex(TypeError, r"GetTags\(\) takes 2 positional"):
            _master_arguments("GetTags", (0.5, True, "nearest"))
        r = Run(self.h5)  # merged=True is served by master: the same result as master's call
        t = float(r.times("B")[0])
        got, ref = r.GetB(t, merged=True), Run(self.h5, legacy=True).GetB(t, merged=True)
        self.assertEqual(got.keys(), ref.keys())
        for k in ref:
            np.testing.assert_array_equal(got[k][0].values, ref[k][0].values)
            assert_all_equal(got[k][1], ref[k][1], k)

    def test_all_times_with_format_skips_stems_of_the_other_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            os.symlink(os.path.join(self.h5, "EM_B.h5"), os.path.join(tmp, "EM_B.h5"))
            os.symlink(os.path.join(self.vtk, "EM_B.vtkhdf"), os.path.join(tmp, "EM_B.vtkhdf"))
            os.symlink(os.path.join(self.h5, "ions_charge_density.h5"),
                       os.path.join(tmp, "ions_charge_density.h5"))
            self.assertEqual(set(Run(tmp, format="vtkhdf").all_times()), {"B"})
            self.assertEqual(set(Run(tmp, format="h5").all_times()), {"B", "Ni"})

    def test_dumps_appended_after_opening_are_visible(self):
        src = os.path.join(self.h5, "EM_B.h5")
        with tempfile.TemporaryDirectory() as tmp:
            dst = os.path.join(tmp, "EM_B.h5")
            _copy_h5(src, dst, ntimes=1)
            r = Run(tmp)
            self.assertEqual(len(r.times("B")), 1)
            _append_h5_time(src, dst, 1)
            times = r.times("B")
            self.assertEqual(len(times), 2)
            ref = Run(self.h5).GetB(float(times[1]))
            np.testing.assert_array_equal(r.GetB(float(times[1])).x.finest().data,
                                          ref.x.finest().data)

    def test_a_held_field_is_reused_only_while_its_file_is_unchanged(self):
        """A loaded field this Run returned, still held by the caller, is reused at the same
        time after an append; after the file is rewritten in place or replaced, the getter
        reads the new file."""
        import h5py

        src = os.path.join(self.h5, "EM_B.h5")
        for how in ("rewritten", "replaced"):
            with self.subTest(how=how), tempfile.TemporaryDirectory() as tmp:
                dst = os.path.join(tmp, "EM_B.h5")
                _copy_h5(src, dst, ntimes=1)
                r = Run(tmp)
                t = float(r.times("B")[0])
                held = r.GetB(t)
                held.load()
                old = held.x.finest().data
                _append_h5_time(src, dst, 1)
                self.assertIs(r.GetB(t).snapshot(t), held.snapshot(t))  # append: reused
                new = os.path.join(tmp, "new.h5") if how == "replaced" else dst
                if how == "replaced":
                    _copy_h5(dst, new)
                with h5py.File(new, "r+") as f:
                    f["t"].visititems(lambda _, o: o.write_direct(2 * o[()])
                                      if isinstance(o, h5py.Dataset) else None)
                if how == "replaced":
                    os.replace(new, dst)
                np.testing.assert_array_equal(r.GetB(t).x.finest().data, 2 * old)

    def test_mass_of_a_particle_only_output(self):
        import h5py

        with tempfile.TemporaryDirectory() as tmp:
            with h5py.File(os.path.join(tmp, "ions_pop_protons_domain.h5"), "w") as f:
                f.attrs["pop_mass"] = 2.0
                f.attrs["pop_charge"] = 1.0
                f.create_group("t")
            r = Run(tmp)
            self.assertEqual(r.GetMass("protons"), 2.0)
            self.assertEqual(r.GetCharge("protons"), 1.0)

    def test_plan_of_J_on_h5_follows_GetJ(self):
        # GetJ on h5 reads the Yee B through its own opening: a held, loaded GetB does not
        # make it free
        r = Run(self.h5)
        t = float(r.times("B")[0])
        B = r.GetB(t)
        B.load()
        self.assertEqual(r.plan("B", t).nbytes, 0)
        self.assertGreater(r.plan("J", t).nbytes, 0)

    def test_pickle_keeps_the_simulation(self):
        """A pickled h5 field keeps its (lazy) simulation, before and after it was loaded;
        the loaded Simulation (local user functions) is never pickled."""
        import pickle

        t = float(Run(self.h5).times("B")[1])
        for load_first in (False, True):
            B = Run(self.h5).GetB(t)
            dl = tuple(Run(self.h5).GetB(t).sim.dl)
            if load_first:
                self.assertEqual(tuple(B.sim.dl), dl)
            for f in (B, B.x, B * 2):
                with self.subTest(load_first=load_first, f=f.quantity):
                    back = pickle.loads(pickle.dumps(f))
                    self.assertIsNotNone(back.sim)
                    self.assertEqual(tuple(back.sim.dl), dl)


if __name__ == "__main__":
    unittest.main()
