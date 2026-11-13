"""Differential table: every Run getter on phareh5 files, legacy reader vs new reader
(DESIGN 11.4c). The new reader is the default for h5; this table is what changed for
users of master's reader (``Run(path, legacy=True)`` or ``PHARESEE_LEGACY=1`` restores it).

For each getter, on ``build/sample_data/2d_static/h5`` at t = 0.007 (3 levels), the test
records for the legacy (``Run(path, legacy=True)``) and the new (``Run(path, legacy=False)``) result:
  class        module tail + class name, or "raises <Exception>"
  times        dtype kind of ``result.times()`` ("U" = strings on master, "f" = floats)
  dtype        dtype of the first patch's ``dataset`` (per-patch views on the new path)
  shape        shape of that ``dataset`` (level 0, patch 0)
  values       max |new - legacy| / max |legacy| over primal nodes, patches paired by box
               (GetDivB: legacy cell values averaged to the nodes around which they sit;
               GetMagneticFlux: on the common level-0 nodes)
It prints the table, then asserts it equals ``APPROVED`` (values: below the approved bound).
A change in either reader shows up here as a failing row to review, not a silent drift.
"""

import os
import unittest
import warnings

import numpy as np

from pyphare.pharesee.run import Run
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES

os.environ.setdefault("PHARESEE_CACHE", "0")

H5 = os.path.join(SAMPLES, "2d_static", "h5")
T = 0.007

# getter -> (call, components compared)
GETTERS = {
    "GetB": (lambda r: r.GetB(T), "xyz"),
    "GetE": (lambda r: r.GetE(T), "xyz"),
    "GetNi": (lambda r: r.GetNi(T), ("value",)),
    "GetMassDensity": (lambda r: r.GetMassDensity(T), ("value",)),
    "GetN": (lambda r: r.GetN(T, "protons"), ("value",)),
    "GetVi": (lambda r: r.GetVi(T), "xyz"),
    "GetChargeFlux": (lambda r: r.GetChargeFlux(T), "xyz"),
    "GetFlux": (lambda r: r.GetFlux(T, "alpha"), "xyz"),
    "GetJ": (lambda r: r.GetJ(T), "xyz"),
    "GetJ(stencil=yee)": (lambda r: r.GetJ(T, stencil="yee") if r.reader == "new"
                          else r.GetJ(T), "xyz"),
    "GetDivB": (lambda r: r.GetDivB(T), ("value",)),
    "GetPe": (lambda r: r.GetPe(T), ("value",)),
    "GetB(all_primal=False)": (lambda r: r.GetB(T, all_primal=False), ("Bx", "By", "Bz")),
    "GetPressure": (lambda r: r.GetPressure(T, "protons"), None),
    "GetPi": (lambda r: r.GetPi(T), None),
    "GetTags": (lambda r: r.GetTags(T), None),
    "GetRanks": (lambda r: r.GetRanks(T), None),
    "GetParticles": (lambda r: r.GetParticles(T, "protons"), None),
    "GetParticleCount": (lambda r: r.GetParticleCount(T), None),
    "GetB(merged=True)": (lambda r: r.GetB(T, merged=True), None),
}

# the approved differences; values = upper bound of the relative difference
APPROVED = {
    "GetB": ("hierarchy.VectorField", "amr.VectorField", "U", "f", "float64", "float32",
             (105, 85), (105, 85), 1.2e-7),
    "GetE": ("hierarchy.VectorField", "amr.VectorField", "U", "f", "float64", "float32",
             (105, 85), (105, 85), 1.2e-7),
    "GetNi": ("hierarchy.ScalarField", "amr.ScalarField", "U", "f", "float32", "float32",
              (105, 85), (105, 85), 0.0),
    "GetMassDensity": ("hierarchy.ScalarField", "amr.ScalarField", "U", "f", "float32",
                       "float32", (105, 85), (105, 85), 0.0),
    "GetN": ("hierarchy.ScalarField", "amr.ScalarField", "U", "f", "float32", "float32",
             (105, 85), (105, 85), 0.0),
    "GetVi": ("hierarchy.VectorField", "amr.VectorField", "U", "f", "float32", "float32",
              (105, 85), (105, 85), 0.0),
    "GetChargeFlux": ("hierarchy.VectorField", "amr.VectorField", "U", "f", "float32",
                      "float32", (105, 85), (105, 85), 0.0),
    "GetFlux": ("hierarchy.VectorField", "amr.VectorField", "U", "f", "float32", "float32",
                (105, 85), (105, 85), 0.0),
    # node stencil vs master's Yee stencil: a (1,2,1)/4 smoothing, tens of % on noisy data
    "GetJ": ("hierarchy.VectorField", "amr.VectorField", "U", "f", "float64", "float32",
             (105, 85), (105, 85), 1.0),
    "GetJ(stencil=yee)": ("hierarchy.VectorField", "amr.VectorField", "U", "f", "float64",
                          "float32", (105, 85), (105, 85), 1.2e-7),
    # master: cells (dual/dual, no NaN ring); new: nodes, the cell values averaged
    "GetDivB": ("hierarchy.ScalarField", "amr.ScalarField", "U", "f", "float64", "float32",
                (104, 84), (105, 85), 1e-6),
    "GetPe": ("hierarchy.ScalarField", "amr.ScalarField", "U", "f", "float64", "float32",
              (105, 85), (105, 85), 1.2e-7),
    # all_primal=False on the new h5 path returns master's object
    "GetB(all_primal=False)": ("hierarchy.PatchHierarchy", "hierarchy.PatchHierarchy", "U",
                               "U", "float32", "float32", (105, 84), (105, 84), 0.0),
    # not on the new reader: master's reader serves them on the new h5 path too, with
    # master's inputs (test_run_amr.RaisesTest compares the values)
    "GetPressure": ("hierarchy.PatchHierarchy", "hierarchy.PatchHierarchy", "U", "U",
                    "float64", "float64", (105, 85), (105, 85), None),
    "GetPi": ("hierarchy.PatchHierarchy", "hierarchy.PatchHierarchy", "U", "U", "float32",
              "float32", (105, 85), (105, 85), None),
    # master's tags hierarchy has no dataset on its patches in this run
    "GetTags": ("hierarchy.PatchHierarchy", "hierarchy.PatchHierarchy", "U", "U", None,
                None, None, None, None),
    "GetRanks": ("hierarchy.ScalarField", "hierarchy.ScalarField", "U", "U", "float64",
                 "float64", (104, 84), (104, 84), None),
    "GetParticles": ("raises FileNotFoundError", "raises FileNotFoundError", None, None, None, None,
                     None, None, None),
    "GetParticleCount": ("raises FileNotFoundError", "raises FileNotFoundError", None, None, None, None,
                         None, None, None),
    # merged=True on h5: master's interpolators on the new path too
    "GetB(merged=True)": ("dict", "dict", None, None, None, None, None, None, None),
}
# non-field getters: name -> (legacy summary, new summary)
APPROVED_META = {
    "times('B')": ("float64 [0.0, 0.007, 0.014, 0.02]", "float64 [0.0, 0.007, 0.014, 0.02]"),
    "times('EM_B')": ("raises KeyError", "float64 [0.0, 0.007, 0.014, 0.02]"),
    # tags.h5 of this run holds no dataset: the new reader leaves it out of all_times()
    "all_times() keys": ("14 keys", "13 keys"),
    "GetDl(1)": ("[0.2, 0.25]", "[0.2, 0.25]"),
    "GetDomainSize()": ("[40.0, 40.0]", "[40.0, 40.0]"),
    "GetMass('alpha')": ("4.0", "4.0"),
    "GetCharge('alpha')": ("2.0", "2.0"),
    "GetMagneticFlux": ("(401, 321) float64", "(101, 81) float64"),
    "unknown keyword GetB(foo=1)": ("accepted", "raises TypeError"),
    "time not dumped GetB(0.5)": ("raises KeyError", "raises TimeNotFoundError"),
}
# behaviours of master's per-patch API: name -> (legacy, new)
APPROVED_BEHAVIOUR = {
    # stored quantities (moments) had real ghost values on master; the new views show the
    # primal nodes inside a NaN ring for every quantity
    "ghost nodes of stored quantities (Ni)": ("real values", "NaN"),
    "ghost nodes of B (all_primal)": ("NaN", "NaN"),
    # FieldDataView.compare treats NaN == NaN at the same place
    "hierarchy_compare(Ni, Ni again)": ("True", "True"),
    # master's all_primal B has a NaN ring too, and its compare says NaN != NaN
    "hierarchy_compare(B, B again)": ("False", "True"),
    "patch.id level 0 patch 0": ("p0#0", "p0#0"),
    "B.x patch_datas keys": ("['x']", "['x']"),
    "B.x quantities()": ("['x']", "['x']"),
    "B.x.global_max('x') == B.global_max('x')": ("True", "True"),
    "B.x['x']": ("raises KeyError", "raises KeyError"),
    # master: identity (False, silently); new: TypeError pointing to np.equal/np.isclose
    "Ni == 1.0": ("False", "raises TypeError"),
    "Ni == Ni": ("True", "True"),
    # selection_box: master's patches (integer Box rule, edge-touching kept) on both
    "selection_box=(5, 5, 20, 13.5) patches": ("3", "3"),
    "selection_box=(0, 0, 10, 10) patches": ("3", "3"),
    "GetB([])": ("raises IndexError", "raises ValueError"),
}
MAGNETIC_FLUX_BOUND = 0.05  # of max |Az| (interpolation, see test_run_amr.py)


def _class_name(obj):
    if isinstance(obj, dict):
        return "dict"
    mod = type(obj).__module__
    tag = "amr" if mod.startswith("pyphare.pharesee.amr") else mod.split(".")[2]
    return f"{tag}.{type(obj).__name__}"


def _first_dataset(h):
    """dataset of the first patch (lowest level) that has one; None if none has."""
    time = h.times()[0]
    for lvl in h.levels(time).values():
        for patch in lvl.patches:
            for pd in patch.patch_datas.values():
                return np.asarray(pd.dataset[:])
    return None


def _by_box(h, comps):
    """{(level, lower): {comp: (primal nodes float64, cell-centred?)}} of any hierarchy."""
    time = h.times()[0]
    out = {}
    for il, lvl in h.levels(time).items():
        for p in lvl.patches:
            entry = {}
            for c in comps:
                pd = p.patch_datas[c]
                ng = int(pd.ghosts_nbr[0])
                a = np.asarray(pd.dataset[:], dtype=np.float64)
                cells = a.shape == tuple(int(n) + 2 * ng for n in p.box.shape)
                if cells:  # dual in every direction: average the 2^d cells around each node
                    a = 0.25 * (a[:-1, :-1] + a[1:, :-1] + a[:-1, 1:] + a[1:, 1:])
                    ng -= 1
                entry[c] = a[ng:-ng, ng:-ng] if ng > 0 else a
            out[(il, tuple(int(v) for v in p.box.lower))] = entry
    return out


def _rel_diff(legacy, new, comps):
    if new is None or legacy is None or isinstance(legacy, dict) or comps is None:
        return None
    a, b = _by_box(legacy, comps), _by_box(new, comps)
    if a.keys() != b.keys():
        return float("inf")
    scale = max(np.nanmax(np.abs(a[k][c])) for k in a for c in comps) or 1.0
    worst = 0.0
    for k in a:
        for c in comps:
            if a[k][c].shape != b[k][c].shape:
                return float("inf")
            worst = max(worst, float(np.nanmax(np.abs(a[k][c] - b[k][c]))))
    return worst / scale


def _describe(obj):
    """(class, times kind, dtype, shape) of a getter result."""
    if isinstance(obj, Exception):
        return (f"raises {type(obj).__name__}", None, None, None)
    if isinstance(obj, dict):
        return ("dict", None, None, None)
    ds = _first_dataset(obj)
    if ds is None:
        return (_class_name(obj), np.asarray(obj.times()).dtype.kind, None, None)
    return (_class_name(obj), np.asarray(obj.times()).dtype.kind, str(ds.dtype),
            tuple(ds.shape))


def _call(f, r):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return f(r)
    except Exception as e:  # noqa: BLE001 - any error is a row of the table
        return e


def _summary(value):
    if isinstance(value, Exception):
        return f"raises {type(value).__name__}"
    if isinstance(value, np.ndarray):
        return f"{value.dtype} {[round(float(v), 10) for v in np.sort(value)]}"
    if isinstance(value, dict):
        return f"{len(value)} keys"
    if isinstance(value, tuple):
        return f"{value[0].shape} {value[0].dtype}"
    return str(value)


def build_table():
    """{getter: (legacy class, new class, legacy times, new times, legacy dtype, new dtype,
    legacy shape, new shape, values)} and {meta getter: (legacy, new)}."""
    old = Run(H5, legacy=True)
    new = Run(H5, legacy=False)
    rows = {}
    for name, (f, comps) in GETTERS.items():
        a, b = _call(f, old), _call(f, new)
        da, db = _describe(a), _describe(b)
        values = None
        if not isinstance(a, Exception) and not isinstance(b, Exception):
            values = _rel_diff(a, b, comps if not isinstance(b, dict) else None)
        rows[name] = (da[0], db[0], da[1], db[1], da[2], db[2], da[3], db[3], values)
    meta = {}
    calls = {
        "times('B')": lambda r: r.times("B"),
        "times('EM_B')": lambda r: r.times("EM_B"),
        "all_times() keys": lambda r: r.all_times(),
        "GetDl(1)": lambda r: list(map(float, r.GetDl(1))),
        "GetDomainSize()": lambda r: list(map(float, r.GetDomainSize())),
        "GetMass('alpha')": lambda r: float(r.GetMass("alpha")),
        "GetCharge('alpha')": lambda r: float(r.GetCharge("alpha")),
        "GetMagneticFlux": lambda r: r.GetMagneticFlux(T),
        "unknown keyword GetB(foo=1)": lambda r: r.GetB(T, foo=1) and "accepted",
        "time not dumped GetB(0.5)": lambda r: r.GetB(0.5),
    }
    for name, f in calls.items():
        meta[name] = (_summary(_call(f, old)), _summary(_call(f, new)))
    flux = (_call(calls["GetMagneticFlux"], old), _call(calls["GetMagneticFlux"], new))
    return rows, meta, flux


def _ring_kind(h, comp):
    """"NaN" / "real values" for the ghost ring of the first patch of level 0."""
    pd = h.level(0, h.times()[0]).patches[0].patch_datas[comp]
    ng = int(pd.ghosts_nbr[0])
    ring = np.ones(pd.dataset.shape, bool)
    ring[tuple(slice(ng, -ng) for _ in range(pd.ndim))] = False
    return "NaN" if np.isnan(np.asarray(pd.dataset[:])[ring]).all() else "real values"


def _npatch(h):
    return str(sum(len(lvl.patches) for lvl in h.levels(h.times()[0]).values()))


def build_behaviour():
    """{behaviour: (legacy summary, new summary)} of master's per-patch API."""
    from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare

    def keys(B):
        return str(list(B.x.level(0, B.times()[0]).patches[0].patch_datas.keys()))

    calls = {
        "ghost nodes of stored quantities (Ni)": lambda r: _ring_kind(r.GetNi(T), "value"),
        "ghost nodes of B (all_primal)": lambda r: _ring_kind(r.GetB(T), "x"),
        "hierarchy_compare(Ni, Ni again)":
            lambda r: str(bool(hierarchy_compare(r.GetNi(T), r.GetNi(T)))),
        "hierarchy_compare(B, B again)":
            lambda r: str(bool(hierarchy_compare(r.GetB(T), r.GetB(T)))),
        "patch.id level 0 patch 0":
            lambda r: r.GetNi(T).level(0, r.GetNi(T).times()[0]).patches[0].id,
        "B.x patch_datas keys": lambda r: keys(r.GetB(T)),
        "B.x quantities()": lambda r: str(r.GetB(T).x.quantities()),
        "B.x.global_max('x') == B.global_max('x')":
            lambda r: str(bool(r.GetB(T).x.global_max("x") == r.GetB(T).global_max("x"))),
        "B.x['x']": lambda r: r.GetB(T).x["x"],
        "Ni == 1.0": lambda r: str(r.GetNi(T) == 1.0),
        "Ni == Ni": lambda r: (lambda n: str(n == n))(r.GetNi(T)),
        "selection_box=(5, 5, 20, 13.5) patches":
            lambda r: _npatch(r.GetB(T, selection_box=(5, 5, 20, 13.5))),
        "selection_box=(0, 0, 10, 10) patches":
            lambda r: _npatch(r.GetB(T, selection_box=(0, 0, 10, 10))),
        "GetB([])": lambda r: r.GetB([]),
    }
    old, new = Run(H5, legacy=True), Run(H5, legacy=False)
    return {name: (_summary(_call(f, old)), _summary(_call(f, new)))
            for name, f in calls.items()}


def format_table(rows, meta):
    head = ("getter", "legacy class", "new class", "times", "dtype", "dataset shape",
            "max rel diff")
    lines = [" | ".join(head)]
    for name, r in rows.items():
        val = "-" if r[8] is None else f"{r[8]:.3g}"
        lines.append(" | ".join([name, r[0], r[1], f"{r[2]}->{r[3]}", f"{r[4]}->{r[5]}",
                                 f"{r[6]}->{r[7]}", val]))
    lines.append("")
    for name, (a, b) in meta.items():
        lines.append(f"{name}: legacy {a} | new {b}")
    return "\n".join(lines)


@unittest.skipUnless(os.path.isdir(H5), "needs build/sample_data")
class DifferentialTableTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._env = os.environ.pop("PHARESEE_LEGACY", None)
        cls.rows, cls.meta, cls.flux = build_table()
        cls.behaviour = build_behaviour()
        print("\n" + format_table(cls.rows, cls.meta))
        for name, (a, b) in cls.behaviour.items():
            print(f"{name}: legacy {a} | new {b}")

    @classmethod
    def tearDownClass(cls):
        if cls._env is not None:
            os.environ["PHARESEE_LEGACY"] = cls._env

    def test_field_getters_match_the_approved_table(self):
        self.assertEqual(set(self.rows), set(APPROVED))
        for name, row in self.rows.items():
            with self.subTest(name):
                want = APPROVED[name]
                self.assertEqual(row[:8], want[:8])
                if want[8] is None:
                    self.assertIsNone(row[8])
                else:
                    self.assertLessEqual(row[8], want[8])

    def test_per_patch_behaviours_match_the_approved_table(self):
        self.assertEqual(self.behaviour, APPROVED_BEHAVIOUR)

    def test_node_current_differs_from_master_by_the_stencil_only(self):
        """GetJ (node) differs from master by the (1,2,1)/4 smoothing: much more than
        rounding, much less than the field itself."""
        self.assertGreater(self.rows["GetJ"][8], 1e-3)

    def test_other_getters_match_the_approved_table(self):
        self.assertEqual(self.meta, APPROVED_META)
        (Al, (xl, _)), (An, (xn, _)) = self.flux
        step = round((xn[1] - xn[0]) / (xl[1] - xl[0]))
        d = np.abs(Al[::step, ::step] - An).max() / np.abs(Al).max()
        self.assertLess(d, MAGNETIC_FLUX_BOUND)


if __name__ == "__main__":
    unittest.main()
