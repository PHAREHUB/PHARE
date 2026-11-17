"""Small helpers shared by the amr test files. Each one replaces a few lines that every file
used to write by hand; none hides what a test checks.

Cleanup  : tmp_dir, clean_env (setUp or setUpClass, undone automatically), PlotTestCase
Data     : vec, snapshot, lazy_snapshot, random_snapshot, fake_backend, fake_field
Checks   : comps, assert_all_equal, assert_matches_master (node means vs master)
Files    : h5py, float64_copy, rewrite_in_place
"""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.field import AMRField
from pyphare.pharesee.amr.store import LevelData, Snapshot
from pyphare_tests.test_pharesee.amr.testing import FakeBackend, fill

# --------------------------------------------------------------------------- cleanup


def _on_cleanup(owner, func, *args):
    """Run func(*args) after the test (owner = self) or after the class (owner = cls)."""
    if isinstance(owner, type):
        owner.addClassCleanup(func, *args)
    else:
        owner.addCleanup(func, *args)


def tmp_dir(owner):
    """New temporary directory (Path), removed after the test / class."""
    tmp = tempfile.TemporaryDirectory()
    _on_cleanup(owner, tmp.cleanup)
    return Path(tmp.name)


def clean_env(owner, *unset, **values):
    """Remove the variables ``unset`` and set ``values`` in os.environ; restored after."""
    patcher = mock.patch.dict(os.environ)
    patcher.start()
    _on_cleanup(owner, patcher.stop)
    for name in unset:
        os.environ.pop(name, None)
    os.environ.update(values)


class PlotTestCase(unittest.TestCase):
    """TestCase closing every matplotlib figure after each test."""

    def tearDown(self):
        import matplotlib.pyplot as plt

        plt.close("all")


# --------------------------------------------------------------------------- data


def vec(fx, fy, fz):
    """The component dict of a vector field, for ``testing.fill``."""
    return {"x": fx, "y": fy, "z": fz}


def snapshot(hier, f, names=None, dtype=np.float64, time=0.0):
    """In-memory Snapshot of ``f`` (callable or {name: callable}, see ``testing.fill``)."""
    return Snapshot(hier, tuple(fill(hier, f, names=names, dtype=dtype, time=time)))


def lazy_snapshot(hier, f, names=None, dtype=np.float64):
    """(Snapshot read lazily through a FakeBackend, that backend) of ``f``."""
    be = FakeBackend(fill(hier, f, names=names, dtype=dtype))
    return Snapshot(hier, tuple(LevelData.from_backend(be, g, 0.0) for g in hier.levels)), be


def random_snapshot(hier, names=("value",), dtype=np.float64, seed=0):
    """In-memory Snapshot of standard-normal values, one rng for all levels and names."""
    rng = np.random.default_rng(seed)
    return Snapshot(hier, tuple(
        LevelData.from_arrays(g, 0.0, {n: rng.standard_normal(g.nnodes).astype(dtype)
                                       for n in names}) for g in hier.levels))


def fake_backend(hier, f, names=None, times=None, dtype=np.float64, quantity="fake"):
    """FakeBackend (counts reads) serving ``f`` on ``hier``.

    times None : one dump at t=0 of ``f`` (as for ``testing.fill``)
    times list : one dump per t of ``f(t)``, so ``f`` is then a function of t returning
                 what ``testing.fill`` takes
    """
    if times is None:
        data = fill(hier, f, names=names, dtype=dtype)
    else:
        data = {float(t): fill(hier, f(t), names=names, dtype=dtype, time=t) for t in times}
    return FakeBackend(data, quantity=quantity)


def fake_field(hier, f, **kw):
    """(AMRField, FakeBackend) of ``fake_backend(hier, f, **kw)``."""
    be = fake_backend(hier, f, **kw)
    return AMRField.from_backend(be), be


# --------------------------------------------------------------------------- checks


def comps(snap, name):
    """Per-level buffers of one component: [level0 array, level1 array, ...]."""
    return [ld.buffer(name) for ld in snap.levels]


def assert_all_equal(got, want, msg=""):
    """np.testing.assert_array_equal on each pair of two equally long iterables."""
    got, want = list(got), list(want)
    assert len(got) == len(want), (len(got), len(want), msg)
    for k, (g, w) in enumerate(zip(got, want)):
        np.testing.assert_array_equal(g, w, err_msg=f"{msg} [item {k}]")


# --------------------------------------------------------------------------- node means vs master

_U64 = 2.0**-53  # unit roundoff of float64


def _gamma(n, u):
    """Bound of n successive roundings (Higham, Accuracy and Stability of Numerical
    Algorithms, 2nd ed., Lemma 3.1)."""
    return n * u / (1.0 - n * u)


def node_mean_reference(raw_pd, raw_name):
    """Exact node values of one Yee or cell-centred component, from the file's cells.

    ``raw_pd``: master's patch data read with ``all_primal=False`` (the file's values, ghosts
    included); ``raw_name``: its name ("Bz", "mhdVx"). Each node takes the k = 1, 2, 4 or 8
    cells of master's stencil (run/utils.py slices_to_primal); their sum is computed exactly
    (math.fsum) and divided by k, a power of two. Returns ``(ref, scale, k)``: the exact mean
    rounded once to float64, NaN where master leaves NaN (outside its ``inner`` slice);
    ``scale`` the mean |cell value| of the k cells, which bounds rounding errors of the sum.
    """
    import math

    from pyphare.pharesee.run.utils import slices_to_primal

    ds = np.asarray(raw_pd.dataset[:])
    ndim = raw_pd.box.ndim
    inner, chunks = slices_to_primal(raw_name, nb_ghosts=int(raw_pd.ghosts_nbr[0]), ndim=ndim)
    shape = [n + (c == "dual") for n, c in zip(ds.shape, raw_pd.centerings)]
    k = len(chunks)
    ops = np.stack([ds[c].astype(np.float64) for c in chunks])  # float32 -> float64: exact
    rows = ops.reshape(k, -1).T.tolist()
    ref, scale = np.full(shape, np.nan), np.full(shape, np.nan)
    ref[inner] = np.reshape([math.fsum(r) / k for r in rows], ops.shape[1:])
    scale[inner] = np.reshape([math.fsum(map(abs, r)) / k for r in rows], ops.shape[1:])
    return ref, scale, k


def node_mean_tolerances(ref, scale, k, file_dtype):
    """Largest error each reader's arithmetic can make at each node, vs the exact mean.

    u = unit roundoff of the file's dtype, d = log2(k) dual axes.
    New reader (amr/io/phareh5.py yee_to_nodes):
      d = 0: a copy, exact.
      d = 1: one sum in the file's dtype, then * 0.5 (exact): one rounding, u |mean|.
      d >= 2: pairwise sums in float64 (d roundings on any path), the last one rounded
      again to the file's dtype: u |mean| + (1 + u) gamma_d(2^-53) scale.
    Master (run/utils.py _compute_to_primal): each cell / k (exact), then k additions to 0 in
    float64 (the first exact) and no rounding to the file's dtype: gamma_{k-1}(2^-53) scale.
    Both get 2^-53 |ref| more, the rounding of ref itself. The relative bounds hold for
    normal numbers; a result or a cell / k in the subnormal range of the file's dtype
    (float32: below 1.2e-38, met at t=0 of the 2d_static sample) rounds instead to within
    half the smallest subnormal; eta = the smallest subnormal (half of it is 0 in float64):
    new + 2 eta (the scaling and the final rounding), master + k eta (one per cell / k, done
    in the file's dtype). Returns (new_tol, master_tol).
    """
    info = np.finfo(file_dtype)
    u, eta = float(info.eps) / 2, float(info.smallest_subnormal)
    d = int(np.log2(k))
    a = np.abs(ref)
    if d == 0:
        new = 0.0 * a
    else:
        new = u * a + (1 + u) * _gamma(d if d > 1 else 0, _U64) * scale + 2 * eta
    master = _gamma(k - 1, _U64) * scale + (k * eta if k > 1 else 0.0)
    return new + _U64 * a, master + _U64 * a


def assert_matches_master(new_ds, master_ds, file_dtype, msg, raw=None):
    """One patch component of the new reader vs master's (``all_primal=True``) dataset.

    ``raw`` = (master's ``all_primal=False`` patch data, its name) for a quantity master
    averages from Yee or cell centring to nodes; both datasets are then checked against the
    exact mean of the file's cells within each reader's rounding bound
    (``node_mean_tolerances``), not against each other: two correct orders of summation
    differ by roundoff, which is large in relative terms where the mean cancels to ~0.
    ``raw=None``: the file already holds the nodes (moments); the values must be identical.
    The new reader keeps the file's dtype.
    """
    new_ds, master_ds = np.asarray(new_ds), np.asarray(master_ds)
    assert new_ds.dtype == file_dtype, (new_ds.dtype, file_dtype, msg)
    assert new_ds.shape == master_ds.shape, (new_ds.shape, master_ds.shape, msg)
    if raw is None:
        np.testing.assert_array_equal(new_ds, master_ds, err_msg=msg)
        return
    ref, scale, k = node_mean_reference(*raw)
    tol_new, tol_master = node_mean_tolerances(ref, scale, k, file_dtype)
    for name, got, tol in (("new", new_ds, tol_new), ("master", master_ds, tol_master)):
        np.testing.assert_array_equal(np.isnan(got), np.isnan(ref), err_msg=f"{msg} {name}: NaN")
        fin = ~np.isnan(ref)
        err = np.abs(got[fin].astype(np.float64) - ref[fin])
        worst = int(np.argmax(err - tol[fin])) if err.size else 0
        assert err.size == 0 or err[worst] <= tol[fin][worst], (
            f"{msg} {name}: |error| {err[worst]:.3e} > bound {tol[fin][worst]:.3e} "
            f"(value {got[fin][worst]!r}, exact mean {ref[fin][worst]!r}, k={k})")


# getters whose all_primal=True values master averages from the file's Yee or cell centring
AVERAGED_BY_MASTER = ("GetB", "GetE", "GetMHDrho", "GetMHDV", "GetMHDP", "GetMHDrhoV",
                      "GetMHDEtot")


def check_against_master(master, new, getter, args, t, label=""):
    """``getter(t, *args)`` of the new reader vs master (``Run(path, legacy=True)``): every
    component of every patch of every level, ghost nodes included, with
    ``assert_matches_master``. Returns the number of patch components checked."""

    def key(box):
        return tuple(box.lower.tolist()) + tuple(box.upper.tolist())

    ref, got = getattr(master, getter)(t, *args), getattr(new, getter)(t, *args)
    raw = (getattr(master, getter)(t, *args, all_primal=False)
           if getter in AVERAGED_BY_MASTER else None)
    n = 0
    for il, mlvl in ref.levels(t).items():
        new_patches = {key(p.box): p for p in got.level(il, t).patches}
        assert len(new_patches) == len(mlvl.patches), (label, getter, il)
        raw_patches = {key(p.box): p for p in raw.level(il, t).patches} if raw else {}
        for mp in mlvl.patches:
            msg = f"{label} {getter}{args} t={t} L{il} {mp.box}"
            for c, mpd in mp.patch_datas.items():
                raw_pd, dtype = None, mpd.dataset.dtype  # moments: master keeps the file's
                if raw is not None:
                    pds = raw_patches[key(mp.box)].patch_datas
                    name = next(iter(pds)) if len(pds) == 1 else next(q for q in pds
                                                                      if q.endswith(c))
                    raw_pd, dtype = (pds[name], name), pds[name].dataset.dtype
                assert_matches_master(new_patches[key(mp.box)].patch_datas[c].dataset,
                                      mpd.dataset[:], dtype, f"{msg} {c}", raw=raw_pd)
                n += 1
    return n


# --------------------------------------------------------------------------- files


def h5py():
    """The h5py module (ruff bans a module-level h5py import, see ruff.toml)."""
    import h5py

    return h5py


def _float32_datasets(group):
    """Paths (from the file's root) of the float32 datasets under ``group``."""
    h5, out = h5py(), []
    group.visititems(lambda k, o: out.append(f"{group.name}/{k}".lstrip("/"))
                     if isinstance(o, h5.Dataset) and o.dtype == np.float32 else None)
    return out


def float64_copy(src_dir, dst_dir):
    """Copy the ``.h5`` run ``src_dir`` to ``dst_dir`` with every float32 dataset under /t
    rewritten as float64, as PHARE writes with PHARE_DIAG_DOUBLES (the PR CI builds).

    A float32 widened to float64 has 29 zero low bits: sums of a few of them are exact and
    every order of summation agrees, which would hide float64 roundoff. Each value v becomes
    v (1 + 2^-25 h(|v|)), h in [0, 1): full float64 mantissas, a relative change below 3e-8,
    and odd in v, so cells that cancel exactly in the run (symmetric initial conditions)
    still cancel, as in a real float64 output."""
    import shutil

    h5 = h5py()
    shutil.copytree(src_dir, dst_dir)
    for name in sorted(os.listdir(dst_dir)):
        if not name.endswith(".h5"):
            continue
        with h5.File(os.path.join(dst_dir, name), "r+") as f:
            for k in _float32_datasets(f["t"]):
                v, attrs = f[k][()].astype(np.float64), dict(f[k].attrs)
                h = np.modf(np.abs(v) * 2.0**20 * 1.6180339887)[0]
                del f[k]
                f.create_dataset(k, data=v * (1.0 + 2.0**-25 * h))
                f[k].attrs.update(attrs)


def rewrite_in_place(path, source):
    """Overwrite ``path`` with the bytes of ``source``, keeping the inode (a re-run into the
    same folder: open, truncate, write)."""
    with open(source, "rb") as src, open(path, "r+b") as dst:
        data = src.read()
        dst.truncate(0)
        dst.write(data)
