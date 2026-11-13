"""Small helpers shared by the amr test files. Each one replaces a few lines that every file
used to write by hand; none hides what a test checks.

Cleanup  : tmp_dir, clean_env (setUp or setUpClass, undone automatically), PlotTestCase
Data     : vec, snapshot, lazy_snapshot, random_snapshot, fake_backend, fake_field
Checks   : comps, assert_all_equal
Files    : h5py, rewrite_in_place
"""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from pyphare.pharesee.amr.field import AMRField
from pyphare.pharesee.amr.store import LevelData, Snapshot
from pyphare.pharesee.amr.testing import FakeBackend, fill

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


# --------------------------------------------------------------------------- files


def h5py():
    """The h5py module (ruff bans a module-level h5py import, see ruff.toml)."""
    import h5py

    return h5py


def rewrite_in_place(path, source):
    """Overwrite ``path`` with the bytes of ``source``, keeping the inode (a re-run into the
    same folder: open, truncate, write)."""
    with open(source, "rb") as src, open(path, "r+b") as dst:
        data = src.read()
        dst.truncate(0)
        dst.write(data)
