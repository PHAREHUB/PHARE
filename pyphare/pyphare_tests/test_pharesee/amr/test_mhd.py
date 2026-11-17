"""``Run`` on a real MHD run: ``build/sample_data/2d_mhd/{h5,vtk}`` (gen_sample_data.py:
Orszag-Tang, 48 x 40 cells of 1/48 x 1/40, one static refined level, 3 dumps), the only
PHARE MHD output read on PR CI (the MHD functional tests are nightly).

Checks, at every dump: every ``GetMHD*`` getter (and GetB) on both formats; h5 vs vtkhdf
(stored quantities: 2.5e-7 x max|level|, the bound of test_run_amr.py; the node current
at nodes with same-level neighbours: 1e-5 x max|J|); the new h5 reader vs master's
(``Run(path, legacy=True)``) at every node of every patch, ghosts included: both average the
same cells; each is checked against the exact mean within the rounding bound of its own
arithmetic (_common.assert_matches_master), on the file and on a float64 copy of it (the
precision of the PR CI builds).
"""

import os
import unittest

import numpy as np

from pyphare.pharesee.amr import field as amr
from pyphare.pharesee.amr.ops.deriv import one_sided_mask
from pyphare.pharesee.run import Run
from pyphare_tests.test_pharesee.amr._common import (
    check_against_master,
    float64_copy,
    tmp_dir,
)
from pyphare_tests.test_pharesee.amr.sample_data import require, run_dirs
from pyphare_tests.test_pharesee.amr.test_run_amr import LegacyEnv, _level_scale

os.environ.setdefault("PHARESEE_CACHE", "0")

CASE = "2d_mhd"
GETTERS = ("GetMHDrho", "GetMHDV", "GetMHDP", "GetMHDrhoV", "GetMHDEtot", "GetB")
SCALARS = ("GetMHDrho", "GetMHDP", "GetMHDEtot")


class MHDRunTest(LegacyEnv, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require(cases=(CASE,))
        cls.h5, cls.vtk = run_dirs(CASE)

    def test_every_mhd_getter_on_both_formats(self):
        master = Run(self.h5, legacy=True)
        for r in (Run(self.vtk), Run(self.h5, legacy=False)):
            times = r.times("B")
            np.testing.assert_array_equal(times, master.times("B"))
            self.assertEqual(len(times), 3)
            self.assertEqual(set(r.all_times()), set(master.all_times()))
            np.testing.assert_array_equal(r.GetDl(), master.GetDl())
            np.testing.assert_array_equal(r.GetDl(1), master.GetDl(1))
            np.testing.assert_array_equal(r.GetDomainSize(), master.GetDomainSize())
            for t in times:
                for g in GETTERS:
                    with self.subTest(reader=repr(r), t=t, getter=g):
                        f = getattr(r, g)(t)
                        cls = amr.ScalarField if g in SCALARS else amr.VectorField
                        self.assertIsInstance(f, cls)
                        self.assertEqual(f.times().tolist(), [t])
                        self.assertEqual(f.meta.model, "mhd")
                        lds = f.levels_data()
                        self.assertEqual([ld.level for ld in lds], [0, 1])
                        for ld in lds:
                            for n in ld.names:
                                self.assertTrue(np.isfinite(ld.buffer(n)).all())

    def test_formats_agree(self):
        rh, rv = Run(self.h5, legacy=False), Run(self.vtk)
        for t in rh.times("B"):
            for g in GETTERS:
                a, b = getattr(rh, g)(t), getattr(rv, g)(t)
                for la, lb in zip(a.levels_data(), b.levels_data(), strict=True):
                    self.assertEqual(la.geom.key(), lb.geom.key())  # same boxes, same dl
                    scale = _level_scale(lb)
                    for n in la.names:
                        d = np.abs(la.buffer(n).astype(np.float64) - lb.buffer(n))
                        self.assertLessEqual(float(d.max()), 2.5e-7 * scale,
                                             f"t={t} {g} L{la.level} {n}")
            a, b = rh.GetJ(t, stencil="node"), rv.GetJ(t)
            for la, lb in zip(a.levels_data(), b.levels_data(), strict=True):
                inner = ~one_sided_mask(la.geom)
                for n in la.names:
                    d = np.abs(la.buffer(n).astype(np.float64) - lb.buffer(n))[inner]
                    self.assertLessEqual(float(d.max()), 1e-5 * _level_scale(lb),
                                         f"t={t} J L{la.level} {n}")

    def test_new_h5_matches_master_at_patch_nodes(self):
        f64 = os.path.join(tmp_dir(self), "h5")
        float64_copy(self.h5, f64)
        for label, path in (("file", self.h5), ("float64", f64)):
            rl, rn = Run(path, legacy=True), Run(path, legacy=False)
            for t in rn.times("B"):
                for g in GETTERS:
                    self.assertGreater(check_against_master(rl, rn, g, (), t, label), 0)


if __name__ == "__main__":
    unittest.main()
