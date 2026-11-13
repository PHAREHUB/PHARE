#!/usr/bin/env python3
"""pharesee.amr tests on the 3-rank sample run <sample root>_mpi3/2d_tagged/{h5,vtk}.

That run (many patches per level, real rank ids) is written by gen_sample_data.py --mpi3
under mpirun -n 3, so ctest registers this file only when configured with testMPI.
The data is required (sample_data.require): an error when absent, never a skip.
"""

import os
import unittest
import warnings

import numpy as np

from pyphare.pharesee.amr.boxes import PhysBox
from pyphare.pharesee.run import Run
from pyphare_tests.test_pharesee.amr import test_compat, test_io_phareh5, test_io_vtkhdf
from pyphare_tests.test_pharesee.amr.sample_data import SAMPLES_MPI3, require
from pyphare_tests.test_pharesee.amr.test_compat import open_field
from pyphare_tests.test_pharesee.amr.test_run_amr import (
    LegacyEnv,
    _dirs,
    _patches_by_box,
)

os.environ.setdefault("PHARESEE_CACHE", "0")  # never write the index cache under ~/.cache


class Mpi3Data:
    """Mixin: the 3-rank run is required by every class of this file."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        require(SAMPLES_MPI3)


# ============================================================================ compat


class TestThreeRanksVtkhdfVsH5(Mpi3Data, unittest.TestCase):
    """compat views on both formats of the 3-rank run (test_compat.TestVtkhdfVsH5)."""

    _compare = test_compat.TestVtkhdfVsH5._compare

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


class TestRealIdsAndCompare(Mpi3Data, unittest.TestCase):
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


# ============================================================================ readers


class TestPhareh5RealDataMpi3(test_io_phareh5.TestRealData):
    """test_io_phareh5.TestRealData on the 3-rank run."""

    ROOT, CASES = SAMPLES_MPI3, ("2d_tagged",)
    MIN_COMPARED = 25  # values at t=0 only there (issue #1526): 30 checks


class TestVtkhdfRealFilesMpi3(test_io_vtkhdf.RealFilesTest):
    """test_io_vtkhdf.RealFilesTest on the 3-rank run."""

    ROOT, CASES = SAMPLES_MPI3, ("2d_tagged",)


# ============================================================================ Run boxes


class ReadCounter:
    """Wraps a backend's read_into and read_derived: {level: sorted file ids read}."""

    def __init__(self, backend):
        self.ids = {}
        for name in ("read_into", "read_derived"):
            orig = getattr(backend, name, None)
            if orig is not None:
                setattr(backend, name, self._wrap(orig, name == "read_derived"))

    def _wrap(self, orig, derived):
        def wrapped(*args):
            _, level, file_ids = args[1:4] if derived else args[0:3]
            self.ids.setdefault(level, set()).update(int(i) for i in file_ids)
            return orig(*args)

        return wrapped


def _expected_ids(backend, time, box):
    """File ids a box read must touch, per level (PhysBox rule, a level without patches
    ends the hierarchy)."""
    meta = backend.meta()
    out = {}
    for li in sorted(backend.index(time), key=lambda li: li.level):
        dl = meta.dl0 / 2**li.level
        sel = np.nonzero(box.selects(li.lower, li.upper, dl, meta.domain_cells0 * 2**li.level))
        if len(sel[0]) == 0:
            break
        out[li.level] = {int(i) for i in sel[0]}
    return out


class BoxTest(Mpi3Data, LegacyEnv, unittest.TestCase):
    BOX = PhysBox((9.3, 4.6), (21.7, 17.2)).to_dict()

    def _runs(self):
        h5, vtk = _dirs("2d_tagged", SAMPLES_MPI3)
        return (Run(vtk), Run(h5, legacy=False))

    def test_box_reads_only_patches_touching_it(self):
        box = PhysBox.parse(self.BOX, 2)
        for r in self._runs():
            B = r.GetB(0.1, box=self.BOX)
            backend = next(iter(r._backends.values()))
            counter = ReadCounter(backend)
            B.load()
            expected = _expected_ids(backend, 0.1, box)
            self.assertEqual(counter.ids, expected)
            total = sum(len(li.lower) for li in backend.index(0.1))
            self.assertLess(sum(map(len, expected.values())), total)
            self.assertEqual(B.extent, box)

    def test_current_with_box_reads_the_widened_box(self):
        box = PhysBox.parse(self.BOX, 2)
        vtk, h5 = self._runs()
        for r, stencil in ((vtk, "node"), (h5, "node"), (h5, "yee")):
            r.times("B")  # opens the EM_B backend
            backend = next(iter(r._backends.values()))
            counter = ReadCounter(backend)
            J = r.GetJ(0.1, box=self.BOX, stencil=stencil)
            J.load()
            self.assertEqual(J.extent, box)
            self.assertEqual((J.quantity, J.stencil), ("J", stencil))
            # node stencil: one level-0 cell of halo; Yee stencil: per patch, no halo
            halo = 1 if stencil == "node" else 0
            dl0 = backend.meta().dl0
            wide = PhysBox(tuple(np.array(box.lo) - halo * dl0),
                           tuple(np.array(box.hi) + halo * dl0))
            self.assertEqual(counter.ids, _expected_ids(backend, 0.1, wide), stencil)

    def test_current_with_box_equals_whole_domain_inside(self):
        box = PhysBox.parse(self.BOX, 2)
        for r in self._runs():
            whole = _patches_by_box(r.GetJ(0.1))
            part = _patches_by_box(r.GetJ(0.1, box=self.BOX))
            self.assertLess(len(part), len(whole))
            dl0 = r.GetDl()
            checked = 0
            for (level, lower, upper), comps in part.items():
                dl = dl0 / 2**level
                x = (lower[0] + np.arange(upper[0] - lower[0] + 2)) * dl[0]
                y = (lower[1] + np.arange(upper[1] - lower[1] + 2)) * dl[1]
                inside = ((x > box.lo[0]) & (x < box.hi[0]))[:, None] & (
                    (y > box.lo[1]) & (y < box.hi[1]))[None, :]
                for c in "xyz":
                    np.testing.assert_array_equal(comps[c][inside],
                                                  whole[(level, lower, upper)][c][inside])
                checked += int(inside.sum())
            self.assertGreater(checked, 1000)

    @staticmethod
    def _boxes(new, ref, t):
        for ld in new.levels_data():
            g = ld.geom
            mine = {(tuple(lo), tuple(hi)) for lo, hi in zip(g.lower.tolist(), g.upper.tolist())}
            theirs = {(tuple(p.box.lower.tolist()), tuple(p.box.upper.tolist()))
                      for p in ref.level(g.level, f"{t:.10f}").patches}
            yield g, mine, theirs

    @staticmethod
    def _patch_list(f, t=None):
        levels = f.levels(t) if t is not None else f.levels()
        return sorted((il, tuple(int(v) for v in p.box.lower), tuple(int(v) for v in p.box.upper))
                      for il, lvl in levels.items() for p in lvl.patches)

    def test_selection_box_selects_exactly_masters_patches(self):
        """Master's rule incl. float truncation (#1451) and edge-touching patches, on both
        formats; box= keeps the PhysBox rule (fewer patches on these boxes)."""
        from pyphare.core.box import Box

        h5, vtk = _dirs("2d_tagged", SAMPLES_MPI3)
        legacy, new, rv = Run(h5, legacy=True), Run(h5, legacy=False), Run(vtk)
        boxes = [(5, 5, 20, 13.5), (0, 0, 10, 10), (9.3, 4.6, 21.7, 17.2), (0, 0, 10.5, 10.5),
                 Box((10, 5), (18, 12)), (19.75, 0, 20.25, 40), (0, 13.5, 40, 13.5)]
        fewer = 0
        for t in (0.0, 0.2):
            for sb in boxes:
                with self.subTest(t=t, sb=str(sb)):
                    ref = self._patch_list(legacy.GetB(t, selection_box=sb), f"{t:.10f}")
                    for r in (new, rv):
                        with warnings.catch_warnings(record=True) as w:
                            warnings.simplefilter("always")
                            got = r.GetB(t, selection_box=sb)
                        self.assertTrue(any(issubclass(x.category, DeprecationWarning)
                                            and "box=" in str(x.message) for x in w))
                        self.assertEqual(self._patch_list(got), ref)
                        with warnings.catch_warnings(record=True) as w:
                            warnings.simplefilter("always")
                            J = r.GetJ(t, selection_box=sb)
                        self.assertEqual(w[0].filename, __file__)  # points at this call
                        self.assertEqual(self._patch_list(J), ref)
                    if not isinstance(sb, Box):
                        lo, hi = sb[:2], sb[2:]
                        pb = self._patch_list(new.GetB(t, box=dict(x=(lo[0], hi[0]),
                                                                    y=(lo[1], hi[1]))))
                        self.assertLessEqual(set(pb), set(ref))
                        fewer += len(ref) - len(pb)
        self.assertGreater(fewer, 0)  # the two rules differ on these boxes

    def test_box_rejects_master_box_and_tuples(self):
        from pyphare.core.box import Box

        r = self._runs()[0]
        with self.assertRaisesRegex(TypeError, "box="):
            r.GetB(0.1, box=Box((1, 1), (4, 4)))
        with self.assertRaisesRegex(TypeError, "box="):
            r.GetB(0.1, box=(1.0, 1.0, 4.0, 4.0))
        with self.assertRaisesRegex(TypeError, "not both"):
            r.GetB(0.1, box=self.BOX, selection_box=(1, 1, 4, 4))

    def test_plan_reads_nothing_and_counts_box_bytes(self):
        r = self._runs()[0]
        r.GetB(0.0)  # open the backend
        backend = r._backend(next(iter(r._backends)))
        counter = ReadCounter(backend)
        whole = r.plan("B", 0.1)
        part = r.plan("B", 0.1, box=self.BOX)
        J = r.plan("J", 0.1, box=self.BOX)
        self.assertEqual(counter.ids, {})
        expected = sum(backend.nbytes(0.1, li.level, np.arange(len(li.lower)))
                       for li in backend.index(0.1))
        self.assertEqual(whole.nbytes, expected)
        self.assertLess(part.nbytes, whole.nbytes)
        self.assertLessEqual(part.nbytes, J.nbytes)
        self.assertEqual(whole.verdict, "fits")
        self.assertIn("fits", repr(whole))


if __name__ == "__main__":
    unittest.main()
