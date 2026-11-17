#!/usr/bin/env python3
"""Tests of the byte counts behind ``Run.plan``: sizes of the patches a read would touch."""

import unittest

import numpy as np

from pyphare.pharesee.amr.geometry import LevelGeometry
from pyphare.pharesee.amr.plan import format_bytes
from pyphare.pharesee.amr.store import LevelData
from pyphare_tests.test_pharesee.amr.testing import FakeBackend, fill, make_amr


class TestReadBytes(unittest.TestCase):
    def test_format(self):
        self.assertEqual(format_bytes(1536 * 2**20), "1.50 GiB")

    def test_sizes_missing_patches_only(self):
        hier = make_amr(2, 2, seed=1)
        src = fill(hier, {"x": lambda x, y: x, "y": lambda x, y: y}, dtype=np.float32)
        fb = FakeBackend(src)
        lds = [LevelData.from_backend(fb, g, 0.0) for g in hier.levels]
        self.assertEqual(sum(ld.read_nbytes() for ld in lds),
                         sum(g.nnodes for g in hier.levels) * 4 * 2)
        lds[0].ensure([0, 1])
        self.assertEqual(lds[0].read_nbytes(),
                         (hier.levels[0].nnodes - hier.levels[0].size[:2].sum()) * 8)
        self.assertEqual(src[0].read_nbytes(), 0)  # in memory

    def test_subset_level_uses_file_ids(self):
        # level restricted to some patches: sizes must be those of the FILE's patches
        hier = make_amr(2, 1, seed=4)
        g = hier.levels[0]
        src = fill(hier, {"x": lambda x, y: x}, dtype=np.float64)
        fb = FakeBackend(src)
        lo, up = g.file_order_boxes()
        size = np.prod(up - lo + 2, axis=1)  # nodes per patch, file order
        ids = np.argsort(size)[-3:][::-1].copy()  # the three largest patches
        sub = LevelGeometry.from_file_order(0, lo[ids], up[ids], g.dl, g.domain_cells,
                                            g.periodic)
        ld = LevelData.from_backend(fb, sub, 0.0, file_ids=ids)
        expected = int(size[ids].sum()) * 8
        wrong = int(size[:3].sum()) * 8  # what positions 0..2 of the file would give
        self.assertNotEqual(expected, wrong, "layout must have patches of different sizes")
        self.assertEqual(ld.read_nbytes(), expected)
        ld.ensure([sub.store_index[0]])
        self.assertEqual(ld.read_nbytes(), int(size[ids[1:]].sum()) * 8)
        self.assertEqual(ld.read_nbytes([sub.store_index[2]]), int(size[ids[2]]) * 8)


if __name__ == "__main__":
    unittest.main()
