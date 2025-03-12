"""
  This file exists independently from test_advance.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest

import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box3D

from tests.simulator.test_advance import AdvanceTestBase

matplotlib.use("Agg")  # for systems without GUI

ndim = 3
interp_orders = [1, 2, 3]
ppc = 5


def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]


@ddt
class AdvanceTest(AdvanceTestBase):
    @data(
        *per_interp({}),
        *per_interp({"L0": [Box3D(10, 19)]}),
        *per_interp({"L0": [Box3D(5, 9), Box3D(10, 14)]}),
    )
    @unpack
    def test_overlapped_particledatas_have_identical_particles(
        self, interp_order, refinement_boxes
    ):
        self._test_overlapped_particledatas_have_identical_particles(
            ndim,
            interp_order,
            refinement_boxes,
            ppc=ppc,
            cells=40,
            largest_patch_size=20,
        )

    @data(*interp_orders)
    def test_L0_particle_number_conservation(self, interp):
        self._test_L0_particle_number_conservation(ndim, interp, ppc=ppc, cells=30)


if __name__ == "__main__":
    unittest.main()
