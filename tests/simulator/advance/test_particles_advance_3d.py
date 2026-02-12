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
    @data(*interp_orders)
    def test_L0_particle_number_conservation(self, interp):
        self._test_L0_particle_number_conservation(ndim, interp, ppc=ppc, cells=20)

    @data(
        *per_interp(({"L0": {"B0": Box3D(10, 14)}})),
    )
    @unpack
    def test_domain_particles_on_refined_level(self, interp_order, refinement_boxes):
        self._test_domain_particles_on_refined_level(
            ndim, interp_order, refinement_boxes, nbr_part_per_cell=ppc, cells=20
        )


if __name__ == "__main__":
    unittest.main()
