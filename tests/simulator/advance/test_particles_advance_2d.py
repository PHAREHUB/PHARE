"""
This file exists independently from test_advance.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core.box import Box2D

from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest

ph.NO_GUI()

ndim = 2
interp_orders = [1, 2, 3]
ppc = 25


def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]


@ddt
class AdvanceTest2D(HybridAdvanceTest):
    @data(*interp_orders)
    def test_L0_particle_number_conservation(self, interp):
        self._test_L0_particle_number_conservation(ndim, interp, ppc=ppc)

    @data(
        *per_interp(({"L0": {"B0": Box2D(10, 14)}})),
    )
    @unpack
    def test_domain_particles_on_refined_level(self, interp_order, refinement_boxes):
        self._test_domain_particles_on_refined_level(
            ndim, interp_order, refinement_boxes, nbr_part_per_cell=ppc
        )


if __name__ == "__main__":
    unittest.main()
