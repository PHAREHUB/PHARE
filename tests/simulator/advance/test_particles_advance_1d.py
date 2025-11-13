"""
  This file exists independently from test_advance.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest

import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box1D

from tests.simulator.test_advance import AdvanceTestBase

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]
#


def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]


@ddt
class AdvanceTest(AdvanceTestBase):
    @data(*interp_orders)
    def test_L0_particle_number_conservation(self, interp):
        self._test_L0_particle_number_conservation(ndim, interp)

    @data(
        *per_interp(({"L0": {"B0": Box1D(10, 14)}})),
    )
    @unpack
    def test_domain_particles_on_refined_level(self, interp_order, refinement_boxes):
        self._test_domain_particles_on_refined_level(
            ndim, interp_order, refinement_boxes
        )


if __name__ == "__main__":
    unittest.main()
