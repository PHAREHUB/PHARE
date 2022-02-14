"""
  This file exists independently from test_advance.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box1D, nDBox
from tests.simulator.test_advance import AdvanceTestBase

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]


def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]

@ddt
class AdvanceTest(AdvanceTestBase):


    @data(
      *per_interp({}),
      *per_interp({"L0": [Box1D(10, 20)]}),
      *per_interp({"L0": [Box1D(2, 12), Box1D(13, 25)]}),
    )
    @unpack
    def test_overlapped_particledatas_have_identical_particles(self, interp_order, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_overlapped_particledatas_have_identical_particles(ndim, interp_order, refinement_boxes)



    @data(*interp_orders)
    def test_L0_particle_number_conservation(self, interp):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_L0_particle_number_conservation(ndim, interp)



    @data(
        *per_interp(({"L0": {"B0": Box1D(10, 14)}})),
    )
    @unpack
    def test_domain_particles_on_refined_level(self, interp_order, refinement_boxes):
        self._test_domain_particles_on_refined_level(ndim, interp_order, refinement_boxes)




if __name__ == "__main__":
    unittest.main()

