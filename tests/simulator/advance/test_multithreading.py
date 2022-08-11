"""
  This file exists independently from test_advance.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box1D, Box2D, nDBox
from tests.simulator.test_advance import AdvanceTestBase

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

interp_orders = [1, 2, 3]

def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]

@ddt
class AdvanceTest(AdvanceTestBase):

    @data(
       *per_interp(({"L0": {"B0": Box1D( 5, 20)}, "L1": {"B0": Box1D(12, 38)}})),
    )
    @unpack
    def test_multithreading_1d(self, interp_order, refinement_boxes):
        ndim = 1
        self._test_multhreading(ndim, interp_order, 5, refinement_boxes)


    @data(
       *per_interp(({"L0": {"B0": Box2D( 5, 20)}, "L1": {"B0": Box2D(12, 38)}})),
    )
    @unpack
    def test_multithreading_2d(self, interp_order, refinement_boxes):
        ndim = 2
        self._test_multhreading(ndim, interp_order, 5, refinement_boxes)


if __name__ == "__main__":
    unittest.main()
