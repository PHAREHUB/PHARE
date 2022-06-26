"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box1D, nDBox
from tests.simulator.test_initialization import InitializationTest

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]


@ddt
class InitializationTest(InitializationTest):

    @data(*interp_orders)
    def test_B_is_as_provided_by_user(self, interp_order):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_B_is_as_provided_by_user(ndim, interp_order)

    @data(*interp_orders)
    def test_bulkvel_is_as_provided_by_user(self, interp_order):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_bulkvel_is_as_provided_by_user(ndim, interp_order)

    @data(*interp_orders)
    def test_density_is_as_provided_by_user(self, interp_order):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_density_is_as_provided_by_user(ndim, interp_order)

    @data(*interp_orders)
    def test_density_decreases_as_1overSqrtN(self, interp_order):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_density_decreases_as_1overSqrtN(ndim, interp_order)


if __name__ == "__main__":
    unittest.main()
