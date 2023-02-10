"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest

import matplotlib
from ddt import data, ddt

from tests.simulator.test_initialization import InitializationTest

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]
ppc = 100

@ddt
class InitializationTest(InitializationTest):

    @data(*interp_orders)
    def test_B_is_as_provided_by_user(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_B_is_as_provided_by_user(ndim, interp_order, nbr_part_per_cell=ppc)

    @data(*interp_orders)
    def test_bulkvel_is_as_provided_by_user(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_bulkvel_is_as_provided_by_user(ndim, interp_order)

    @data(*interp_orders)
    def test_density_is_as_provided_by_user(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_density_is_as_provided_by_user(ndim, interp_order)

    @data(*interp_orders)
    def test_density_decreases_as_1overSqrtN(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_density_decreases_as_1overSqrtN(ndim, interp_order)


if __name__ == "__main__":
    unittest.main()
