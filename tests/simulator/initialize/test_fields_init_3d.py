"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest

import numpy as np
import matplotlib
from ddt import data, ddt

from tests.simulator.test_initialization import InitializationTest

matplotlib.use("Agg")  # for systems without GUI

ndim = 3
interp_orders = [1, 2, 3]
ppc, cells = 10, 20


@ddt
class Initialization3DTest(InitializationTest):
    @data(*interp_orders)
    def test_B_is_as_provided_by_user(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_B_is_as_provided_by_user(ndim, interp_order, ppc=ppc, cells=cells)

    @data(*interp_orders)
    def test_bulkvel_is_as_provided_by_user(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_bulkvel_is_as_provided_by_user(
            ndim, interp_order, ppc=ppc, cells=cells
        )

    @data(*interp_orders)
    def test_density_is_as_provided_by_user(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_density_is_as_provided_by_user(ndim, interp_order, cells=cells)

    @data(*interp_orders)  # uses too much RAM - to isolate somehow
    def test_density_decreases_as_1overSqrtN(self, interp_order):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_density_decreases_as_1overSqrtN(
            ndim, interp_order, np.asarray([20, 50, 75]), cells=10
        )


if __name__ == "__main__":
    unittest.main()
