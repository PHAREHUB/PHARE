import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box2D, nDBox
from tests.simulator.test_initialization import InitializationTest

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]
ppc = 100

@ddt
class InitializationTest(InitializationTest):

    def test_B_is_as_provided_by_user(self):
        print(f"\n{self._testMethodName}_{ndim}d")
        for interp_order in interp_orders:
            self._test_B_is_as_provided_by_user(ndim, interp_order, nbr_part_per_cell=ppc)

    def test_bulkvel_is_as_provided_by_user(self):
        print(f"\n{self._testMethodName}_{ndim}d")
        for interp_order in interp_orders:
            self._test_bulkvel_is_as_provided_by_user(ndim, interp_order)

    # def test_density_is_as_provided_by_user(self): ## need to verify ppc / results - NEEDS RAM
    #     print(f"\n{self._testMethodName}_{ndim}d")
    #     for interp_order in interp_orders:
    #         self._test_density_is_as_provided_by_user(ndim, interp_order)

    def test_density_decreases_as_1overSqrtN(self):
        print(f"\n{self._testMethodName}_{ndim}d")
        for interp_order in interp_orders:
            self._test_density_decreases_as_1overSqrtN(ndim, interp_order)


if __name__ == "__main__":
    unittest.main()