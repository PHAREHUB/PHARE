"""
This file exists independently from test_initialization.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import data, ddt, unpack

import pyphare.pharein as ph

from tests.simulator.initialize.test_init_mhd import MHDInitializationTest
from tests.simulator.initialize.test_init_hybrid import HybridInitializationTest

ph.NO_GUI()

ndim = 1
interp_orders = [1, 2, 3]


def permute_hybrid():
    return [
        dict(super_class=HybridInitializationTest, interp_order=interp_order)
        for interp_order in interp_orders
    ]


def permute_mhd():  # interp_order hax todo
    return [dict(super_class=MHDInitializationTest, interp_order=2)]


def permute(hybrid=True, mhd=False):
    return (permute_hybrid() if hybrid else []) + (permute_mhd() if mhd else [])


@ddt
class Initialization1DTest(MHDInitializationTest, HybridInitializationTest):
    @data(*permute())
    @unpack
    def test_B_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self.__class__ = super_class  # cast to super class
        self._test_B_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute())
    @unpack
    def test_bulkvel_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self.__class__ = super_class  # cast to super class
        self._test_bulkvel_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute())
    @unpack
    def test_density_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self.__class__ = super_class  # cast to super class
        self._test_density_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute())
    @unpack
    def test_density_decreases_as_1overSqrtN(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self.__class__ = super_class  # cast to super class
        self._test_density_decreases_as_1overSqrtN(ndim, **kwargs)


if __name__ == "__main__":
    unittest.main()
