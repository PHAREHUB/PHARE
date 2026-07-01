"""
This file exists independently from test_initialization.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core import phare_utilities as phut

from tests.simulator.initialize.test_init_mhd import (
    MHDInitializationTest,
    mhd_component_cases,
)
from tests.simulator.initialize.test_init_hybrid import HybridInitializationTest

ph.NO_GUI()

ndim = 1
interp_orders = [1, 2, 3]


def mhd_split_cases(mhd_ndim):
    # 1D: component paths only (the vector-potential curl is degenerate in 1D)
    return [
        dict(name=n, model_kwargs=k, b0=b0, b1=b1, potential=False)
        for (n, k, b0, b1) in mhd_component_cases(mhd_ndim)
    ]


def permute_hybrid():
    return [
        dict(super_class=HybridInitializationTest, interp_order=interp_order)
        for interp_order in interp_orders
    ]


def permute_mhd():
    return [
        dict(super_class=MHDInitializationTest, hall=False, interp_order=1),
    ]


def permute(hybrid=True, mhd=True):
    return (permute_hybrid() if hybrid else []) + (permute_mhd() if mhd else [])


@ddt
class Initialization1DTest(MHDInitializationTest, HybridInitializationTest):
    @data(*permute())
    @unpack
    def test_B_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_B_is_as_provided_by_user(ndim, **kwargs)

    @data(*mhd_split_cases(ndim))
    def test_mhd_B0_B1_split_is_as_provided_by_user(self, case):
        print(f"\n{self._testMethodName}_{ndim}d {case['name']}")
        self._test_mhd_split_is_as_provided(
            ndim,
            case["name"],
            case["model_kwargs"],
            case["b0"],
            case["b1"],
            case["potential"],
        )


@ddt
class HybridInitialization1DTest(HybridInitializationTest):
    @data(*permute_hybrid())
    @unpack
    def test_density_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_density_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute_hybrid())
    @unpack
    def test_bulkvel_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_bulkvel_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute_hybrid())
    @unpack
    def test_density_decreases_as_1overSqrtN(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_density_decreases_as_1overSqrtN(ndim, **kwargs)


if __name__ == "__main__":
    unittest.main()
