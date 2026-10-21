"""
This file exists independently from test_initialization.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import unittest
import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core import phare_utilities as phut

from tests.simulator.initialize.test_init_mhd import (
    MHDInitializationTest,
    mhd_split_cases,
)
from tests.simulator.initialize.test_init_hybrid import HybridInitializationTest

ph.NO_GUI()

ndim = 2
interp_orders = [1, 2, 3]


def permute_hybrid():
    return [
        dict(super_class=HybridInitializationTest, interp_order=interp_order)
        for interp_order in interp_orders
    ]


def permute_mhd():
    return [
        dict(super_class=MHDInitializationTest, hall=False, interp_order=None),
    ]


def permute(hybrid=True, mhd=False):
    return (permute_hybrid() if hybrid else []) + (permute_mhd() if mhd else [])


@ddt
class Initialization2DTest(MHDInitializationTest, HybridInitializationTest):
    @data(*permute())
    @unpack
    def test_B_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_B_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute())
    @unpack
    def test_bulkvel_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_bulkvel_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute())
    @unpack
    def test_density_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_density_is_as_provided_by_user(ndim, **kwargs)

    @data(*permute())
    @unpack
    def test_density_decreases_as_1overSqrtN(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_density_decreases_as_1overSqrtN(
            ndim, nbr_particles=np.asarray([50, 500, 1000, 2222]), cells=100, **kwargs
        )

    @data(*mhd_split_cases(ndim))
    def test_mhd_B0_B1_split_is_as_provided_by_user(self, case):
        print(f"\n{self._testMethodName}_{ndim}d {case['name']}")
        phut.cast_to(self, MHDInitializationTest)
        self._test_mhd_split_is_as_provided(
            ndim,
            case["name"],
            case["model_kwargs"],
            case["external_field"],
            case["b0"],
            case["b1"],
        )

    def test_mhd_time_dependent_uniform_B0(self):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, MHDInitializationTest)
        self._test_mhd_time_dependent_uniform_B0(ndim)

    @data("TVDRK2", "TVDRK3", "SSPRK4_5")
    def test_mhd_B0_follows_the_rk_stages(self, mhd_timestepper):
        print(f"\n{self._testMethodName}_{ndim}d {mhd_timestepper}")
        phut.cast_to(self, MHDInitializationTest)
        self._test_mhd_B0_follows_the_rk_stages(ndim, mhd_timestepper)


if __name__ == "__main__":
    unittest.main()
