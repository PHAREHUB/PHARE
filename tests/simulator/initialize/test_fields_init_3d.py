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
    mhd_component_cases,
    mhd_potential_cases,
)
from tests.simulator.initialize.test_init_hybrid import HybridInitializationTest

ph.NO_GUI()

ndim = 3
interp_orders = [1, 2, 3]
cells = 20


def permute_hybrid():
    return [
        dict(super_class=HybridInitializationTest, interp_order=interp_order)
        for interp_order in interp_orders
    ]


def permute_mhd():
    return [dict(super_class=MHDInitializationTest, hall=False)]


def permute(hybrid=True, mhd=False):
    return (permute_hybrid() if hybrid else []) + (permute_mhd() if mhd else [])


def mhd_split_cases(mhd_ndim, include_potential):
    out = [
        dict(name=n, model_kwargs=k, b0=b0, b1=b1, potential=False)
        for (n, k, b0, b1) in mhd_component_cases(mhd_ndim)
    ]
    if include_potential:
        out += [
            dict(name=n, model_kwargs=k, b0=b0, b1=b1, potential=True)
            for (n, k, b0, b1) in mhd_potential_cases(mhd_ndim)
        ]
    return out


@ddt
class Initialization3DTest(MHDInitializationTest, HybridInitializationTest):
    @data(*permute())
    @unpack
    def test_B_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_B_is_as_provided_by_user(ndim, cells=cells, **kwargs)

    @data(*permute())
    @unpack
    def test_bulkvel_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_bulkvel_is_as_provided_by_user(ndim, cells=cells, **kwargs)

    @data(*permute())
    @unpack
    def test_density_is_as_provided_by_user(self, super_class, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_density_is_as_provided_by_user(ndim, cells=cells, **kwargs)

    @data(*permute())
    @unpack
    def test_density_decreases_as_1overSqrtN(self, super_class, **kwargs):
        print(f"\n{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_density_decreases_as_1overSqrtN(
            ndim, nbr_particles=np.asarray([20, 50, 75]), cells=10, **kwargs
        )

    @data(*mhd_split_cases(ndim, include_potential=True))
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

    def test_mhd_potential_init_is_divergence_free(self):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_mhd_potential_is_divergence_free(ndim)

    def test_mhd_time_dependent_uniform_B0(self):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_mhd_time_dependent_uniform_B0(ndim)

    def test_mhd_time_dependent_potential_B0(self):
        print(f"\n{self._testMethodName}_{ndim}d")
        self._test_mhd_time_dependent_potential_B0(ndim)


if __name__ == "__main__":
    unittest.main()
