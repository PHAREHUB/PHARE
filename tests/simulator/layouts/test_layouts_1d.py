#!/usr/bin/env python3
"""
1D tile-based particle layout tests - see tests/simulator/layouts/__init__.py
for the shared comparison machinery.
"""

import unittest
from ddt import data, ddt, unpack

from tests.simulator.layouts import (
    permute,
    LayoutInitL0Test as _LayoutInitL0Test,
    LayoutInitL1Test as _LayoutInitL1Test,
    LayoutAdvanceL0Test as _LayoutAdvanceL0Test,
    LayoutAdvanceL1Test as _LayoutAdvanceL1Test,
    LayoutAdvanceL2Test as _LayoutAdvanceL2Test,
)

ndim_list = [1]
interp_orders = [1]


@ddt
class LayoutInitL0Test(_LayoutInitL0Test):
    @data(*permute(ndim_list, interp_orders))
    @unpack
    def test_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, cmp_layout)


@ddt
class LayoutInitL1Test(_LayoutInitL1Test):
    @data(*permute(ndim_list, interp_orders))
    @unpack
    def test_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, cmp_layout)


@ddt
class LayoutAdvanceL0Test(_LayoutAdvanceL0Test):
    @data(*permute(ndim_list, interp_orders))
    @unpack
    def test_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, cmp_layout)


@ddt
class LayoutAdvanceL1Test(_LayoutAdvanceL1Test):
    @data(*permute(ndim_list, interp_orders))
    @unpack
    def test_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, cmp_layout)


@ddt
class LayoutAdvanceL2Test(_LayoutAdvanceL2Test):
    @data(*permute(ndim_list, interp_orders))
    @unpack
    def test_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, cmp_layout)


if __name__ == "__main__":
    from pyphare.simulator.simulator import startMPI

    startMPI()
    unittest.main()
