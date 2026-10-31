"""
This file exists independently from test_advance.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import itertools
import unittest

from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core import phare_utilities as phut
from pyphare.core.box import Box1D
from pyphare.cpp import supported_model_layouts, supported_mhd_model_layouts
from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest
from tests.simulator.advance.test_advance_mhd import MHDAdvanceTest

ph.NO_GUI()
ndim = 1
interp_orders = [1, 2, 3]


def permute_hybrid(boxes=None):
    if boxes is None:
        boxes = {}
    return [
        dict(
            super_class=HybridAdvanceTest,
            interp_order=interp_order,
            refinement_boxes=boxes,
            model_layout=layout,
        )
        for interp_order, layout in itertools.product(
            interp_orders, supported_model_layouts()
        )
    ]


def permute_mhd(boxes=None):
    if boxes is None:
        boxes = {}
    return [
        dict(
            super_class=MHDAdvanceTest,
            hall=False,
            refinement_boxes=boxes,
            interp_order=None,
            model_layout=layout,
        )
        for layout in supported_mhd_model_layouts()
    ]


def permute(boxes=None, hybrid=True, mhd=False):
    return (permute_hybrid(boxes) if hybrid else []) + (
        permute_mhd(boxes) if mhd else []
    )


@ddt
class AdvanceTest1D(HybridAdvanceTest, MHDAdvanceTest):
    @data(
        *permute({}),
        *permute({"L0": [Box1D(10, 19)]}),
        *permute({"L0": [Box1D(8, 20)]}),
    )
    @unpack
    def test_overlaped_fields_are_equal(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        time_step_nbr = 3
        time_step = 0.001

        datahier = self.getHierarchy(
            ndim, qty="eb", time_step=time_step, time_step_nbr=time_step_nbr, **kwargs
        )
        self._test_overlaped_fields_are_equal(datahier, time_step_nbr, time_step)

    @data(
        *permute({"L0": {"B0": Box1D(10, 19)}}),
        *permute({"L0": {"B0": Box1D(10, 14), "B1": Box1D(15, 19)}}),
        *permute({"L0": {"B0": Box1D(6, 23)}}),
        *permute({"L0": {"B0": Box1D(2, 12), "B1": Box1D(13, 25)}}),
        *permute({"L0": {"B0": Box1D(5, 20)}, "L1": {"B0": Box1D(15, 19)}}),
        *permute(
            {
                "L0": {"B0": Box1D(5, 20)},
                "L1": {"B0": Box1D(12, 38)},
                "L2": {"B0": Box1D(30, 52)},
            }
        ),
    )
    @unpack
    def test_field_coarsening_via_subcycles(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_field_coarsening_via_subcycles(ndim, **kwargs)

    @unittest.skip("should change to work on moments")
    @data(  # only supports a hierarchy with 2 levels
        *permute({"L0": [Box1D(5, 9)]}),
        *permute({"L0": [Box1D(5, 24)]}),
        *permute({"L0": [Box1D(5, 9), Box1D(20, 24)]}),
    )
    @unpack
    def test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
        self, super_class, **kwargs
    ):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
            ndim, **kwargs
        )


if __name__ == "__main__":
    unittest.main()
