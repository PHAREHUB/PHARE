"""
This file exists independently from test_advance.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core.box import Box2D
from pyphare.core import phare_utilities as phut

from tests.simulator.advance.test_advance_mhd import MHDAdvanceTest
from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest

ph.NO_GUI()
ndim = 2
interp_orders = [1, 2, 3]
ppc = 25


def permute_hybrid(boxes={}):
    return [
        dict(
            super_class=HybridAdvanceTest,
            interp_order=interp_order,
            refinement_boxes=boxes,
            nbr_part_per_cell=ppc,
        )
        for interp_order in interp_orders
    ]


def permute_mhd(boxes={}):
    return [dict(super_class=MHDAdvanceTest, hall=False, refinement_boxes=boxes)]


def permute(boxes={}, hybrid=True, mhd=False):
    return (permute_hybrid(boxes) if hybrid else []) + (
        permute_mhd(boxes) if mhd else []
    )


@ddt
class AdvanceTest2D(HybridAdvanceTest, MHDAdvanceTest):
    @data(
        *permute({}),
        *permute({"L0": [Box2D(10, 19)]}),
        *permute({"L0": [Box2D(8, 20)]}),
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

    @unittest.skip("maybe invalid now?")
    @data(
        *permute(({"L0": {"B0": Box2D(10, 14)}})),
        *permute(({"L0": {"B0": Box2D(10, 14), "B1": Box2D(15, 19)}})),
        *permute(({"L0": {"B0": Box2D(6, 23)}})),
        *permute(({"L0": {"B0": Box2D(2, 12), "B1": Box2D(13, 25)}})),
        *permute(({"L0": {"B0": Box2D(5, 20)}, "L1": {"B0": Box2D(15, 19)}})),
        *permute(
            (
                {
                    "L0": {"B0": Box2D(5, 20)},
                    "L1": {"B0": Box2D(12, 38)},
                    "L2": {"B0": Box2D(30, 52)},
                }
            )
        ),
    )
    @unpack
    def test_field_coarsening_via_subcycles(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_field_coarsening_via_subcycles(ndim, dl=0.3, **kwargs)

    @unittest.skip("should change to work with moments")
    @data(  # only supports a hierarchy with 2 levels
        *permute(({"L0": [Box2D(0, 4)]})),
        *permute(({"L0": [Box2D(10, 14)]})),
        *permute(({"L0": [Box2D(0, 4), Box2D(10, 14)]})),
        *permute(({"L0": [Box2D(0, 4), Box2D(5, 9), Box2D(10, 14)]})),
        *permute(({"L0": [Box2D(20, 24)]})),
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
