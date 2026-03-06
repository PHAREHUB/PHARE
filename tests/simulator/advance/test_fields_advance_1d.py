"""
This file exists independently from test_advance.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core.box import Box1D
from pyphare.core import phare_utilities as phut

from tests.simulator.advance.test_advance_mhd import MHDAdvanceTest
from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest

ph.NO_GUI()
ndim = 1
interp_orders = [1, 2, 3]


def permute_hybrid(boxes={}):
    return [
        dict(
            super_class=HybridAdvanceTest,
            interp_order=interp_order,
            refinement_boxes=boxes,
        )
        for interp_order in interp_orders
    ]


def permute_mhd(boxes={}):  # interp_order hax todo
    return [dict(super_class=MHDAdvanceTest, interp_order=2, refinement_boxes=boxes)]


def permute(boxes={}, hybrid=True, mhd=False):
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
        *permute({}),
        *permute({"L0": [Box1D(10, 19)]}),
    )
    @unpack
    def test_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts(
        self, super_class, interp_order, **kwargs
    ):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        time_step_nbr = 3
        time_step = 0.001
        from pyphare.pharein.simulation import check_patch_size

        largest_patch_size, smallest_patch_size = check_patch_size(
            ndim, interp_order=interp_order, cells=[60] * ndim
        )
        datahier = self.getHierarchy(
            ndim,
            interp_order,
            qty="eb",
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=smallest_patch_size,
            time_step=time_step,
            time_step_nbr=time_step_nbr,
            **kwargs,
        )
        self._test_overlaped_fields_are_equal(datahier, time_step_nbr, time_step)

    @data(
        *permute(({"L0": {"B0": Box1D(10, 19)}})),
        *permute(({"L0": {"B0": Box1D(10, 14), "B1": Box1D(15, 19)}})),
        *permute(({"L0": {"B0": Box1D(6, 23)}})),
        *permute(({"L0": {"B0": Box1D(2, 12), "B1": Box1D(13, 25)}})),
        *permute(({"L0": {"B0": Box1D(5, 20)}, "L1": {"B0": Box1D(15, 19)}})),
        *permute(
            (
                {
                    "L0": {"B0": Box1D(5, 20)},
                    "L1": {"B0": Box1D(12, 38)},
                    "L2": {"B0": Box1D(30, 52)},
                }
            )
        ),
    )
    @unpack
    def test_field_coarsening_via_subcycles(self, super_class, **kwargs):
        print(f"{self._testMethodName}_{ndim}d")
        phut.cast_to(self, super_class)
        self._test_field_coarsening_via_subcycles(ndim, **kwargs)

    @unittest.skip("should change to work on moments")
    @data(  # only supports a hierarchy with 2 levels
        *permute(({"L0": [Box1D(5, 9)]})),
        *permute(({"L0": [Box1D(5, 24)]})),
        *permute(({"L0": [Box1D(5, 9), Box1D(20, 24)]})),
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
