"""
  This file exists independently from test_advance.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box2D, nDBox
from tests.simulator.test_advance import AdvanceTestBase

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]
ppc = 25

def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]

@ddt
class AdvanceTest(AdvanceTestBase):

    @data(
      *per_interp({}),
      *per_interp({"L0": [Box2D(10, 19)]}),
      *per_interp({"L0": [Box2D(8, 20)]}),
    )
    @unpack
    def test_overlaped_fields_are_equal(self, interp_order, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr=3
        time_step=0.001
        diag_outputs=f"phare_overlaped_fields_are_equal_{ndim}_{self.ddt_test_id()}"
        datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                  time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim, nbr_part_per_cell=ppc)
        self._test_overlaped_fields_are_equal(datahier, time_step_nbr, time_step)


    @data(
      *per_interp({}),
      *per_interp({"L0": [Box2D(10, 19)]}),
    )
    @unpack
    def test_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts(self, interp_order, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr=3
        time_step=0.001
        from pyphare.pharein.simulation import check_patch_size
        diag_outputs=f"phare_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts/{ndim}/{interp_order}/{self.ddt_test_id()}"
        largest_patch_size, smallest_patch_size = check_patch_size(ndim, interp_order=interp_order, cells=[60] * ndim)
        datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                  smallest_patch_size=smallest_patch_size, largest_patch_size=smallest_patch_size,
                                  time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim, nbr_part_per_cell=ppc)
        self._test_overlaped_fields_are_equal(datahier, time_step_nbr, time_step)


    @data(
       *per_interp(({"L0": {"B0": Box2D(10, 14)}})),
       *per_interp(({"L0": {"B0": Box2D(10, 14), "B1": Box2D(15, 19)}})),
       *per_interp(({"L0": {"B0": Box2D(6, 23)}})),
       *per_interp(({"L0": {"B0": Box2D( 2, 12), "B1": Box2D(13, 25)}})),
       *per_interp(({"L0": {"B0": Box2D( 5, 20)}, "L1": {"B0": Box2D(15, 19)}})),
       *per_interp(({"L0": {"B0": Box2D( 5, 20)}, "L1": {"B0": Box2D(12, 38)}, "L2": {"B0": Box2D(30, 52)} })),
    )
    @unpack
    def test_field_coarsening_via_subcycles(self, interp_order, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes, dl=.3)


    @data( # only supports a hierarchy with 2 levels
       # *per_interp(({"L0": [Box2D(0, 4)]})),   # fails?
       *per_interp(({"L0": [Box2D(10, 14)]})),
       *per_interp(({"L0": [Box2D(0, 4), Box2D(10, 14)]})),
       *per_interp(({"L0": [Box2D(0, 4), Box2D(5, 9), Box2D(10, 14)]})),
       *per_interp(({"L0": [Box2D(20, 24)]})),
       # *per_interp(({"L0": [Box2D(30, 34)]})), # fails?
    )
    @unpack
    def test_field_level_ghosts_via_subcycles_and_coarser_interpolation(self, interp_order, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_field_level_ghosts_via_subcycles_and_coarser_interpolation(ndim, interp_order, refinement_boxes)


if __name__ == "__main__":
    unittest.main()

