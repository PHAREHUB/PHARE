import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box2D, nDBox
from tests.simulator.test_advance import AdvanceTestBase

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]


@ddt
class AdvanceTest(AdvanceTestBase):

    @data(
      {"L0": [Box2D(10, 19)]},
      {"L0": [Box2D(8, 20)]},
    )
    def test_overlaped_fields_are_equal(self, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr=3
        time_step=0.001
        diag_outputs=f"phare_overlaped_fields_are_equal_{ndim}_{self.ddt_test_id()}"
        for interp_order in [1, 2, 3]:
            datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                      time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim)
            self._test_overlaped_fields_are_equal(time_step, time_step_nbr, datahier)


    @data({"L0": [Box2D(10, 19)]})
    def test_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts(self, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr=3
        time_step=0.001
        from pyphare.pharein.simulation import check_patch_size
        diag_outputs=f"phare_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts_{ndim}_{self.ddt_test_id()}"
        for interp_order in [1, 2, 3]:
            largest_patch_size, smallest_patch_size = check_patch_size(ndim, interp_order=interp_order, cells=[30] * ndim)
            datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                      smallest_patch_size=smallest_patch_size, largest_patch_size=smallest_patch_size,
                                      time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim)
            self._test_overlaped_fields_are_equal(time_step, time_step_nbr, datahier)


    @data( # POSSIBLE PROBLEMO!
       ({"L0": {"B0": Box2D(10, 19)}}),
       ({"L0": {"B0": Box2D(10, 14), "B1": Box2D(15, 19)}}),
       ({"L0": {"B0": Box2D(6, 23)}}),
       ({"L0": {"B0": Box2D( 2, 12), "B1": Box2D(13, 25)}}),
       ({"L0": {"B0": Box2D( 5, 20)}, "L1": {"B0": Box2D(15, 19)}}),
       # ({"L0": {"B0": Box2D( 5, 20)}, "L1": {"B0": Box2D(12, 38)}, "L2": {"B0": Box2D(30, 52)} }), # particle cell jump > 1
    )
    def test_field_coarsening_via_subcycles(self, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        for interp_order in [1, 2, 3]:
            self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)


    ## needs update in test_refine_field.py::refine
    # @data( # only supports a hierarchy with 2 levels
    #    ({"L0": [Box2D(5, 9)]}),
    #    ({"L0": [Box2D(5, 24)]}),
    #    ({"L0": [Box2D(5, 9), Box2D(20, 24)]}),
    # )
    # def test_field_level_ghosts_via_subcycles_and_coarser_interpolation(self, refinement_boxes):
    #     for interp in [1, 2, 3]:
    #         self._test_field_level_ghosts_via_subcycles_and_coarser_interpolation(ndim, interp, refinement_boxes)


if __name__ == "__main__":
    unittest.main()

