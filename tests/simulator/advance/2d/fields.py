import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box2D, nDBox
from tests.simulator.test_advance import AdvanceTest

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]


@ddt
class Advance2dTest(AdvanceTest):


    @data(
      {"L0": [Box2D(10, 19)]},
      {"L0": [Box2D(8, 20)]},
    )
    def test_overlaped_fields_are_equal(self, refinement_boxes):
        time_step_nbr=3
        time_step=0.001
        diag_outputs=f"phare_overlaped_fields_are_equal_{ndim}_{self.ddt_test_id()}"
        for interp_order in [1, 2, 3]:
            datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                      time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim)
            self._test_overlaped_fields_are_equal(time_step, time_step_nbr, datahier)


    @data({"L0": [Box2D(10, 19)]})
    def test_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts(self, refinement_boxes):
        time_step_nbr=3
        time_step=0.001
        from pyphare.pharein.simulation import check_patch_size
        diag_outputs=f"phare_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts_{ndim}_{self.ddt_test_id()}"
        for interp_order in [1, 2, 3]:
            largest_patch_size, smallest_patch_size = check_patch_size(interp_order=interp_order, cells=[30] * ndim)
            datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                      smallest_patch_size=smallest_patch_size, largest_patch_size=smallest_patch_size,
                                      time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim)
            self._test_overlaped_fields_are_equal(time_step, time_step_nbr, datahier)


if __name__ == "__main__":
    unittest.main()

