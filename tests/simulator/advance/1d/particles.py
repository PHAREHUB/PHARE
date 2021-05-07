import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box1D, nDBox
from tests.simulator.test_advance import AdvanceTest

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]


@ddt
class Advance2dTest(AdvanceTest):

    @data(
      {"L0": [Box1D(10, 20)]},
      {"L0": [Box1D(2, 12), Box1D(13, 25)]},
    )
    def test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(self, refinement_boxes):
        for interp_order in [1, 2, 3]:
            self._test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(ndim, interp_order, refinement_boxes)


    @data(
      {"L0": [Box1D(10, 20)]},
      {"L0": [Box1D(2, 12), Box1D(13, 25)]},
    )
    def test_overlapped_particledatas_have_identical_particles(self, refinement_boxes):
        for interp_order in [1, 2, 3]:
            self._test_overlapped_particledatas_have_identical_particles(ndim, interp_order, refinement_boxes)


    @data(
       ({"L0": {"B0": Box1D(10, 19)}}),
       ({"L0": {"B0": Box1D(10, 14), "B1": Box1D(15, 19)}}),
       ({"L0": {"B0": Box1D(6, 23)}}),
       ({"L0": {"B0": Box1D( 2, 12), "B1": Box1D(13, 25)}}),
       ({"L0": {"B0": Box1D( 5, 20)}, "L1": {"B0": Box1D(15, 19)}}),
       ({"L0": {"B0": Box1D( 5, 20)}, "L1": {"B0": Box1D(12, 38)}, "L2": {"B0": Box1D(30, 52)} }),
    )
    def test_field_coarsening_via_subcycles(self, refinement_boxes):
        for interp_order in [1, 2, 3]:
            self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)



    @data( # only supports a hierarchy with 2 levels
       ({"L0": [Box1D(5, 9)]}),
       ({"L0": [Box1D(5, 24)]}),
       ({"L0": [Box1D(5, 9), Box1D(20, 24)]}),
    )
    def test_field_level_ghosts_via_subcycles_and_coarser_interpolation(self, refinement_boxes):
        for interp in [1, 2, 3]:
            self._test_field_level_ghosts_via_subcycles_and_coarser_interpolation(ndim, interp, refinement_boxes)


    def test_L0_particle_number_conservation(self):
        self._test_L0_particle_number_conservation(ndim)


if __name__ == "__main__":
    unittest.main()

