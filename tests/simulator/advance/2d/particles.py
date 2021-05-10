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
      {"L0": [Box2D(10, 20)]},
      {"L0": [Box2D(2, 12), Box2D(13, 25)]},
    )
    def test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(self, refinement_boxes):
        for interp_order in [1, 2, 3]:
            self._test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(ndim, interp_order, refinement_boxes)


    @data(
      {"L0": [Box2D(10, 20)]},
      {"L0": [Box2D(2, 12), Box2D(13, 25)]},
    )
    def test_overlapped_particledatas_have_identical_particles(self, refinement_boxes):
        for interp_order in [1, 2, 3]:
            self._test_overlapped_particledatas_have_identical_particles(ndim, interp_order, refinement_boxes)


    def test_L0_particle_number_conservation(self):
        self._test_L0_particle_number_conservation(ndim)


if __name__ == "__main__":
    unittest.main()

