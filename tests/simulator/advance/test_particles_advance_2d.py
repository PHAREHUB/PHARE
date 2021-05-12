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
      *per_interp({"L0": [Box2D(5, 9), Box2D(10, 14)]}),
    )
    @unpack
    def test_overlapped_particledatas_have_identical_particles(self, interp_order, refinement_boxes):
        self._test_overlapped_particledatas_have_identical_particles(
            ndim, interp_order, refinement_boxes, ppc=ppc, cells=40, smallest_patch_size=5, largest_patch_size=20)


    def test_L0_particle_number_conservation(self):
        self._test_L0_particle_number_conservation(ndim, ppc=ppc)


if __name__ == "__main__":
    unittest.main()

