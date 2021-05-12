import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box2D, nDBox
from tests.simulator.test_initialization import InitializationTest

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]
ppc = 10

def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]

@ddt
class InitializationTest(InitializationTest):

    @data(*interp_orders)
    def test_nbr_particles_per_cell_is_as_provided(self, interp_order):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_nbr_particles_per_cell_is_as_provided(ndim, interp_order)


    @data(
        *per_interp(({"L0": {"B0": Box2D(10, 14)}})),
        *per_interp(({"L0": {"B0": Box2D(10, 14)}, "L1": {"B0": Box2D(22, 26)}})),
        *per_interp(({"L0": {"B0": Box2D(2, 6), "B1": Box2D(7, 11)}})),
    )
    @unpack
    def test_levelghostparticles_have_correct_split_from_coarser_particle(
        self, interp_order, refinement_boxes
    ):
        print(f"\n{self._testMethodName}_{ndim}d")
        now = self.datetime_now()
        self._test_levelghostparticles_have_correct_split_from_coarser_particle(
            self.getHierarchy(
                interp_order,
                refinement_boxes,
                "particles",
                ndim=ndim,
                cells=30, nbr_part_per_cell=ppc,
                diag_outputs=f"phare_outputs/test_levelghost/{ndim}/{interp_order}/{self.ddt_test_id()}",
            )
        )
        print(f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds")


    @data(
        *per_interp(({"L0": {"B0": Box2D(10, 14)}})),
        *per_interp(({"L0": {"B0": Box2D(5, 20)}, "L1": {"B0": Box2D(15, 35)}})),
        *per_interp(({"L0": {"B0": Box2D(2, 12), "B1": Box2D(13, 25)}})),
    )
    @unpack
    def test_domainparticles_have_correct_split_from_coarser_particle(
        self, interp_order, refinement_boxes
    ):
        print(f"\n{self._testMethodName}_{ndim}d")
        now = self.datetime_now()
        self._test_domainparticles_have_correct_split_from_coarser_particle(
            ndim, interp_order, refinement_boxes, nbr_part_per_cell=ppc
        )
        print(f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds")



    @data({"cells": 40, "smallest_patch_size": 20, "largest_patch_size": 20, "nbr_part_per_cell" : ppc})
    def test_no_patch_ghost_on_refined_level_case(self, simInput):
        print(f"\n{self._testMethodName}_{ndim}d")
        now = self.datetime_now()
        self._test_patch_ghost_on_refined_level_case(ndim, False, **simInput)
        print(f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds")

    @data({"cells": 40, "smallest_patch_size": 5, "largest_patch_size": 5, "nbr_part_per_cell" : ppc})
    def test_has_patch_ghost_on_refined_level_case(self, simInput):
        print(f"\n{self._testMethodName}_{ndim}d")
        now = self.datetime_now()
        self._test_patch_ghost_on_refined_level_case(ndim, True, **simInput)
        print(f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds")


if __name__ == "__main__":
    unittest.main()
