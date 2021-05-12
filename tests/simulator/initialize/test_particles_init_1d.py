import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box1D, nDBox
from tests.simulator.test_initialization import InitializationTest

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]
ppc = 25

def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]

@ddt
class InitializationTest(InitializationTest):

    @data(*interp_orders)
    def test_nbr_particles_per_cell_is_as_provided(self, interp_order):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_nbr_particles_per_cell_is_as_provided(ndim, interp_order)


    @data(
        *per_interp(({"L0": {"B0": Box1D(10, 14)}})),
        *per_interp(({"L0": {"B0": Box1D(5, 20)}, "L1": {"B0": Box1D(15, 35)}})),
        *per_interp(({"L0": {"B0": Box1D(2, 12), "B1": Box1D(13, 25)}})),
    )
    @unpack
    def test_levelghostparticles_have_correct_split_from_coarser_particle(
        self, interp_order, refinement_boxes
    ):
        print(f"{self._testMethodName}_{ndim}d")

        self._test_levelghostparticles_have_correct_split_from_coarser_particle(
            self.getHierarchy(
                interp_order,
                refinement_boxes,
                "particles",
                ndim=ndim,
                cells=30,
                diag_outputs=f"phare_outputs/test_levelghost/{ndim}/{interp_order}/{self.ddt_test_id()}",
            )
        )



    @data(
        *per_interp(({"L0": {"B0": Box1D(10, 14)}})),
        *per_interp(({"L0": {"B0": Box1D(5, 20)}, "L1": {"B0": Box1D(15, 35)}})),
        *per_interp(({"L0": {"B0": Box1D(2, 12), "B1": Box1D(13, 25)}})),
    )
    @unpack
    def test_domainparticles_have_correct_split_from_coarser_particle(
        self, interp_order, refinement_boxes
    ):
        print(f"{self._testMethodName}_{ndim}d")

        self._test_domainparticles_have_correct_split_from_coarser_particle(
            ndim, interp_order, refinement_boxes
        )




    @data({"cells": 40, "smallest_patch_size": 20, "largest_patch_size": 20})
    def test_no_patch_ghost_on_refined_level_case(self, simInput):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_patch_ghost_on_refined_level_case(ndim, False, **simInput)

    @data({"cells": 40, "smallest_patch_size": 5, "largest_patch_size": 5})
    def test_has_patch_ghost_on_refined_level_case(self, simInput):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_patch_ghost_on_refined_level_case(ndim, True, **simInput)


if __name__ == "__main__":
    unittest.main()
