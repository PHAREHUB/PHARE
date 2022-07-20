"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box1D, nDBox
from tests.simulator.test_initialization import InitializationTest

import matplotlib

from pyphare.cpp import cpp_lib
cpp = cpp_lib()

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

    @data({"cells": 40, "interp_order": 1})
    def test_has_patch_ghost_on_refined_level_case(self, simInput):
        print(f"{self._testMethodName}_{ndim}d")
        from pyphare.pharein.simulation import check_patch_size
        diag_outputs=f"phare_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts_{ndim}_{self.ddt_test_id()}"
        _, smallest_patch_size = check_patch_size(ndim, **simInput)
        simInput["smallest_patch_size"] = smallest_patch_size
        simInput["largest_patch_size"] = smallest_patch_size
        self._test_patch_ghost_on_refined_level_case(ndim, True, **simInput)



    @data("berger", "tile")
    def test_amr_clustering(self, clustering):
        interp_order = 1
        test_id = self.ddt_test_id()
        local_out = f"test_amr_clustering/mpi/{cpp.mpi_size()}/{test_id}"
        datahier = self.getHierarchy(interp_order, {"L0": {"B0": [(10, ), (20, )]}}, "particles", clustering=clustering, diag_outputs=local_out)




if __name__ == "__main__":
    unittest.main()
