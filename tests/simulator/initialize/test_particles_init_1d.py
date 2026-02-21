"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest

import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box1D

from tests.simulator.test_initialization import InitializationTest

matplotlib.use("Agg")  # for systems without GUI

ndim = 1
interp_orders = [1, 2, 3]
ppc = 25


def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]


@ddt
class Initialization1DTest(InitializationTest):
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
                ndim,
                interp_order,
                refinement_boxes,
                "particles",
                cells=30,
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

    @data("berger", "tile")
    def test_amr_clustering(self, clustering):
        dim = 1
        interp_order = 1
        self.getHierarchy(
            dim,
            interp_order,
            {"L0": {"B0": [(10,), (20,)]}},
            "particles",
            clustering=clustering,
        )


if __name__ == "__main__":
    unittest.main()
