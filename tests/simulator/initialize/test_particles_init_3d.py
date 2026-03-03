"""
  This file exists independently from test_initialization.py to isolate dimension
    test cases and allow each to be overridden in some way if required.
"""

import unittest

import matplotlib
from ddt import data, ddt, unpack
from pyphare.core.box import Box3D

from tests.simulator.test_initialization import InitializationTest

matplotlib.use("Agg")  # for systems without GUI

ndim = 3
interp_orders = [1, 2, 3]
ppc, cells = 10, 30


def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]


@ddt
class Initialization3DTest(InitializationTest):
    @data(*interp_orders)
    def test_nbr_particles_per_cell_is_as_provided(self, interp_order):
        print(f"{self._testMethodName}_{ndim}d")
        self._test_nbr_particles_per_cell_is_as_provided(
            ndim, interp_order, ppc, cells=cells
        )

    @data(
        *per_interp(({"L0": {"B0": Box3D(10, 14)}})),
        *per_interp(({"L0": {"B0": Box3D(10, 14)}, "L1": {"B0": Box3D(22, 26)}})),
        *per_interp(({"L0": {"B0": Box3D(2, 6), "B1": Box3D(7, 11)}})),
    )
    @unpack
    def test_levelghostparticles_have_correct_split_from_coarser_particle(
        self, interp_order, refinement_boxes
    ):
        print(f"\n{self._testMethodName}_{ndim}d")
        now = self.datetime_now()
        self._test_levelghostparticles_have_correct_split_from_coarser_particle(
            self.getHierarchy(
                ndim,
                interp_order,
                refinement_boxes,
                "particles",
                cells=cells,
                nbr_part_per_cell=ppc,
            )
        )
        print(
            f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds"
        )

    @data(
        *per_interp(({"L0": {"B0": Box3D(10, 14)}})),
        *per_interp(({"L0": {"B0": Box3D(5, 14)}, "L1": {"B0": Box3D(15, 19)}})),
        *per_interp(({"L0": {"B0": Box3D(2, 12), "B1": Box3D(13, 25)}})),
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
        print(
            f"\n{self._testMethodName}_{ndim}d took {self.datetime_diff(now)} seconds"
        )


if __name__ == "__main__":
    unittest.main()
