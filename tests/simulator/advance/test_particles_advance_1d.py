"""
This file exists independently from test_advance.py to isolate dimension
  test cases and allow each to be overridden in some way if required.
"""

import unittest
import itertools
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.core.box import Box1D
from pyphare.cpp import supported_particle_layouts

from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest

ph.NO_GUI()

ndim = 1
interp_orders = [1, 2, 3]


def permute(boxes={}):
    return [
        dict(
            interp_order=interp_order,
            refinement_boxes=boxes,
            sim_setup_kwargs=dict(layout=layout),
        )
        for interp_order, layout in itertools.product(
            interp_orders, supported_particle_layouts()
        )
    ]


@ddt
class AdvanceTest1D(HybridAdvanceTest):
    @data(*permute())
    @unpack
    def test_L0_particle_number_conservation(self, interp_order, **kwargs):
        self._test_L0_particle_number_conservation(ndim, interp_order, **kwargs)

    @data(
        *permute({"L0": {"B0": Box1D(10, 14)}}),
    )
    @unpack
    def test_domain_particles_on_refined_level(self, interp_order, **kwargs):
        self._test_domain_particles_on_refined_level(ndim, interp_order, **kwargs)


if __name__ == "__main__":
    unittest.main()
