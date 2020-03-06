#!/usr/bin/env python3
#
# formatted with black

import unittest, os, phare.pharein as ph
from datetime import datetime, timezone
from ddt import ddt, data
from tests.simulator import test_simulator as tst
from tests.simulator.py import basicSimulatorArgs, makeBasicModel

out = "phare_outputs/valid/refinement_boxes/"
diags = {"diag_options": {"format": "phareh5", "options": {"dir": out}}}


@ddt
class SimulatorRefineBoxInputsB(unittest.TestCase):
    def dup(dic):
        dic.update(diags.copy())
        return dic

    valid1D = [
        dup({"refinement_boxes": {"L0": {"B0": [(10,), (50,)]}}}),
        dup({"refinement_boxes": {"L0": {"B0": [(5,), (55,)]}}}),
    ]

    invalid1D = [
        dup({"max_nbr_levels": 1, "refinement_boxes": {"L0": {"B0": [(10,), (50,)]}}}),
        dup({"cells": 55, "refinement_boxes": {"L0": {"B0": [(5,), (65,)]}}}),
        dup({"smallest_patch_size": 100, "largest_patch_size": 64,}),
    ]


    def _create_simulator(self, dim, interp, **input):
        ph.globals.sim = None
        ph.Simulation(**basicSimulatorArgs(dim, interp, **input))
        makeBasicModel()
        ph.populateDict()
        hier = tst.make_hierarchy()
        sim = tst.make_simulator(hier)
        sim.initialize()
        return [hier, sim, tst.make_diagnostic_manager(sim, hier)]


    def _do_dim(self, dim, input, valid: bool = False):
        for interp in range(1, 4):
            try:
                hier, sim, dman = self._create_simulator(dim, interp, **input)
                self.assertTrue(valid)
                dman.dump()
                [tst.unmake(x) for x in [hier, sim, dman]]
                tst.reset()
            except ValueError as e:
                self.assertTrue(not valid)

    @data(*valid1D)
    def test_1d_valid(self, input):
        self._do_dim(1, input, True)

    @data(*invalid1D)
    def test_1d_invalid(self, input):
        self._do_dim(1, input)


if __name__ == "__main__":
    unittest.main()
