#!/usr/bin/env python3
#
# formatted with black

import unittest, os, phare.pharein as ph
from datetime import datetime, timezone
from ddt import ddt, data
from phare import populateDict
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

    def setUp(self):
        ph.globals.sim = None

    def tearDown(self):
        tst.unmake_simulator()

    def _do_refinement_boxes(self, dim, interp, **rbs):
        ph.Simulation(**basicSimulatorArgs(dim, interp, **rbs))
        makeBasicModel()
        populateDict()
        tst.make_simulator().initialize()
        tst.make_diagnostic_manager().dump()

    def _do_dim(self, dim, rbs, valid: bool = False):
        dd = out + str(dim) + "/"
        for interp in range(1, 4):
            rbs["diag_options"]["options"]["dir"] = dd + str(interp)
            try:
                self.setUp()
                self._do_refinement_boxes(dim, interp, **rbs)
                self.assertTrue(valid)
            except ValueError as e:
                self.assertTrue(not valid)
            self.tearDown()

    @data(*valid1D)
    def test_1d_valid(self, rbs):
        self._do_dim(1, rbs, True)

    @data(*invalid1D)
    def test_1d_invalid(self, rbs):
        self._do_dim(1, rbs)


if __name__ == "__main__":
    unittest.main()
