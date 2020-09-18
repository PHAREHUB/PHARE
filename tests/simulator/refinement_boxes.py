#!/usr/bin/env python3
#
# formatted with black

from pybindlibs import cpp

import unittest, os, pyphare.pharein as ph
from datetime import datetime, timezone
from ddt import ddt, data
from tests.diagnostic import dump_all_diags
from tests.simulator import NoOverwriteDict, populate_simulation
from pyphare.simulator.simulator import Simulator,startMPI
import shutil

out = "phare_outputs/valid/refinement_boxes/"
diags = {"diag_options": {"format": "phareh5", "options": {"dir": out, "mode":"overwrite" }}}


@ddt
class SimulatorRefineBoxInputs(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(SimulatorRefineBoxInputs, self).__init__(*args, **kwargs)
        self.simulator = None
        # so we can delete previous diags only on mpi_rank 0
        startMPI()


    def dup(dic):
        dic = NoOverwriteDict(dic)
        dic.update(diags.copy())
        dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
        return dic
    """
      The first set of boxes "B0": [(10,), (14,)]
      Are configured to force there to be a single patch on L0
      This creates a case with MPI that there are an unequal number of
      Patches across MPI domains. This case must be handled and not hang due
      to collective calls not being handled properly.
    """
    valid1D = [
        dup({"refinement_boxes": {"L0": {"B0": [(10,), (14,)]}}}),
        dup({"refinement_boxes": {"L0": {"B0": [(5,), (55,)]}}}),
    ]

    invalid1D = [
        dup({"max_nbr_levels": 1, "refinement_boxes": {"L0": {"B0": [(10,), (50,)]}}}),
        dup({"cells": 55, "refinement_boxes": {"L0": {"B0": [(5,), (65,)]}}}),
        dup({"smallest_patch_size": 100, "largest_patch_size": 64,}),
        dup({"refined_particle_nbr": 1}),
    ]

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()



    def _do_dim(self, dim, input, valid: bool = False):
        for interp in range(1, 4):

            try:
                print("START {}".format(cpp.mpi_rank()))
                self.simulator = Simulator(populate_simulation(dim, interp, **input))
                self.simulator.initialize()

                self.assertTrue(valid)

                self.simulator.diagnostics().dump(self.simulator.currentTime(), self.simulator.timeStep())

                self.simulator = None

                # delete previous diags / can't truncate
                if cpp.mpi_rank() == 0 and os.path.exists("phare_outputs"):
                    shutil.rmtree("phare_outputs")
                    print("RM")
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
