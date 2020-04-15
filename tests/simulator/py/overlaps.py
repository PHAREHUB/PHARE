#!/usr/bin/env python3
#
# formatted with black

import os
from phare import cpp

from tests.diagnostic import dump_all_diags
from tests.simulator.py import InitValueValidation

from phare.pp.diagnostics import Diagnostics

from phare.pp.diagnostic.data.overlap import Overlap
from phare.pp.diagnostic.data.em_overlap import EMOverlap, getEMOverlapsFrom
from phare.pp.diagnostic.data.fluid_overlap import FluidOverlap, getFluidOverlapsFrom
from phare.pp.diagnostic.data.particle_overlap import (
    ParticleOverlap,
    ParticleOverlapComparator,
    getParticleOverlapsFrom,
)
from phare.pp.diagnostic.data.level_ghost_particle_overlap import (
    LevelGhostParticleOverlap,
)
from phare.pp.diagnostic.data.patch_ghost_particle_overlap import (
    PatchGhostParticleOverlap,
)

from datetime import datetime, timezone
from ddt import ddt, data


diag_out_dir = "py/overlaps/"


class SimulatorOverlaps(InitValueValidation):
    min_interp = 1
    max_interp = 3

    nbr_overlaps = {
        # 12 intralevel on two levels, 6 periodic on lvl 0, 6 = xyz for E/B
        "em": 18,
        # 2 pops, 1 vecfields + 1 field per pop , + ion density/bulkV
        "fluid": 36,
        # 10 =  2 pops * (2 intralevel, 1 periodic patchghost) + 2 levelghost
        "particles": 10,
    }

    def test_1d_valid(self):
        min_interp = SimulatorOverlaps.min_interp
        max_interp = SimulatorOverlaps.max_interp + 1
        dim = 1
        for interp in range(min_interp, max_interp):
            out = diag_out_dir + str(dim) + "_" + str(interp)
            diagsList = [Diagnostics(out + "_" + str(i)).diags for i in range(0, 2)]
            self._checkEM(diagsList, overlapNbr=SimulatorOverlaps.nbr_overlaps["em"])
            self._checkFluid(
                diagsList, overlapNbr=SimulatorOverlaps.nbr_overlaps["fluid"]
            )
            self._checkParticles(
                diagsList, overlapNbr=SimulatorOverlaps.nbr_overlaps["particles"]
            )

    def _compare(self, overlaps0, overlaps1, overlapNbr):
        self.assertTrue(len(overlaps0) == len(overlaps1) == overlapNbr)

    def _checkFluid(self, diagsList, overlapNbr):
        self._compare(*[getFluidOverlapsFrom(diags) for diags in diagsList], overlapNbr)

    def _checkEM(self, diagsList, overlapNbr):
        self._compare(*[getEMOverlapsFrom(diags) for diags in diagsList], overlapNbr)

    def _checkParticles(self, diagsList, overlapNbr):
        self._compare(
            *[getParticleOverlapsFrom(diags) for diags in diagsList], overlapNbr
        )

    ## uncomment to generate files
    # def setUp(self):
    #     dim = 1
    #     min_interp = SimulatorOverlaps.min_interp
    #     max_interp = SimulatorOverlaps.max_interp + 1
    #     for interp in range(min_interp, max_interp):
    #         self._generate_hdf5_files(dim, interp)

    # def _generate_hdf5_files(self, dim, interp):
    #     """generate hdf5 files"""
    #     out = diag_out_dir + str(dim) + "_" + str(interp)
    #     sim0 = SimulatorOverlaps.valid1D[0]
    #     sim0["diag_options"]["options"]["dir"] = out + "_0"
    #     sim1 = SimulatorOverlaps.valid1D[1]
    #     sim1["diag_options"]["options"]["dir"] = out + "_1"
    #     sim1["origin"] = [0.25 for d in range(dim)]
    #     [self.runAndDump(dim=1, interp=interp, input=simput) for simput in [sim0, sim1]]

    # def add_to_dict(dic):
    #     dic.update(InitValueValidation.diag_options(diag_out_dir))
    #     dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
    #     return dic

    # valid1D = [add_to_dict({}), add_to_dict({})]


if __name__ == "__main__":
    import unittest

    unittest.main()
