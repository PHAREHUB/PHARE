#!/usr/bin/env python3

import numpy as np
from pathlib import Path
import matplotlib as mpl

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.run import Run
from tests.simulator import SimulatorTest
from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5

diag_dir = "phare_outputs/advance/test_domain_particles_on_refined_level_2/1/1/1/3/"
diag_dir = "phare_outputs/advance/test_field_coarsening_via_subcycles_18/1/1/3/3/"
run = Run(diag_dir)
times = get_times_from_h5(diag_dir + "ions_pop_protons_domain.h5")


pop = "protons"


def p(datahier, time=0):
    for ilvl in datahier.levels():
        print("ilvl", ilvl)
        for patch in datahier.level(ilvl, time).patches:
            print(patch.box)
            for pd_key, pd in patch.patch_datas.items():
                if pd_key.endswith("_domain"):
                    print(pd.size())
                    for icell in pd.dataset.iCells:
                        assert icell in patch.box


for time in times:
    print("time", time)
    p(run.GetParticles(time, pop), time)
