import copy

import time
import h5py
import datetime
import unittest
import numpy as np
from pathlib import Path
from datetime import timedelta

from ddt import ddt, data, unpack

from pyphare.cpp import cpp_lib

cpp = cpp_lib()

import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags
from pyphare.pharesee.hierarchy.patchdata import ParticleData
from pyphare.pharesee.hierarchy.fromh5 import get_all_available_quantities_from_h5

# ./build/tests/simulator/phare_outputs/restarts/test/test_restarts_1/1/1/1/00000.00400/restore.000000/nodes.0000001/proc.0000000


def setup_model(ppc=100):
    model = ph.MaxwellianFluidModel(protons={"nbr_part_per_cell": ppc})
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return model


timestep = 0.001
out = "phare_outputs/restarts/test/test_restarts_1/1/1/1/00000.00400"
simArgs = dict(
    time_step_nbr=2,
    time_step=timestep,
    cells=100,
    dl=0.3,
    init_options=dict(dir=out, mode="overwrite"),
)


def dup(dic={}):
    dic.update(copy.deepcopy(simArgs))
    return dic


def traverse_h5_for_groups_recursive(h5content: "H5Content", group, path=""):
    if "level_0000" in path:
        for key in group.attrs:
            h5content.attr[f"{path}/{key}"] = group.attrs[key]
    if isinstance(group, h5py._hl.group.Group):
        for key in group:
            kpath = f"{path}/{key}"
            traverse_h5_for_groups_recursive(h5content, group[key], path=kpath)
    else:
        if "level_0000" not in path:
            return
        h5content.data[path] = group


class H5Content:
    def __init__(self, path):
        self.file = h5py.File(path, "r")
        self.data = {}
        self.attr = {}
        traverse_h5_for_groups_recursive(self, self.file)


@ddt
class RestartsParserTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RestartsParserTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RestartsParserTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_restart_parser(self):
        h5_filepath = "phare_outputs/restarts/test/test_restarts_1/1/1/1/00000.00400/restore.000000/nodes.0000001/proc.0000000"

        h5 = H5Content(h5_filepath)
        print(
            h5.file[
                "/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/EMAvg_B_x##default/d_ghost_box"
            ][:]
        )
        for k in h5.data:
            print(k)

        sim = ph.Simulation(**dup())
        model = setup_model()
        Simulator(sim).initialize()


if __name__ == "__main__":
    unittest.main()
