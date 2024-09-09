import copy

import time
import h5py
import datetime
import unittest
import numpy as np


import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags
from pyphare.pharesee.hierarchy.patchdata import ParticleData
from pyphare.pharesee.hierarchy.fromh5 import get_all_available_quantities_from_h5
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare


def setup_model(sim, ppc=100):
    model = ph.MaxwellianFluidModel(
        protons={"mass": 1, "charge": 1, "nbr_part_per_cell": ppc},
        alpha={"mass": 4.0, "charge": 1, "nbr_part_per_cell": ppc},
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    dump_all_diags(model.populations)
    return model


timestep = 0.001
out = "phare_outputs/restarts/test/test_restarts_1/1/1/1/00000.00400"
diags = "phare_outputs/restarts"

simInitArgs = dict(
    largest_patch_size=100,
    time_step_nbr=2,
    time_step=timestep,
    cells=200,
    dl=0.3,
    init_options=dict(dir=out, mode="overwrite"),
    diag_options=dict(format="phareh5", options=dict(dir=diags, mode="overwrite")),
)


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
        # h5_filepath = "phare_outputs/restarts/test/test_restarts_1/1/1/1/00000.00400/restore.000000/nodes.0000001/proc.0000000"
        # h5 = H5Content(h5_filepath)
        # print(
        #     "h5.file[]",
        #     h5.file[
        #         "/PHARE_hierarchy/level_0000/level_0000-patch_0000000-block_0000000/EM_B_y##default/d_box"
        #     ][:],
        # )

        sim = ph.Simulation(**copy.deepcopy(simInitArgs))
        model = setup_model(sim)
        Simulator(sim).initialize().reset()

        datahier0 = get_all_available_quantities_from_h5(diags)
        datahier1 = get_all_available_quantities_from_h5(out)

        eq = hierarchy_compare(datahier0, datahier1)
        if not eq:
            print(eq)
        assert eq


if __name__ == "__main__":
    unittest.main()
