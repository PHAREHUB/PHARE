#

import unittest
import itertools
import numpy as np

from pathlib import Path
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.simulator.simulator import startMPI
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags

cpp = cpp_lib()

ndims = [2]  # 3d todo / 1d not supported
interp_orders = [1]  # , 2, 3


def setup_model(ppc=100):
    def density(*xyz):
        return 1.0

    def by(*xyz):
        from pyphare.pharein.global_vars import sim

        L = sim.simulation_domain()
        _ = lambda i: 0.1 * np.sin(2 * np.pi * xyz[i] / L[i])
        return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

    def bz(*xyz):
        from pyphare.pharein.global_vars import sim

        L = sim.simulation_domain()
        _ = lambda i: 0.1 * np.sin(2 * np.pi * xyz[i] / L[i])
        return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

    def bx(*xyz):
        return 1.0

    def vx(*xyz):
        return 0.0

    def vy(*xyz):
        from pyphare.pharein.global_vars import sim

        L = sim.simulation_domain()
        _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
        return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

    def vz(*xyz):
        from pyphare.pharein.global_vars import sim

        L = sim.simulation_domain()
        _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
        return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

    def vthx(*xyz):
        return 0.01

    def vthy(*xyz):
        return 0.01

    def vthz(*xyz):
        return 0.01

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
    }

    model = ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "mass": 1,
            "charge": 1,
            "density": density,
            **vvv,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1337},
        },
        alpha={
            "mass": 4,
            "charge": 1,
            "density": density,
            **vvv,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 2334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return model


out = "phare_outputs/vtk_diagnostic_test/"
simArgs = {
    "time_step_nbr": 1,
    "final_time": 0.003,
    "boundary_types": "periodic",
    "cells": 40,
    "dl": 0.3,
    "diag_options": {
        "format": "pharevtkhdf",
        "options": {"dir": out, "mode": "overwrite"},
    },
}


def permute(dic):
    dic.update(simArgs.copy())
    return [
        dict(
            ndim=ndim,
            interp=interp_order,
            simInput=dic,
        )
        for ndim, interp_order in itertools.product(ndims, interp_orders)
    ]


@ddt
class VTKDiagnosticsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(VTKDiagnosticsTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None
        super().tearDown()

    @data(*permute({}))
    @unpack
    def test_dump_diags(self, ndim, interp, simInput):
        print("test_dump_diags dim/interp:{}/{}".format(ndim, interp))
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(ndim)]
        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        local_out = self.unique_diag_dir_for_test_case(f"{out}/test_vtk", ndim, interp)
        simInput["diag_options"]["options"]["dir"] = local_out
        simInput["diag_options"]["format"] = "pharevtkhdf"
        simulation = ph.Simulation(**simInput)
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(setup_model().populations)
        self.simulator = Simulator(simulation).run().reset()

        # maybe use vtk python module to generate artifacts here


if __name__ == "__main__":
    startMPI()
    unittest.main()
