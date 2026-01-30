#

import unittest
import itertools
import numpy as np
from copy import deepcopy
from ddt import data, ddt, unpack

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run

from pyphare.core import phare_utilities as phut
from pyphare.simulator.simulator import startMPI
from pyphare.simulator.simulator import Simulator
from pyphare.pharein.simulation import supported_dimensions
from pyphare.pharesee.hierarchy import hierarchy_utils as hootils

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags

ppc_per_dim = [100, 25, 10]


def config(ndim, interp, **simInput):
    ppc = ppc_per_dim[ndim - 1]
    sim = ph.Simulation(**simInput)

    def density(*xyz):
        return 1.0

    def bx(*xyz):
        return 0.0

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def by(*xyz):
        L = sim.simulation_domain()[0]
        v1, v2 = -1, 1.0
        return v1 + (v2 - v1) * (S(xyz[0], L * 0.25, 1) - S(xyz[0], L * 0.75, 1))

    def bz(*xyz):
        return 0.5

    def b2(*xyz):
        return bx(xyz[0]) ** 2 + by(xyz[0]) ** 2 + bz(xyz[0]) ** 2

    def T(*xyz):
        K = 1
        return 1 / density(xyz[0]) * (K - b2(xyz[0]) * 0.5)

    def vx(*xyz):
        return 2.0

    def vy(*xyz):
        return 0.0

    def vz(*xyz):
        return 0.0

    def vthxyz(*xyz):
        return T(xyz[0])

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
    }

    ph.MaxwellianFluidModel(
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
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return sim


out = "phare_outputs/vtk_diagnostic_test"
simArgs = {
    "time_step_nbr": 1,
    "final_time": 0.001,
    "boundary_types": "periodic",
    "cells": 40,
    "dl": 0.3,
    "diag_options": {
        "format": "pharevtkhdf",
        "options": {"dir": out, "mode": "overwrite"},
    },
}


def permute(dic):
    interp_orders = [1]  # no ghosts on vtk data so not so relevant
    dic.update(simArgs.copy())
    return [
        dict(
            ndim=ndim,
            interp=interp_order,
            simInput=deepcopy(dic),
        )
        for ndim, interp_order in itertools.product(
            supported_dimensions(), interp_orders
        )
    ]


@ddt
class VTKDiagnosticsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(VTKDiagnosticsTest, self).__init__(*args, **kwargs)
        self.simulator = None
        ph.global_vars.sim = None

    def _run(self, ndim, interp, simInput, diag_dir="", **kwargs):
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        local_out = self.unique_diag_dir_for_test_case(
            f"{out}{'/'+diag_dir if diag_dir else ''}", ndim, interp
        )
        self.register_diag_dir_for_cleanup(local_out)
        simInput["diag_options"]["options"]["dir"] = local_out
        simulation = config(ndim, interp, **simInput)
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(simulation.model.populations)
        Simulator(simulation).run().reset()
        ph.global_vars.sim = None
        return local_out

    @data(*permute({}))
    @unpack
    def test_dump_diags(self, ndim, interp, simInput):
        print("test_dump_diags dim/interp:{}/{}".format(ndim, interp))

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        local_out = self._run(ndim, interp, simInput)

        if cpp.mpi_rank() == 0:
            try:
                from pyphare.pharesee.phare_vtk import plot as plot_vtk

                plot_vtk(local_out + "/EM_B.vtkhdf", f"B{ndim}d.vtk.png")
                plot_vtk(local_out + "/EM_E.vtkhdf", f"E{ndim}d.vtk.png")
            except ModuleNotFoundError:
                print("WARNING: vtk python module not found - cannot make plots")


if __name__ == "__main__":
    startMPI()
    unittest.main()
