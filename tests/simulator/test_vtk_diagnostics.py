#

import unittest
import itertools
import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run

from pyphare.core import phare_utilities as phut
from pyphare.simulator.simulator import startMPI
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_utils as hootils

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


out = "phare_outputs/vtk_diagnostic_test"
simArgs = {
    "time_step_nbr": 10,
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

    def _run(self, ndim, interp, simInput, diag_dir="", **kwargs):
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        local_out = self.unique_diag_dir_for_test_case(
            f"{out}{'/'+diag_dir if diag_dir else ''}", ndim, interp
        )
        simInput["diag_options"]["options"]["dir"] = local_out
        simulation = ph.Simulation(**simInput)
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(setup_model().populations)
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

        try:
            from pyphare.pharesee.phare_vtk import plot as plot_vtk

            plot_vtk(local_out + "/EM_B.vtkhdf", f"B{ndim}d.vtk.png")
            plot_vtk(local_out + "/EM_E.vtkhdf", f"E{ndim}d.vtk.png")
        except ModuleNotFoundError:
            print("WARNING: vtk python module not found - cannot make plots")

    @data(*permute({}))
    @unpack
    def test_compare_l0_primal_to_phareh5(self, ndim, interp, simInput):
        print("test_compare_l0_primal_to_phareh5 dim/interp:{}/{}".format(ndim, interp))

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}

        vtk_diags = self._run(ndim, interp, simInput, "test_vtk")

        simInput["diag_options"]["format"] = "phareh5"
        phareh5_diags = self._run(ndim, interp, simInput, "test_h5")

        # not binary == with more than 2 cores
        atol = 0 if cpp.mpi_size() <= 2 else 1e-17

        for time in [0, simInput["final_time"]]:
            # choose all primal value generally
            phareh5_hier = Run(phareh5_diags).GetVi(time)
            vtk_hier = Run(vtk_diags).GetVi(time)

            eqr = hootils.hierarchy_compare(phareh5_hier, vtk_hier, atol)
            if not eqr:
                print(eqr)
            self.assertTrue(eqr)


if __name__ == "__main__":
    startMPI()
    unittest.main()
