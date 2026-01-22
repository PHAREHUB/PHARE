#

import unittest
import itertools
import numpy as np
from copy import deepcopy
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


def setup_model(sim, ppc=100):
    L = 0.5

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        return dBy1 + dBy2

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 0.7
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

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
    dic.update(simArgs.copy())
    return [
        dict(
            ndim=ndim,
            interp=interp_order,
            simInput=deepcopy(dic),
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
        self.register_diag_dir_for_cleanup(local_out)
        simInput["diag_options"]["options"]["dir"] = local_out
        simulation = ph.Simulation(**simInput)
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(setup_model(simulation).populations)
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
        atol = 0 if cpp.mpi_size() <= 2 else 1e-16

        for time in [0, simInput["final_time"]]:
            # choose all primal value generally
            phareh5_hier = Run(phareh5_diags).GetVi(time)
            vtk_hier = Run(vtk_diags).GetVi(time)

            eqr = hootils.hierarchy_compare(phareh5_hier, vtk_hier, atol)
            if not eqr:
                print(eqr)
            self.assertTrue(eqr)

    @data(*permute({}))
    @unpack
    def test_missing_level_case(self, ndim, interp, simInput):
        simInput.update(
            dict(
                refinement="tagging",
                max_nbr_levels=2,
                tagging_threshold=0.99,  # prevent level,
            )
        )
        simInput["restart_options"] = dict(
            dir="phare_outputs/vtk_padding_test", mode="overwrite", timestamps=[0.001]
        )
        vtk_diags = self._run(ndim, interp, simInput, "test_vtk")

        simInput.update(
            dict(
                final_time=0.002,
                tagging_threshold=0.01,  # prefer level
            )
        )

        simInput["restart_options"]["timestamps"] = []
        simInput["restart_options"]["restart_time"] = 0.001
        del simInput["diag_options"]["options"]["mode"]  # do not truncate diags
        vtk_diags = self._run(ndim, interp, simInput, "test_vtk")

        hier0 = Run(vtk_diags).GetVi(0)
        hier1 = Run(vtk_diags).GetVi(0.001)
        hier2 = Run(vtk_diags).GetVi(0.002)

        self.assertTrue(1 not in hier0.levels())
        self.assertTrue(1 not in hier1.levels())
        self.assertTrue(1 in hier2.levels())


if __name__ == "__main__":
    startMPI()
    unittest.main()
