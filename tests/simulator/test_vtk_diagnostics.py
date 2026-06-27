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


def config(sim):
    ppc = ppc_per_dim[sim.ndim - 1]

    def density(*xyz):
        return 1.0

    def bx(*xyz):
        return 0.0

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def by(*xyz):
        L = ph.global_vars.sim.simulation_domain()[0]
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

    def vxalpha(*xyz):
        return 3.0

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
    vvvalpha = {
        "vbulkx": vxalpha,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthxyz,
        "vthy": vthxyz,
        "vthz": vthxyz,
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
            "mass": 4.0,
            "charge": 1,
            "density": density,
            **vvvalpha,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 2334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return sim


simArgs = {
    "time_step_nbr": 1,
    "final_time": 0.001,
    "boundary_types": "periodic",
    "cells": 40,
    "dl": 0.3,
    "diag_options": {
        "format": "pharevtkhdf",
        "options": {"mode": "overwrite", "dir": "phare_outputs/vtk_diagnostic_test"},
    },
}


def permute(dic):
    ndims = supported_dimensions()
    interp_orders = [1]
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
        simulation = config(self.simulation(interp_order=interp, **simInput))
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(simulation.model.populations)
        Simulator(simulation).run().reset()
        return simulation.diag_options["options"]["dir"]

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

                plot_vtk(local_out + "/EM_B.vtkhdf", f"B{ndim}d_interp{interp}.vtk.png")
                plot_vtk(local_out + "/EM_E.vtkhdf", f"E{ndim}d_interp{interp}.vtk.png")
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
        simInput["diag_options"]["options"]["dir"] += "/phareh5"
        phareh5_diags = self._run(ndim, interp, simInput, "test_h5")

        # not binary == with more than 2 cores
        atol = 0 if cpp.mpi_size() <= 2 else [1e-15, 1e-14, 1e-14][ndim - 1]
        time = 0
        # choose all primal value generally
        phareh5_hier = Run(phareh5_diags).GetVi(time)
        vtk_hier = Run(vtk_diags).GetVi(time)

        self.assertTrue(phareh5_hier.data_files)
        self.assertTrue(vtk_hier.data_files)
        self.assertTrue(vtk_hier.data_files != phareh5_hier.data_files)
        eqr = hootils.hierarchy_compare(vtk_hier, phareh5_hier, atol)
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
        simInput["restart_options"] = dict(mode="overwrite", timestamps=[0.001])
        vtk_diags = self._run(ndim, interp, simInput, "test_padding_vtk")

        simInput.update(
            dict(
                final_time=0.002,
                tagging_threshold=0.01,  # prefer level
            )
        )

        simInput["restart_options"]["timestamps"] = []
        simInput["restart_options"]["restart_time"] = 0.001
        del simInput["diag_options"]["options"]["mode"]  # do not truncate diags
        vtk_diags = self._run(ndim, interp, simInput, "test_padding_vtk")

        hier0 = Run(vtk_diags).GetVi(0)
        hier1 = Run(vtk_diags).GetVi(0.001)
        hier2 = Run(vtk_diags).GetVi(0.002)

        self.assertTrue(1 not in hier0.levels())
        self.assertTrue(1 not in hier1.levels())
        self.assertTrue(1 in hier2.levels())


def config_3d(**simInput):
    L = 0.5
    ppc = ppc_per_dim[2]

    sim = ph.Simulation(**simInput)
    from pyphare.pharein.diagnostics import all_timestamps

    timestamps = all_timestamps(sim)

    def density(x, y, z):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y, z):
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

    def bx(x, y, z):
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

    def bz(x, y, z):
        return 0.0

    def b2(x, y, z):
        return bx(x, y, z) ** 2 + by(x, y, z) ** 2 + bz(x, y, z) ** 2

    def T(x, y, z):
        K = 0.7
        temp = 1.0 / density(x, y, z) * (K - b2(x, y, z) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y, z):
        return 0.0

    def vy(x, y, z):
        return 0.0

    def vz(x, y, z):
        return 0.0

    def vthx(x, y, z):
        return np.sqrt(T(x, y, z))

    def vthy(x, y, z):
        return np.sqrt(T(x, y, z))

    def vthz(x, y, z):
        return np.sqrt(T(x, y, z))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": ppc,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
            "init": {"seed": cpp.mpi_rank() + 12},
        },
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.LoadBalancer(active=True, mode="nppc", tol=0.05, every=1000)

    pop = "protons"
    # ph.FluidDiagnostics(quantity="bulkVelocity", write_timestamps=timestamps)
    # ph.FluidDiagnostics(
    #     quantity="density", write_timestamps=timestamps, population_name=pop
    # )

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            component_properties={
                "x": [ph.FieldDiagnosticSlice((20, 0, 0), (20, 40, 40))],
                # "y": [ph.FieldDiagnosticSlice((20, 0, 0), (20, 40, 40))],
                # "z": [ph.FieldDiagnosticSlice((20, 0, 0), (20, 40, 40))],
            },
        )

    # ph.InfoDiagnostics(quantity="particle_count", write_timestamps=timestamps)
    print("diags", sim.diagnostics, [d for d in sim.diagnostics])
    return sim


# @ddt
# class VTKDiagnostics3dTest(SimulatorTest):
#     def __init__(self, *args, **kwargs):
#         super(VTKDiagnostics3dTest, self).__init__(*args, **kwargs)
#         self.simulator = None
#         ph.global_vars.sim = None

#     def _run(self, ndim, interp, simInput, diag_dir="", **kwargs):
#         for key in ["cells", "dl", "boundary_types"]:
#             simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
#         local_out = self.unique_diag_dir_for_test_case(
#             f"{out}{'/'+diag_dir if diag_dir else ''}", ndim, interp
#         )
#         self.register_diag_dir_for_cleanup(local_out)
#         simInput["diag_options"]["options"]["dir"] = local_out
#         simulation = config_3d(**simInput)
#         self.assertTrue(len(simulation.cells) == ndim)
#         Simulator(simulation).run().reset()
#         ph.global_vars.sim = None
#         return local_out

#     def test_slice_3d(self):
#         print("test_slice_3d dim/interp:{}/{}".format(3, 1))

#         local_out = self._run(ndim=3, interp=1, simInput=deepcopy(simArgs))

#         if cpp.mpi_rank() == 0:
#             run = Run(local_out)
#             time = 0
#             B = run.GetBSlice(time)
#             B.plot(filename=f"bx_slice.png", qty=f"x")
#             B.magnitude().plot(filename=f"bmag_slice.png")


if __name__ == "__main__":
    startMPI()
    unittest.main()
