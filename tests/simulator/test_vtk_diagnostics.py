#

import os
import glob
import h5py
import unittest
import itertools
import numpy as np
from copy import deepcopy
from ddt import data, ddt, unpack

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.core import phare_utilities as phut
from pyphare.simulator.simulator import startMPI
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags

# only in 2d for now
ppc_per_dim = [100, 25, 10]


def config(sim):
    ppc = ppc_per_dim[sim.ndim - 1]

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
    ndims = [2]
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


def permute_multistep(dic, time_step_nbr=3):
    # several dumps, so the per step Steps/ index can be checked against itself
    perms = permute(dic)
    for perm in perms:
        perm["simInput"].update(
            time_step_nbr=time_step_nbr,
            final_time=time_step_nbr * simArgs["final_time"],
        )
    return perms


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

    @data(*permute_multistep({}))
    @unpack
    def test_steps_index(self, ndim, interp, simInput):
        print("test_steps_index dim/interp:{}/{}".format(ndim, interp))

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        local_out = self._run(ndim, interp, simInput)

        if cpp.mpi_rank() != 0:
            return

        files = glob.glob(os.path.join(local_out, "*.vtkhdf"))
        self.assertGreater(len(files), 0)
        for path in sorted(files):
            with h5py.File(path, "r") as h5:
                self._check_steps_index(h5, ndim, os.path.basename(path))

    def _check_steps_index(self, h5, ndim, name):
        # lower dimensions are duplicated to fit a 3d view, cf HierarchyData::X_TIMES
        x_times = [4, 2, 1][ndim - 1]
        steps = h5["VTKHDF/Steps"]
        nsteps = int(steps.attrs["NSteps"])
        self.assertGreater(nsteps, 2)

        for lvl in [key for key in steps if key.startswith("Level")]:
            offsets = {key: steps[lvl][key] for key in ["AMRBoxOffset", "NumberOfAMRBox"]}
            offsets["PointDataOffset/data"] = steps[lvl]["PointDataOffset/data"]

            for key, dataset in offsets.items():
                # int32 saturates at INT32_MAX rather than wrapping, so a dump of more
                # than 2**31 rows silently pins every later step to the same offset
                self.assertGreaterEqual(
                    dataset.dtype.itemsize,
                    8,
                    f"{name} {lvl}/{key} is {dataset.dtype}, too narrow past INT32_MAX",
                )

            boxes = h5["VTKHDF"][lvl]["AMRBox"][:]
            box_offset = offsets["AMRBoxOffset"][:]
            nbr_boxes = offsets["NumberOfAMRBox"][:]
            data_offset = offsets["PointDataOffset/data"][:]
            self.assertEqual(len(data_offset), nsteps)

            # each step starts where the previous one ended, no step may be skipped
            boxes_so_far, rows_so_far = 0, 0
            for step in range(nsteps):
                self.assertEqual(box_offset[step], boxes_so_far, f"{name} {lvl} step {step}")
                self.assertEqual(data_offset[step], rows_so_far, f"{name} {lvl} step {step}")
                for box in boxes[boxes_so_far : boxes_so_far + nbr_boxes[step]]:
                    # AMRBox holds an inclusive cell box, dumped data is all primal
                    nodes = [box[2 * d + 1] - box[2 * d] + 2 for d in range(ndim)]
                    rows_so_far += int(np.prod(nodes)) * x_times
                boxes_so_far += int(nbr_boxes[step])

            self.assertEqual(rows_so_far, h5["VTKHDF"][lvl]["PointData"]["data"].shape[0])
            self.assertEqual(boxes_so_far, boxes.shape[0])


if __name__ == "__main__":
    startMPI()
    unittest.main()
