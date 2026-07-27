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
    # the Steps/ index only breaks *between* steps, a single dump cannot expose it,
    #  so these permutations run several timesteps hence write several dumps
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
    def test_step_offsets_are_64bit_and_contiguous(self, ndim, interp, simInput):
        """A .vtkhdf file is one file for the whole run: every dump appends its rows to
        /VTKHDF/Level<n>/{PointData/data,AMRBox} and records where its own rows start in
        /VTKHDF/Steps/Level<n>/{PointDataOffset/data,AMRBoxOffset,NumberOfAMRBox}. Those
        three datasets are the only way a reader can tell one timestep from the next, so
        a wrong value there does not corrupt the data, it silently serves the wrong step.

        Two ways that index can go wrong, one assertion each:

        1. too narrow. The offsets used to be created as int32 while the values written
           are std::size_t. HDF5 saturates on conversion instead of wrapping, so once a
           run passes 2**31 rows of point data in one dump, every later step stores
           INT32_MAX and readers replay one frozen timestep for the rest of the run. Seen
           on a 4096^2 / 3072 rank run: frozen from step 63 of 101. No test we can afford
           to run reaches 2**31 rows, so the declared width is asserted directly - it is
           the only observable that carries this defect at a size that runs in seconds.

        2. not contiguous. Independently of any type, the offsets must describe exactly
           the rows the writer appended, so they are rebuilt here from the AMRBox
           geometry stored in the same file and compared to what the writer claimed.
        """
        print(
            "test_step_offsets_are_64bit_and_contiguous dim/interp:{}/{}".format(
                ndim, interp
            )
        )

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        local_out = self._run(ndim, interp, simInput)

        if cpp.mpi_rank() != 0:
            return

        files = glob.glob(os.path.join(local_out, "*.vtkhdf"))
        self.assertGreater(len(files), 0)
        for path in sorted(files):
            name = os.path.basename(path)
            with h5py.File(path, "r") as h5:
                steps = h5["VTKHDF/Steps"]
                nsteps = int(steps.attrs["NSteps"])
                self.assertGreater(nsteps, 2)  # nothing to index with fewer

                for lvl in [key for key in steps if key.startswith("Level")]:
                    offsets = {
                        key: steps[lvl][key]
                        for key in [
                            "AMRBoxOffset",
                            "NumberOfAMRBox",
                            "PointDataOffset/data",
                        ]
                    }
                    self._assert_offsets_are_64bit(offsets, name, lvl)
                    self._assert_step_index_is_contiguous(
                        h5, offsets, nsteps, ndim, name, lvl
                    )

    def _assert_offsets_are_64bit(self, offsets, name, lvl):
        for key, dataset in offsets.items():
            self.assertGreaterEqual(
                dataset.dtype.itemsize,
                8,
                f"{name} {lvl}/{key} is {dataset.dtype}, saturates at INT32_MAX",
            )

    def _assert_step_index_is_contiguous(self, h5, offsets, nsteps, ndim, name, lvl):
        # lower dimensions are duplicated to fit a 3d view, cf HierarchyData::X_TIMES
        x_times = [4, 2, 1][ndim - 1]

        boxes = h5["VTKHDF"][lvl]["AMRBox"][:]
        box_offset = offsets["AMRBoxOffset"][:]
        nbr_boxes = offsets["NumberOfAMRBox"][:]
        data_offset = offsets["PointDataOffset/data"][:]
        self.assertEqual(len(data_offset), nsteps)  # one offset per step, no more

        # each step must start exactly where the previous one stopped writing. Exact
        #  equality and not monotonicity: a level missing at a step legitimately repeats
        #  its offset, which is also what a saturated offset looks like, so a weaker
        #  check could not tell the two apart
        boxes_so_far, rows_so_far = 0, 0
        for step in range(nsteps):
            self.assertEqual(
                box_offset[step], boxes_so_far, f"{name} {lvl} step {step}"
            )
            self.assertEqual(
                data_offset[step], rows_so_far, f"{name} {lvl} step {step}"
            )
            for box in boxes[boxes_so_far : boxes_so_far + nbr_boxes[step]]:
                # AMRBox holds an inclusive cell box, dumped data is all primal
                nodes = [box[2 * d + 1] - box[2 * d] + 2 for d in range(ndim)]
                rows_so_far += int(np.prod(nodes)) * x_times
            boxes_so_far += int(nbr_boxes[step])

        # and the last step must close on the arrays it indexes into
        self.assertEqual(rows_so_far, h5["VTKHDF"][lvl]["PointData"]["data"].shape[0])
        self.assertEqual(boxes_so_far, boxes.shape[0])


if __name__ == "__main__":
    startMPI()
    unittest.main()
