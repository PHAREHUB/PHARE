import glob
import itertools
import unittest
from copy import deepcopy

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core import phare_utilities as phut
from pyphare.pharein.simulation import supported_dimensions
from pyphare.pharesee.hierarchy import hierarchy_utils as hootils
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI
from tests.diagnostic import dump_all_diags
from tests.simulator import SimulatorTest

# only in 2d for now
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
            "mass": 4,
            "charge": 2,
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


def assert_vtkhdf_layout(test, path, ndim):
    import h5py

    x_times = {1: 4, 2: 2, 3: 1}[ndim]
    with h5py.File(path, "r") as h5_file:
        root = h5_file["VTKHDF"]
        nsteps = int(root["Steps"].attrs["NSteps"])
        test.assertEqual(nsteps, root["Steps/Values"].shape[0])

        ilvl = 0
        while f"Level{ilvl}/PointData/data" in root:
            data = root[f"Level{ilvl}/PointData/data"][:]
            boxes = root[f"Level{ilvl}/AMRBox"][:]
            steps = root[f"Steps/Level{ilvl}"]
            n_boxes = steps["NumberOfAMRBox"][:]
            box_offsets = steps["AMRBoxOffset"][:]
            data_offsets = steps["PointDataOffset/data"][:]
            data_ends = np.append(data_offsets[1:], data.shape[0])

            for step in range(nsteps):
                step_boxes = boxes[
                    box_offsets[step] : box_offsets[step] + n_boxes[step]
                ]
                test.assertTrue((step_boxes[:, 2 * ndim :] == 0).all())
                lower = step_boxes[:, 0 : 2 * ndim : 2]
                shapes = step_boxes[:, 1 : 2 * ndim : 2] - lower + 2
                values = data[data_offsets[step] : data_ends[step]]
                test.assertEqual(
                    values.shape[0], np.prod(shapes, axis=1).sum() * x_times
                )
                if n_boxes[step] == 0:
                    continue

                origin = lower.min(axis=0)
                extent = (lower + shapes).max(axis=0) - origin
                nodes = np.full(tuple(extent) + (data.shape[1],), np.nan)
                row = 0
                for box_lower, shape in zip(lower, shapes):
                    size = int(np.prod(shape))
                    block = values[row : row + size]
                    for copy in range(1, x_times):
                        np.testing.assert_array_equal(
                            values[row + copy * size : row + (copy + 1) * size], block
                        )
                    row += size * x_times

                    block = block.reshape(tuple(shape[::-1]) + (data.shape[1],))
                    block = np.moveaxis(
                        block, list(range(ndim)), list(range(ndim))[::-1]
                    )
                    region = tuple(
                        slice(lo - o, lo - o + n)
                        for lo, o, n in zip(box_lower, origin, shape)
                    )
                    shared = ~np.isnan(nodes[region])
                    np.testing.assert_allclose(
                        block[shared], nodes[region][shared], atol=1e-5
                    )
                    nodes[region] = block
            ilvl += 1


@ddt
class VTKDiagnosticsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def _run(self, ndim, interp, simInput, diag_dir="", **kwargs):
        for key in ["cells", "dl"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        simulation = config(self.simulation(interp_order=interp, **simInput))
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(simulation.model.populations)
        Simulator(simulation).run().reset()
        return simulation.diag_options["options"]["dir"]

    @data(*permute({}))
    @unpack
    def test_dump_diags(self, ndim, interp, simInput):
        print(f"test_dump_diags dim/interp:{ndim}/{interp}")

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        local_out = self._run(ndim, interp, simInput)

        if cpp.mpi_rank() == 0:
            import h5py

            for pop, mass, charge in [("protons", 1, 1), ("alpha", 4, 2)]:
                for qty in ["density", "charge_density", "flux"]:
                    file = f"{local_out}/ions_pop_{pop}_{qty}.vtkhdf"
                    with h5py.File(file, "r") as h5_file:
                        self.assertEqual(h5_file.attrs["pop_mass"], mass)
                        self.assertEqual(h5_file.attrs["pop_charge"], charge)

            for path in sorted(glob.glob(f"{local_out}/*.vtkhdf")):
                assert_vtkhdf_layout(self, path, ndim)

            try:
                from pyphare.pharesee.phare_vtk import plot as plot_vtk

                plot_vtk(local_out + "/EM_B.vtkhdf", f"B{ndim}d_interp{interp}.vtk.png")
                plot_vtk(local_out + "/EM_E.vtkhdf", f"E{ndim}d_interp{interp}.vtk.png")
            except ModuleNotFoundError:
                print("WARNING: vtk python module not found - cannot make plots")

    def test_level_without_patches_on_some_ranks(self):
        ndim, interp = 2, 1
        simInput = deepcopy(simArgs)
        simInput["refinement_boxes"] = {"L0": {"B0": [[10] * ndim, [12] * ndim]}}
        local_out = self._run(ndim, interp, simInput)

        if cpp.mpi_rank() == 0:
            import h5py

            with h5py.File(f"{local_out}/EM_B.vtkhdf", "r") as h5_file:
                n_boxes = h5_file["VTKHDF/Steps/Level1/NumberOfAMRBox"][:]
            if cpp.mpi_size() > 1:
                self.assertLess(n_boxes.max(), cpp.mpi_size())

            for path in sorted(glob.glob(f"{local_out}/*.vtkhdf")):
                assert_vtkhdf_layout(self, path, ndim)

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


if __name__ == "__main__":
    startMPI()
    unittest.main()
