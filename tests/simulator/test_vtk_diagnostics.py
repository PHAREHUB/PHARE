import itertools
import unittest
from copy import deepcopy

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core import phare_utilities as phut
from pyphare.simulator.simulator import Simulator, startMPI
from tests.diagnostic import dump_all_diags
from tests.simulator import SimulatorTest

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
            "charge": 2,
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
    "cells": 40,
    "dl": 0.3,
    "diag_options": {
        "format": "pharevtkhdf",
        "options": {"mode": "overwrite", "dir": "phare_outputs/vtk_diagnostic_test"},
    },
}


_ref_layout = "AoSMapped"


def permute(dic, layouts=(_ref_layout,)):
    ndims = [2]
    interp_orders = [1]
    dic.update(simArgs.copy())
    return [
        dict(
            ndim=ndim,
            interp=interp_order,
            layout=layout,
            simInput=deepcopy(dic),
        )
        for ndim, interp_order, layout in itertools.product(
            ndims, interp_orders, layouts
        )
    ]


def cmp_layouts():
    return [l for l in cpp.supported_model_layouts() if l != _ref_layout]


def h5_datasets(h5_file):
    datasets = {}
    h5_file.visititems(
        lambda name, obj: (
            datasets.update({name: obj[()]}) if hasattr(obj, "dtype") else None
        )
    )
    return datasets


@ddt
class VTKDiagnosticsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def _run(self, ndim, interp, simInput, layout=_ref_layout, **kwargs):
        simInput = deepcopy(simInput)
        for key in ["cells", "dl"]:
            simInput[key] = list(phut.np_array_ify(simInput[key], ndim))
        simulation = config(
            self.simulation(interp_order=interp, model_layout=layout, **simInput)
        )
        self.assertTrue(len(simulation.cells) == ndim)
        dump_all_diags(simulation.model.populations)
        Simulator(simulation).run().reset()
        return simulation.diag_options["options"]["dir"]

    @data(*permute({}))
    @unpack
    def test_dump_diags(self, ndim, interp, layout, simInput):
        print(f"test_dump_diags dim/interp:{ndim}/{interp}")

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        local_out = self._run(ndim, interp, simInput, layout)

        if cpp.mpi_rank() == 0:
            import h5py

            for pop, mass, charge in [("protons", 1, 1), ("alpha", 4, 2)]:
                for qty in ["density", "charge_density", "flux"]:
                    file = f"{local_out}/ions_pop_{pop}_{qty}.vtkhdf"
                    with h5py.File(file, "r") as h5_file:
                        self.assertEqual(h5_file.attrs["pop_mass"], mass)
                        self.assertEqual(h5_file.attrs["pop_charge"], charge)

            try:
                from pyphare.pharesee.phare_vtk import plot as plot_vtk

                plot_vtk(local_out + "/EM_B.vtkhdf", f"B{ndim}d_interp{interp}.vtk.png")
                plot_vtk(local_out + "/EM_E.vtkhdf", f"E{ndim}d_interp{interp}.vtk.png")
            except ModuleNotFoundError:
                print("WARNING: vtk python module not found - cannot make plots")

    @data(*permute({}, cmp_layouts()))
    @unpack
    def test_layout_matches_reference(self, ndim, interp, layout, simInput):
        """tiled fields are reduced then converted to primal, as AoSMapped's are converted"""
        print(
            f"test_layout_matches_reference dim/interp/layout:{ndim}/{interp}/{layout}"
        )

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}
        ref_out = self._run(ndim, interp, simInput)
        cmp_out = self._run(ndim, interp, simInput, layout)

        if cpp.mpi_rank() == 0:
            import glob
            import os

            import h5py

            files = sorted(
                os.path.basename(f) for f in glob.glob(f"{ref_out}/*.vtkhdf")
            )
            self.assertTrue(len(files) > 0)
            for file in files:
                with (
                    h5py.File(f"{ref_out}/{file}", "r") as ref,
                    h5py.File(f"{cmp_out}/{file}", "r") as cmp,
                ):
                    ref_ds, cmp_ds = h5_datasets(ref), h5_datasets(cmp)
                    self.assertEqual(ref_ds.keys(), cmp_ds.keys(), file)
                    for name, ref_data in ref_ds.items():
                        np.testing.assert_allclose(
                            cmp_ds[name],
                            ref_data,
                            atol=1e-13,
                            rtol=0,
                            err_msg=f"{file}:{name}",
                        )


if __name__ == "__main__":
    startMPI()
    unittest.main()
