#!/usr/bin/env python3


import os
import copy
import unittest
import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.simulator.simulator import startMPI
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharein.simulation import supported_dimensions
from pyphare.pharesee.hierarchy.fromh5 import h5_filename_from, h5_time_grp_key

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags

cpp = cpp_lib()


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


out = "phare_outputs/diagnostic_test/"
simArgs = {
    "time_step_nbr": 30000,
    "final_time": 30.0,
    "boundary_types": "periodic",
    "cells": 40,
    "dl": 0.3,
    "diag_options": {
        "format": "phareh5",
        "options": {"dir": out, "mode": "overwrite", "fine_dump_lvl_max": 10},
    },
}


def permute(dic):
    args = copy.deepcopy(simArgs)
    args.update(dic)
    dims = supported_dimensions()
    return [[dim, interp, copy.deepcopy(args)] for dim in dims for interp in [1, 2, 3]]


@ddt
class DiagnosticsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(DiagnosticsTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def _check_diags(self, sim, diag_path, times):
        import h5py  # see doc/conventions.md section 2.1.1

        py_attrs = [f"{dep}_version" for dep in ["samrai", "highfive", "pybind"]]
        py_attrs += ["git_hash", "serialized_simulation"]
        particle_files = 0
        for diagname, diagInfo in sim.diagnostics.items():
            h5_filepath = os.path.join(diag_path, h5_filename_from(diagInfo))
            self.assertTrue(os.path.exists(h5_filepath))

            h5_file = h5py.File(h5_filepath, "r")

            self.assertTrue(len(times))
            for time in times:
                self.assertTrue(time in h5_file[h5_time_grp_key])

            h5_py_attrs = h5_file["py_attrs"].attrs.keys()
            for py_attr in py_attrs:
                self.assertIn(py_attr, h5_py_attrs)

            h5_version = h5_file["py_attrs"].attrs["highfive_version"].split(".")
            self.assertTrue(len(h5_version) == 3)
            # semver patch version may contain "-beta" so ignore
            self.assertTrue(all(i.isdigit() for i in h5_version[:2]))

            self.assertTrue(
                ph.simulation.deserialize(
                    h5_file["py_attrs"].attrs["serialized_simulation"]
                ).electrons.closure.Te
                == 0.12
            )

            hier = hierarchy_from(h5_filename=h5_filepath)

            self.assertTrue(hier.sim.electrons.closure.Te == 0.12)

            if h5_filepath.endswith("domain.h5"):
                particle_files += 1
                self.assertTrue("pop_mass" in h5_file.attrs)

                if "protons" in h5_filepath:
                    self.assertTrue(h5_file.attrs["pop_mass"] == 1)
                elif "alpha" in h5_filepath:
                    self.assertTrue(h5_file.attrs["pop_mass"] == 4)
                else:
                    raise RuntimeError("Unknown population")

                self.assertGreater(len(hier.level(0).patches), 0)

                for patch in hier.level(0).patches:
                    self.assertTrue(len(patch.patch_datas.items()))
                    for qty_name, pd in patch.patch_datas.items():
                        splits = pd.dataset.split(ph.global_vars.sim)
                        self.assertTrue(splits.size() > 0)
                        self.assertTrue(pd.dataset.size() > 0)
                        self.assertTrue(
                            splits.size()
                            == pd.dataset.size() * sim.refined_particle_nbr
                        )

        self.assertEqual(particle_files, ph.global_vars.sim.model.nbr_populations())

    @data(
        *permute({"smallest_patch_size": 10, "largest_patch_size": 20}),
        *permute({"smallest_patch_size": 20, "largest_patch_size": 20}),
        *permute({"smallest_patch_size": 20, "largest_patch_size": 40}),
    )
    @unpack
    def test_dump_diags(self, dim, interp, simInput):
        print("test_dump_diags dim/interp:{}/{}".format(dim, interp))

        # configure simulation dim sized values
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(dim)]

        b0 = [[10 for i in range(dim)], [19 for i in range(dim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}

        diag_path = self.unique_diag_dir_for_test_case(f"{out}/test", dim, interp)
        simInput["diag_options"]["options"]["dir"] = diag_path

        simulation = ph.Simulation(**simInput)
        self.register_diag_dir_for_cleanup(diag_path)
        self.assertTrue(len(simulation.cells) == dim)

        dump_all_diags(setup_model().populations)
        self.simulator = Simulator(simulation).initialize().advance().reset()

        self.assertTrue(
            any(
                [
                    diagInfo.quantity.endswith("domain")
                    for diagname, diagInfo in ph.global_vars.sim.diagnostics.items()
                ]
            )
        )
        self._check_diags(simulation, diag_path, ["0.0000000000", "0.0010000000"])

    def test_dump_elapsed_time_diags(self, dim=1, interp=1):
        print("test_dump_elapsed_time_diags dim/interp:{}/{}".format(dim, interp))

        simInput = copy.deepcopy(simArgs)
        # configure simulation dim sized values
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(dim)]

        b0 = [[10 for i in range(dim)], [19 for i in range(dim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}

        diag_path = self.unique_diag_dir_for_test_case(f"{out}/test", dim, interp)
        simInput["diag_options"]["options"]["dir"] = diag_path

        simulation = ph.Simulation(**simInput)
        self.register_diag_dir_for_cleanup(diag_path)
        self.assertTrue(len(simulation.cells) == dim)

        dump_all_diags(setup_model().populations)
        for diagname, diagInfo in simulation.diagnostics.items():
            diagInfo.write_timestamps = []  # disable
            diagInfo.elapsed_timestamps = [0]  # expect init dump
        self.simulator = Simulator(simulation).initialize().reset()

        self._check_diags(simulation, diag_path, ["0.0000000000"])


if __name__ == "__main__":
    startMPI()
    unittest.main()
