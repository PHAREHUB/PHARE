#!/usr/bin/env python3


import os
import copy
import unittest
import itertools
import numpy as np
from time import sleep
from pathlib import Path
from copy import deepcopy
from ddt import data, ddt, unpack

import pyphare.pharein as ph

from pyphare.pharein.simulation import supported_dimensions
from pyphare.pharesee.hierarchy.fromh5 import h5_filename_from, h5_time_grp_key
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags


ppc_per_dim = [100, 25, 10]


def config(sim):
    L = sim.simulation_domain()
    ppc = ppc_per_dim[sim.ndim - 1]

    def density(*xyz):
        return 1.0

    def by(*xyz):
        return np.asarray(
            [0.1 * np.sin(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
        ).prod(axis=0)

    def bz(*xyz):
        return np.asarray(
            [0.1 * np.sin(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
        ).prod(axis=0)

    def bx(*xyz):
        return 1.0

    def vx(*xyz):
        return 0.0

    def vy(*xyz):
        return np.asarray(
            [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
        ).prod(axis=0)

    def vz(*xyz):
        return np.asarray(
            [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
        ).prod(axis=0)

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
    interp_orders = [1, 2, 3]
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


def _h5_time_group_count(h5_filepath):
    import h5py  # see doc/conventions.md section 2.1.1

    with h5py.File(h5_filepath, "r") as h5_file:
        return len(h5_file[h5_time_grp_key].keys())


# Dedicated 1d setup for the dump-cadence tests below. These tests drive adaptive dt,
# which only the MHD solver supports (SolverPPC::computeStableDt is blocked in hybrid),
# so they run an MHD-only sim: a uniform, quiet medium on a handful of cells. They only
# inspect *when* dumps happen, never their contents, so nothing more is needed.
cadence_out = "phare_outputs/adaptive_cadence"


def cadence_config(sim):
    def density(x):
        return 1.0

    def bx(x):
        return 1.0

    def by(x):
        return 0.0

    def bz(x):
        return 0.0

    def v(x):
        return 0.0

    def p(x):
        return 1.0

    ph.MHDModel(density=density, vx=v, vy=v, vz=v, bx=bx, by=by, bz=bz, p=p)
    return sim


def cadence_args(diagdir, **extra):
    args = dict(
        interp_order=None,  # MHD-only: no hybrid model, hence no interpolation order
        model_options=["MHDModel"],
        reconstruction="Linear",
        limiter="VanLeer",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
        hall=False,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        max_mhd_level=1,
        time_step={"mode": "adaptive", "cfl_wave": 0.8},
        final_time=100.0,  # generous: we control the number of steps manually
        boundary_types="periodic",
        cells=20,
        dl=0.3,
        diag_options=dict(
            format="phareh5", options=dict(dir=diagdir, mode="overwrite")
        ),
    )
    args.update(extra)
    return args


@ddt
class DiagnosticsTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(DiagnosticsTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super().tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def _check_diags(self, sim, times):
        import h5py  # see doc/conventions.md section 2.1.1

        diag_path = sim.diag_options["options"]["dir"]
        py_attrs = [f"{dep}_version" for dep in ["samrai", "highfive", "pybind"]]
        py_attrs += ["git_hash", "serialized_simulation"]
        particle_files = 0
        for diagname, diagInfo in sim.diagnostics.items():
            h5_filepath = os.path.join(diag_path, h5_filename_from(diagInfo))
            self.assertTrue(os.path.exists(h5_filepath))

            self.assertTrue(Path(h5_filepath).exists())
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
    def test_dump_diags(self, ndim, interp, simInput):
        print("test_dump_diags ndim/interp:{}/{}".format(ndim, interp))

        # configure simulation ndim sized values
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(ndim)]

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}

        sim = config(self.simulation(interp_order=interp, **simInput))
        self.assertTrue(len(sim.cells) == ndim)

        dump_all_diags(sim.model.populations)
        self.simulator = Simulator(sim).initialize().advance().reset()

        self.assertTrue(
            any(
                [
                    diagInfo.quantity.endswith("domain")
                    for diagname, diagInfo in ph.global_vars.sim.diagnostics.items()
                ]
            )
        )

        self._check_diags(sim, ["0.0000000000", "0.0010000000"])

    def test_dump_elapsed_time_diags(self, ndim=1, interp=1):
        print("test_dump_elapsed_time_diags dim/interp:{}/{}".format(ndim, interp))

        simInput = copy.deepcopy(simArgs)
        # configure simulation ndim sized values
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(ndim)]

        b0 = [[10 for i in range(ndim)], [19 for i in range(ndim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}

        del simInput["diag_options"]["options"]["fine_dump_lvl_max"]  # don't want

        sim = config(self.simulation(interp_order=interp, **simInput))
        self.assertTrue(len(sim.cells) == ndim)

        dump_all_diags(sim.model.populations)
        for diagname, diagInfo in sim.diagnostics.items():
            diagInfo.write_timestamps = []  # disable
            diagInfo.elapsed_timestamps = [0]  # expect init dump
        simulator = Simulator(sim).setup()
        sleep(3)  # wait so time "elapses"
        simulator.initialize().reset()

        self._check_diags(sim, ["0.0000000000"])

    # ---- dump cadence scheduling --------------
    #
    # Check that
    #  - dumping cadence does not stay stuck when several
    #    scheduled times fall within a single coarse step (dt > gap between timestamps);
    #  - a diagnostic scheduled at final_time still occurs under adaptive dt.

    def test_write_timestamps_catch_up_does_not_freeze(self):
        """Regression test for the dump-cadence-freeze bug: when several scheduled write
        timestamps fall within a single coarse step, the scheduler must catch up (advance past
        all of them) instead of consuming one and freezing."""
        n_advances = 4
        # a timestamp grid far finer than any plausible dt forces the catch-up loop to consume
        # several scheduled times per step. final_time stays large enough that n_advances steps
        # don't reach it. Fed directly as write_timestamps (no period option involved).
        sim = cadence_config(
            self.simulation(
                **cadence_args(cadence_out + "/timestamps_catchup", final_time=1.0)
            )
        )
        ph.ElectromagDiagnostics(
            quantity="B", write_timestamps=np.arange(0.0, 1.0, 1e-4)
        )

        self.simulator = Simulator(sim).initialize()
        for _ in range(n_advances):
            self.simulator.advance()
        self.simulator.reset()

        diag_dir = sim.diag_options["options"]["dir"]
        diagInfo = next(iter(sim.diagnostics.values()))
        h5_filepath = os.path.join(diag_dir, h5_filename_from(diagInfo))
        # one dump per step (init + each advance): if the cadence had frozen after the
        # first step, this would be 1 instead of n_advances + 1
        self.assertEqual(_h5_time_group_count(h5_filepath), n_advances + 1)

    def test_dump_scheduled_exactly_at_final_time_is_not_dropped(self):
        """Regression test: under adaptive dt, timeStep() clamps to 0 on the very last step
        (currentTime()==endTime()). If dump() re-derived its timestep from a fresh timeStep()
        query instead of the dt that actually produced the current state, a diagnostic scheduled
        exactly at final_time would silently never fire (0 < 0 is always false)."""
        final_time = 0.2  # small: run() drives the whole thing to completion, few steps expected

        sim = cadence_config(
            self.simulation(
                **cadence_args(cadence_out + "/final_time_dump", final_time=final_time)
            )
        )
        ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([final_time]))

        self.simulator = Simulator(sim)
        self.simulator.run()
        self.simulator = None  # run() already reset()

        diag_dir = sim.diag_options["options"]["dir"]
        diagInfo = next(iter(sim.diagnostics.values()))
        h5_filepath = os.path.join(diag_dir, h5_filename_from(diagInfo))
        self.assertEqual(_h5_time_group_count(h5_filepath), 1)


if __name__ == "__main__":
    startMPI()
    unittest.main()
