#!/usr/bin/env python3
import os

import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest

os.environ["PHARE_SCOPE_TIMING"] = "1"  # turn on scope timing

ph.NO_GUI()
cpp = cpp_lib()


def create_settings(cells, dl, m, nbr_periods, timestep, diag_dir):
    settings = {
        "cells": cells,
        "dl": dl,
        "m": m,
        "nbr_periods": nbr_periods,
        "timestep": timestep,
        "diag_dir": diag_dir,
    }

    settings["lx"] = settings["cells"][0] * settings["dl"][0]
    settings["k"] = 2 * np.pi / settings["lx"]
    settings["kt"] = 2 * np.pi / settings["lx"] * settings["m"]
    settings["w"] = (settings["kt"] ** 2 / 2) * (
        np.sqrt(1 + 4 / settings["kt"] ** 2) + 1
    )
    settings["final_time"] = (2 * np.pi / settings["w"]) * settings["nbr_periods"]
    # np.arrange() looked like it had some precision issues
    settings["nsteps"] = int(np.round(settings["final_time"] / settings["timestep"]))
    settings["timestamps"] = np.linspace(0, settings["final_time"], settings["nsteps"])
    return settings


high_settings = create_settings(
    cells=(128,),
    dl=(0.05,),
    m=1,
    nbr_periods=10,
    timestep=0.0006,
    diag_dir="phare_outputs/whistler/high",
)

low_settings = create_settings(
    cells=(128,),
    dl=(0.8,),
    m=1,
    nbr_periods=10,
    timestep=0.077,
    diag_dir="phare_outputs/whistler/low",
)


def config(settings):
    sim = ph.Simulation(
        time_step=settings["timestep"],
        final_time=settings["final_time"],
        cells=settings["cells"],
        dl=settings["dl"],
        refinement="tagging",
        max_mhd_level=1,
        max_nbr_levels=1,
        hyper_resistivity=0.0,
        resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": settings["diag_dir"], "mode": "overwrite"},
        },
        strict=True,
        nesting_buffer=1,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="constant",
        limiter="",
        riemann="rusanov",
        mhd_timestepper="euler",
        hall=True,
        model_options=["MHDModel"],
    )

    k = settings["k"]
    modes = [1, 2, 4, 8]

    np.random.seed(0)
    phases = np.random.rand(len(modes))

    def density(x):
        return 1.0

    def vx(x):
        return 0.0

    def vy(x):
        return sum(-np.cos(k * x * m + phi) * 1e-2 * k for m, phi in zip(modes, phases))

    def vz(x):
        return sum(np.sin(k * x * m + phi) * 1e-2 * k for m, phi in zip(modes, phases))

    def bx(x):
        return 1.0

    def by(x):
        return sum(np.cos(k * x * m + phi) * 1e-2 for m, phi in zip(modes, phases))

    def bz(x):
        return sum(-np.sin(k * x * m + phi) * 1e-2 for m, phi in zip(modes, phases))

    def p(x):
        return 1.0

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=settings["timestamps"])

    return sim


class DispersionTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(DispersionTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(DispersionTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self, settings):
        # self.register_diag_dir_for_cleanup(settings["diag_dir"])
        Simulator(config(settings)).run().reset()
        return self


if __name__ == "__main__":
    startMPI()
    DispersionTest().test_run(high_settings).tearDown()
    DispersionTest().test_run(low_settings).tearDown()
