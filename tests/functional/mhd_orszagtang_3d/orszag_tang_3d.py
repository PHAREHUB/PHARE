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

final_time = 0.5
time_step = 0.00025
diag_dir = "phare_outputs/orszag_tang"

time_step_nbr = int(final_time / time_step)
start_dump_time = 0.0
dumpfrequency = 200
dt = dumpfrequency * time_step
timestamps = (
    dt * np.arange(int((final_time - start_dump_time) / dt) + 1) + start_dump_time
)


def config():
    cells = (64, 64, 64)
    dl = (1.0 / cells[0], 1.0 / cells[1], 1.0 / cells[2])

    sim = ph.Simulation(
        smallest_patch_size=15,
        # largest_patch_size=25,
        time_step_nbr=time_step_nbr,
        time_step=time_step,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_mhd_level=2,
        max_nbr_levels=2,
        hyper_resistivity=0.0,
        resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
        nesting_buffer=1,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="linear",
        limiter="vanleer",
        riemann="rusanov",
        mhd_timestepper="tvdrk2",
        model_options=["MHDModel"],
    )

    B0 = 1.0 / (np.sqrt(4.0 * np.pi))
    epsp = 0.2

    def density(x, y, z):
        return 25.0 / (36.0 * np.pi)

    def vx(x, y, z):
        return -(1 + epsp * np.sin(2 * np.pi * z)) * np.sin(2.0 * np.pi * y)

    def vy(x, y, z):
        return (1 + epsp * np.sin(2 * np.pi * z)) * np.sin(2.0 * np.pi * x)

    def vz(x, y, z):
        return epsp * np.sin(2 * np.pi * z)

    def bx(x, y, z):
        return -B0 * np.sin(2.0 * np.pi * y)

    def by(x, y, z):
        return B0 * np.sin(4.0 * np.pi * x)

    def bz(x, y, z):
        return 0.0

    def p(x, y, z):
        return 5.0 / (12.0 * np.pi)

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


class OrszagTangTest3d(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(OrszagTangTest3d, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(OrszagTangTest3d, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        Simulator(config()).run().reset()
        return self


if __name__ == "__main__":
    startMPI()
    OrszagTangTest3d().test_run().tearDown()
