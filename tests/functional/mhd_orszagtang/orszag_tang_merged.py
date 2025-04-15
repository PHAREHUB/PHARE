#!/usr/bin/env python3
import os

import numpy as np
import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.simulator.simulator import Simulator, startMPI

os.environ["PHARE_SCOPE_TIMING"] = "1"  # turn on scope timing

ph.NO_GUI()
cpp = cpp_lib()
startMPI()

diag_outputs = "phare_outputs"
time_step_nbr = 357
time_step = 0.0014
dt = 10 * time_step
timestamps = dt * np.arange(int(time_step_nbr / 10) + 1)


def config():
    cells = (128, 128)
    dl = (1.0 / cells[0], 1.0 / cells[1])

    sim = ph.Simulation(
        smallest_patch_size=15,
        largest_patch_size=25,
        time_step_nbr=time_step_nbr,
        time_step=time_step,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_mhd_level=1,
        max_nbr_levels=1,
        hyper_resistivity=0.0,
        resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        strict=True,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="weno3",
        limiter="",
        riemann="rusanov",
        mhd_timestepper="tvdrk3",
        model_options=["MHDModel"],
    )

    B0 = 1.0 / (np.sqrt(4.0 * np.pi))

    def density(x, y):
        return 25.0 / (36.0 * np.pi)

    def vx(x, y):
        return -np.sin(2.0 * np.pi * y)

    def vy(x, y):
        return np.sin(2.0 * np.pi * x)

    def vz(x, y):
        return 0.0

    def bx(x, y):
        return -B0 * np.sin(2.0 * np.pi * y)

    def by(x, y):
        return B0 * np.sin(4.0 * np.pi * x)

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 5.0 / (12.0 * np.pi)

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    # for quantity in ["rho", "V", "B", "P"]:
    #     ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def main():
    Simulator(config()).run()


if __name__ == "__main__":
    main()
