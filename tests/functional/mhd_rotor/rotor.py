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

diag_outputs = "phare_outputs/low"
final_time = 0.03
time_step = 0.0003
time_step_nbr = int(final_time / time_step)

dumpfrequency = 1
dt = dumpfrequency * time_step
timestamps = dt * np.arange(int(time_step_nbr / dumpfrequency) + 1)


def config():
    cells = (100, 100)
    dl = (1.0 / cells[0], 1.0 / cells[1])

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
            "options": {"dir": diag_outputs, "mode": "overwrite"},
            "fine_dump_lvl_max": 10,
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

    B0 = 5 / (np.sqrt(4 * np.pi))
    v0 = 2

    r0 = 0.1
    r1 = 0.115

    def r(x, y):
        return np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2)

    def f(r):
        return (r1 - r) / (r1 - r0)

    def density(x, y):
        r_ = r(x, y)
        f_ = f(r_)

        rho_values = np.where(r_ <= r0, 10.0, np.where(r_ < r1, 1.0 + 9.0 * f_, 1.0))
        return rho_values

    def vx(x, y):
        r_ = r(x, y)
        f_ = f(r_)

        vx_values = np.where(
            r_ <= r0,
            -v0 * (y - 0.5) / r0,
            np.where(r_ < r1, -f_ * v0 * (y - 0.5) / r_, 0.0),
        )
        return vx_values

    def vy(x, y):
        r_ = r(x, y)
        f_ = f(r_)

        vy_values = np.where(
            r_ <= r0,
            v0 * (x - 0.5) / r0,
            np.where(r_ < r1, f_ * v0 * (x - 0.5) / r_, 0.0),
        )
        return vy_values

    def vz(x, y):
        return 0.0

    def bx(x, y):
        return B0

    def by(x, y):
        return 0.0

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 1.0

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def main():
    Simulator(config()).run()


if __name__ == "__main__":
    main()
