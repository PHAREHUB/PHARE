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

cells = (600, 300)
time_step = 0.005 / 16
final_time = 50
timestamps = np.arange(0, final_time + time_step, final_time / 5)
diag_dir = "phare_outputs/mhd_harris"

hall = True


def config():
    L = 0.5

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=(0.10, 0.10),
        refinement="tagging",
        max_mhd_level=1,
        max_nbr_levels=1,
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
        hall=hall,
        model_options=["MHDModel"],
    )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

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

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 1.0 - (bx(x, y) ** 2 + by(x, y) ** 2) / 2.0

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty(plot_dir, "divb", time),
            plot_patches=True,
            vmin=-1e-11,
            vmax=1e-11,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty(plot_dir, "Ranks", time), plot_patches=True
        )
        run.GetMHDrho(time).plot(
            filename=plot_file_for_qty(plot_dir, "rho", time), plot_patches=True
        )
        for c in ["x", "y", "z"]:
            run.GetMHDV(time).plot(
                filename=plot_file_for_qty(plot_dir, f"v{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
            run.GetB(time).plot(
                filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
        run.GetMHDP(time).plot(
            filename=plot_file_for_qty(plot_dir, "p", time), plot_patches=True
        )
        if hall:
            run.GetJ(time).plot(
                filename=plot_file_for_qty(plot_dir, "jz", time),
                qty="z",
                plot_patches=True,
            )


class HarrisTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(HarrisTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(HarrisTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()
        if cpp.mpi_rank() == 0:
            plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
            plot_dir.mkdir(parents=True, exist_ok=True)
            plot(diag_dir, plot_dir)
        cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    startMPI()
    HarrisTest().test_run().tearDown()
