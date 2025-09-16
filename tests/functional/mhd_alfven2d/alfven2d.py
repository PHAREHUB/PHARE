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

time_step = 0.002
final_time = 1.0  # time for one period
timestamps = [0.0, final_time]
diag_dir = "phare_outputs/alfven2d"


def config():
    alpha = 30.0 * np.pi / 180.0
    cosalpha = np.cos(alpha)
    sinalpha = np.sin(alpha)

    cells = (100, 100)
    dl = ((1.0 / cells[0]) * 1 / cosalpha, (1.0 / cells[1]) * 1 / sinalpha)

    sim = ph.Simulation(
        smallest_patch_size=15,
        # largest_patch_size=25,
        time_step=time_step,
        final_time=final_time,
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
        reconstruction="constant",
        limiter="",
        riemann="rusanov",
        mhd_timestepper="euler",
        model_options=["MHDModel"],
    )

    def density(x, y):
        return 1.0

    def vx(x, y):
        return -0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * sinalpha

    def vy(x, y):
        return 0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * cosalpha

    def vz(x, y):
        return 0.1 * np.cos(2 * np.pi * (x * cosalpha + y * sinalpha))

    def bx(x, y):
        return (
            cosalpha
            - 0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * sinalpha
        )

    def by(x, y):
        return (
            sinalpha
            + 0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * cosalpha
        )

    def bz(x, y):
        return 0.1 * np.cos(2 * np.pi * (x * cosalpha + y * sinalpha))

    def p(x, y):
        return 0.1

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/alfven2d_{qty}_t{time}.png"


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
        for c in ["x", "y"]:
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


class AlfvenTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(AlfvenTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(AlfvenTest, self).tearDown()
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


def main():
    Simulator(config()).run()


if __name__ == "__main__":
    startMPI()
    AlfvenTest().test_run().tearDown()
