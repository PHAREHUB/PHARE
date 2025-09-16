#!/usr/bin/env python3
from dataclasses import field
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

final_time = 80
time_step = 0.2
timestamps = [final_time]
diag_dir = "phare_outputs/shock"


def config():
    cells = (800,)
    dl = (1.0,)

    sim = ph.Simulation(
        smallest_patch_size=15,
        # largest_patch_size=25,
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
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
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="constant",
        limiter="",
        riemann="rusanov",
        mhd_timestepper="euler",
        model_options=["MHDModel"],
    )

    def density(x):
        return np.where(x < (cells[0] * dl[0] / 2), 1, 0.125)

    def vx(x):
        return 0.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def bx(x):
        return 0.75

    def by(x):
        return np.where(x < (cells[0] * dl[0] / 2), 1, -1)

    def bz(x):
        return 0.0

    def p(x):
        return np.where(x < (cells[0] * dl[0] / 2), 1, 0.1)

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/shock_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetMHDrho(time).plot(
            filename=plot_file_for_qty(plot_dir, "rho", time), plot_patches=True
        )
        for c in ["x", "y"]:
            run.GetMHDV(time).plot(
                filename=plot_file_for_qty(plot_dir, f"v{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
        run.GetB(time).plot(
            filename=plot_file_for_qty(plot_dir, "by", time),
            plot_patches=True,
            qty="y",
        )
        run.GetMHDP(time).plot(
            filename=plot_file_for_qty(plot_dir, "p", time), plot_patches=True
        )


class ShockTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(ShockTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(ShockTest, self).tearDown()
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
    ShockTest().test_run().tearDown()
