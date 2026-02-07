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

final_time = 0.15
time_step = 0.0003
timestamps = np.arange(0, final_time + time_step, final_time / 5)
diag_dir = "phare_outputs/rotor"


def config():
    cells = (100, 100)
    dl = (1.0 / cells[0], 1.0 / cells[1])

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
            "fine_dump_lvl_max": 10,
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


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/rotor_{qty}_t{time}.png"


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
        run.GetMHDP(time).plot(
            filename=plot_file_for_qty(plot_dir, "p", time), plot_patches=True
        )


class RotorTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RotorTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RotorTest, self).tearDown()
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
    RotorTest().test_run().tearDown()
