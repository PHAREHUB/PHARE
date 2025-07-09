#!/usr/bin/env python3

import os
import numpy as np
import matplotlib as mpl
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest

mpl.use("Agg")

cpp = cpp_lib()


cells = (200, 100)
time_step = 0.005
final_time = 50
timestamps = np.arange(0, final_time + time_step, final_time / 5)
diag_dir = "phare_outputs/harris"


def config():
    L = 0.5

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=2,
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
    )

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

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

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 0.7
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": 100,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 12334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["mass_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["density", "pressure_tensor"]:
        ph.FluidDiagnostics(
            quantity=quantity, write_timestamps=timestamps, population_name="protons"
        )

    ph.InfoDiagnostics(quantity="particle_count")

    ph.LoadBalancer(active=True, auto=True, mode="nppc", tol=0.05)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty(plot_dir, "divb", time),
            plot_patches=True,
            vmin=1e-11,
            vmax=2e-10,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty(plot_dir, "Ranks", time), plot_patches=True
        )
        run.GetN(time, pop_name="protons").plot(
            filename=plot_file_for_qty(plot_dir, "N", time), plot_patches=True
        )
        for c in ["x", "y", "z"]:
            run.GetB(time).plot(
                filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
        run.GetJ(time).plot(
            filename=plot_file_for_qty(plot_dir, "jz", time),
            qty="z",
            plot_patches=True,
            vmin=-2,
            vmax=2,
        )


class HarrisTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(HarrisTest, self).__init__(*args, **kwargs)
        self.simulator = None
        self.plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
        self.plot_dir.mkdir(parents=True, exist_ok=True)

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
            plot(diag_dir, self.plot_dir)
        cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    startMPI()
    HarrisTest().test_run().tearDown()
