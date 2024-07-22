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
from tools.python3 import plotting as m_plotting

mpl.use("Agg")

LOAD_BALANCE = os.getenv("LOAD_BALANCE", "True").lower() in ("true", "1", "t")

cpp = cpp_lib()
startMPI()

cells = (800, 800)
time_step = 0.005
final_time = 50
timestamps = np.arange(0, final_time + time_step, final_time / 5)

if cpp.mpi_rank() == 0:
    print(LOAD_BALANCE, "diag timestamps:", timestamps)

diag_dir = "phare_outputs/harris_lb"
if not LOAD_BALANCE:
    diag_dir = "phare_outputs/harris"

plot_dir = Path(f"{diag_dir}_plots")
plot_dir.mkdir(parents=True, exist_ok=True)


def config():
    L = 0.5

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=2,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="1",
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        restart_options={
            "dir": "checkpoints",
            "mode": "overwrite",
            "timestamps": timestamps,
            # "restart_time": 0.0,
        },
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

    def vxyz(x, y):
        return 0.0

    def vthxyz(x, y):
        return np.sqrt(T(x, y))

    vvv = {**{f"vbulk{c}": vxyz for c in "xyz"}, **{f"vth{c}": vthxyz for c in "xyz"}}

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": 1, "density": density, **vvv}
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    ph.FluidDiagnostics(
        quantity="density", write_timestamps=timestamps, population_name="protons"
    )
    ph.InfoDiagnostics(quantity="particle_count")

    if LOAD_BALANCE:
        ph.LoadBalancer(active=True, auto=True, mode="nppc", tol=0.05)

    return sim


def plot_file_for_qty(qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty("divb", time),
            plot_patches=True,
            vmin=1e-11,
            vmax=2e-10,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty("Ranks", time),
            plot_patches=True,
        )
        run.GetN(time, pop_name="protons").plot(
            filename=plot_file_for_qty("N", time),
            plot_patches=True,
        )
        for c in ["x", "y", "z"]:
            run.GetB(time).plot(
                filename=plot_file_for_qty(f"b{c}", time),
                qty=f"B{c}",
                plot_patches=True,
            )
        run.GetJ(time).plot(
            filename=plot_file_for_qty("jz", time),
            qty="Jz",
            plot_patches=True,
            vmin=-2,
            vmax=2,
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
            plot(diag_dir)
        m_plotting.plot_run_timer_data(diag_dir, cpp.mpi_rank())
        cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    HarrisTest().test_run().tearDown()
