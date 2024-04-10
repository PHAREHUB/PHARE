#!/usr/bin/env python3

from pathlib import Path

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.run import Run

import numpy as np

import matplotlib as mpl

mpl.use("Agg")

from pyphare.cpp import cpp_lib
from tests.simulator import SimulatorTest


cpp = cpp_lib()
startMPI()

time_step = 0.005
final_time = .1
time_step_nbr = int(final_time / time_step)
timestamps = np.arange(0, final_time+.01, 0.05)
diag_dir = "phare_outputs/test_run"
plot_dir = Path(f"{diag_dir}_plots")
plot_dir.mkdir(parents=True, exist_ok=True)

def config():
    L = 0.5

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=(40, 40),
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=3,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="1",
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        }
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
        "nbr_part_per_cell": 44,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    pop = "protons"
    ph.ParticleDiagnostics(
        quantity="domain",
        write_timestamps=timestamps,
        population_name=pop,
    )
    ph.FluidDiagnostics(quantity="density", write_timestamps=timestamps, population_name=pop)

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
        for c in ["x","y","z"]:
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


def assert_file_exists_with_size_at_least(file, size=10000):
    path = Path(file)
    if not path.exists():
        raise FileNotFoundError("file not found: " + file)
    if path.stat().st_size < size:
        raise ValueError("file has unexpected size, possibly corrupt or not written properly: " + file)


class RunTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RunTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RunTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        sim = config()
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(sim).run()
        if cpp.mpi_rank() == 0:
            plot(diag_dir)

        for time in timestamps:
            for q in ["divb","Ranks","N","jz"]:
                assert_file_exists_with_size_at_least(plot_file_for_qty(q, time))

            for c in ["x","y","z"]:
                assert_file_exists_with_size_at_least(plot_file_for_qty(f"b{c}", time))

        cpp.mpi_barrier()


if __name__ == "__main__":
    import unittest
    unittest.main()
