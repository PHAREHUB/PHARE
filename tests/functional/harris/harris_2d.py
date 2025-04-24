#!/usr/bin/env python3

import numpy as np

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI

from tests.simulator import SimulatorTest


ph.NO_GUI()
cpp = cpp_lib()


diag_outputs = "phare_outputs/test/harris/2d"
time_step_nbr = 1000
time_step = 0.001
final_time = time_step * time_step_nbr


def default_timestamps():
    dt = 10 * time_step
    nt = final_time / dt + 1
    timestamps = dt * np.arange(nt)
    return timestamps


def default_setup():
    startMPI()

    return ph.Simulation(
        smallest_patch_size=15,
        largest_patch_size=25,
        time_step_nbr=time_step_nbr,
        time_step=time_step,
        # boundary_types="periodic",
        cells=(200, 400),
        dl=(0.2, 0.2),
        refinement="tagging",
        max_nbr_levels=1,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        strict=True,
    )


def config(sim=None, timestamps=None, seed=12334):
    if sim is None:
        sim = default_setup()
    if timestamps is None:
        timestamps = default_timestamps()

    def density(x, y):
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        return (w5 * x0 * w3) + (-w5 * x0 * w4)

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        v1 = -1
        v2 = 1.0
        return (
            v1
            + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
            + (-w5 * y1 * w3)
            + (+w5 * y2 * w4)
        )

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 1
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
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": seed}},
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    ph.InfoDiagnostics(quantity="particle_count")  # defaults all coarse time steps

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


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

    def test_run(self, diag_dir=None, sim=None):
        diag_dir = diag_dir if diag_dir else diag_outputs
        sim = sim if sim else config()
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(sim).run().reset()
        return self

    def plot(self, timestamps, diag_dir, plot_dir):
        run = self.getRun(diag_dir)
        for time in timestamps:
            run.GetDivB(time).plot(
                filename=plot_file_for_qty(plot_dir, "divb", time),
                plot_patches=True,
                vmin=1e-11,
                vmax=2e-10,
            )
            run.GetRanks(time).plot(
                filename=plot_file_for_qty(plot_dir, "Ranks", time),
                plot_patches=True,
            )
            run.GetN(time, pop_name="protons").plot(
                filename=plot_file_for_qty(plot_dir, "N", time),
                plot_patches=True,
            )
            for c in ["x", "y", "z"]:
                run.GetB(time).plot(
                    filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                    qty=f"{c}",
                    plot_patches=True,
                )
            run.GetJ(time).plot(
                filename=plot_file_for_qty(plot_dir, "jz", time),
                qty="z",
                plot_patches=True,
                vmin=-2,
                vmax=2,
            )

    def scope_timing(self, diag_dir):
        try:
            from tools.python3 import plotting as m_plotting

            m_plotting.plot_run_timer_data(diag_dir, cpp.mpi_rank())
        except ImportError:
            print("Phlop not found - install with: `pip install phlop`")
        except FileNotFoundError:
            print("Phlop installed but not active")


if __name__ == "__main__":
    HarrisTest().test_run().tearDown()
