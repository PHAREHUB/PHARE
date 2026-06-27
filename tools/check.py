#!/usr/bin/env python3

import subprocess
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.core import gridlayout
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest

ph.NO_GUI()

cpp = cpp_lib()


cells = (200, 100)
time_step = 0.005
final_time = .01
timestamps = [0, final_time/2,final_time]

def _run(s):
    args = s.split(" ")
    result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode > 0:
        raise RuntimeError("subprocess error:", result.stderr)
    return result.stdout.strip()

commit_hash = _run("git rev-parse --short HEAD")
commit_date = _run("git log --no-show-signature -1 --format=%cI")

diag_dir = f"phare_outputs/harris/4/{commit_date}_{commit_hash}"
# diag_dir = "phare_outputs/harris"

def config():
    L = 0.5

    sim = ph.Simulation(
        interp_order=1,
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
            "options": {"dir": diag_dir, "mode": "overwrite", "fine_dump_lvl_max": 0},
        },
        strict=True,
        nesting_buffer=1,
        tag_buffer=3,
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
    # for quantity in ["mass_density", "bulkVelocity"]:
    #     ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)
    # for quantity in ["density", "pressure_tensor"]:
    #     ph.FluidDiagnostics(
    #         quantity=quantity, write_timestamps=timestamps, population_name="protons"
    #     )
    ph.InfoDiagnostics(quantity="particle_count")
    ph.LoadBalancer(active=True, auto=True, mode="nppc", tol=0.05)

    return sim


def box_as_filename_string(box):
    L = [str(v) for v in box.lower]
    U = [str(v) for v in box.upper]
    return f"L{"_".join(L)}_U{"_".join(U)}"

def plot_file_for_qty(plot_dir, qty, time, extra=""):
    return f"{plot_dir}/harris_t{"{:.10f}".format(time)}_{qty}_{extra}.png"

def plot_time(run, time, diag_dir, plot_dir, **kwargs):
    qty = "Ey"
    hier = run.GetE(time, all_primal=False)

    for ilvl, lvl in hier.levels(time).items():
        hier.plot(filename=plot_file_for_qty(
                                plot_dir, qty, time, f"L{ilvl}"
                            ),
                            plot_patches=True,
                            # vmin=0,
                            # vmax=+1e-16,
                            # vmin=-2,
                            # vmax=2,
                            qty=qty,
                            levels=(ilvl,),
                            dpi=1000,
                        )

    hier = (list(hier.time_hier.values()))[0]

    for patch in hier[1]:
        if patch.box.lower[0] != 0:
            continue
        print("patch.box", patch.box)
        fig, ax = plt.subplots()
        pdat = patch.patch_datas[qty]
        layout = pdat.layout
        im = ax.pcolormesh(
            layout.yeeCoordsFor(qty, "x", withGhosts=True),
            layout.yeeCoordsFor(qty, "y", withGhosts=True),
            pdat.dataset[:].T,
            cmap=kwargs.get("cmap", "Spectral_r"),
            vmin=kwargs.get("vmin", np.min(pdat.dataset) - 1e-6),
            vmax=kwargs.get("vmax", np.max(pdat.dataset) + 1e-6),
        )
        plt.colorbar(im, ax=ax)
        fig.savefig(plot_file_for_qty(plot_dir, qty, time, f"L1_{box_as_filename_string(patch.box)}"), dpi=kwargs.get("dpi", 200))
        # return


def plot(diag_dir, plot_dir, **kwargs):
    run = Run(diag_dir)

    for time in run.all_times()["B"]:
        plot_time(run, time, diag_dir, plot_dir)




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
        # self.register_diag_dir_for_cleanup(diag_dir)
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
