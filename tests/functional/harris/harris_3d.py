#!/usr/bin/env python3
import os
import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI

os.environ["PHARE_SCOPE_TIMING"] = "0"  # turn on scope timing


ph.NO_GUI()
cpp = cpp_lib()
startMPI()

cells = (50, 50, 50)
dl = (0.2, 0.2, 0.2)

diag_outputs = "phare_outputs/test/harris/3d"
time_step_nbr = 1000
time_step = 0.001
final_time = time_step * time_step_nbr

timestamps = [0, final_time]
if time_step_nbr > 100:
    dt = 100 * time_step
    nt = final_time / dt + 1
    timestamps = dt * np.arange(nt)

plot_dir = Path(f"{diag_outputs}_plots")
plot_dir.mkdir(parents=True, exist_ok=True)


def config():
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        dl=dl,
        cells=cells,
        refinement="tagging",
        max_nbr_levels=2,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        # strict=True,
    )

    def density(x, y, z):
        L = sim.simulation_domain()[1]
        return (
            0.2
            + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y, z):
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

    def bx(x, y, z):
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

    def bz(x, y, z):
        return 0.1

    def b2(x, y, z):
        return bx(x, y, z) ** 2 + by(x, y, z) ** 2 + bz(x, y, z) ** 2

    def T(x, y, z):
        K = 1
        temp = 1.0 / density(x, y, z) * (K - b2(x, y, z) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vxyz(x, y, z):
        return 0.0

    def vthxyz(x, y, z):
        return np.sqrt(T(x, y, z))

    C = "xyz"
    vvv = {
        **{f"vbulk{c}": vxyz for c in C},
        **{f"vth{c}": vthxyz for c in C},
        "nbr_part_per_cell": 50,
    }
    protons = {"charge": 1, "density": density, **vvv, "init": {"seed": 12334}}
    ph.MaxwellianFluidModel(bx=bx, by=by, bz=bz, protons=protons)
    ph.ElectronModel(closure="isothermal", Te=0.0)

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)
    ph.FluidDiagnostics(
        quantity="density", write_timestamps=timestamps, population_name="protons"
    )
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    ph.InfoDiagnostics(quantity="particle_count")  # defaults all coarse time steps

    return sim


def plot_file_for_qty(qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir):
    run = Run(diag_dir)
    for time in timestamps:
        # run.GetDivB(time).plot(
        #     filename=plot_file_for_qty("divb", time),
        #     plot_patches=True,
        #     vmin=1e-11,
        #     vmax=2e-10,
        # )
        # run.GetRanks(time).plot(
        #     filename=plot_file_for_qty("Ranks", time),
        #     plot_patches=True,
        # )
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
        # run.GetJ(time).plot(
        #     filename=plot_file_for_qty("jz", time),
        #     qty="Jz",
        #     plot_patches=True,
        #     vmin=-2,
        #     vmax=2,
        # )


def main():
    Simulator(config()).run()

    # TODO3D nico : I removed the plot because the hierarchy plot3d has been removed
    # TODO3D we need to find a way to plot the data

    # if cpp.mpi_rank() == 0:
    #     plot(diag_outputs)

    cpp.mpi_barrier()


if __name__ == "__main__":
    main()
