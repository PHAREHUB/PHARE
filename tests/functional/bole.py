#
#
#


import os
import sys
import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI

np.set_printoptions(threshold=sys.maxsize)


ph.NO_GUI()


cells = (100, 100, 100)
# cells = (10, 10, 10)
dl = (0.1, 0.1, 0.1)

start_time = 0.2
name = "bowler"
diag_outputs = f"phare_outputs/test/{name}"
time_step_nbr = 1
time_step = 0.001
final_time = time_step * time_step_nbr + start_time

timestamps = []
# if time_step_nbr > 50:
#     dump_step = 1
#     nbr_dump_step = final_time / dump_step
#     timestamps = dump_step * np.arange(nbr_dump_step + 1)

print("timestamps=", timestamps)
plot_dir = Path(f"{diag_outputs}_plots")
plot_dir.mkdir(parents=True, exist_ok=True)


def config():
    sim = ph.Simulation(
        # largest_patch_size=50,
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        dl=dl,
        cells=cells,
        refinement="tagging",
        max_nbr_levels=1,
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_outputs, "mode": "overwrite"},
        },
        # restart_options={
        #     "dir": "checkpoints",
        #     "mode": "overwrite",
        #     "timestamps": [final_time],
        #     "restart_time": start_time,
        # },
    )

    def density(x, y, z):
        # L = sim.simulation_domain()[1]
        return 0.5

    def b(x, y, z, lx, ly, lz):
        L = sim.simulation_domain()[0]
        mid = L / 2

        X = x - mid  # .reshape(lx, ly, lz)
        Y = y - mid  # .reshape(lx, ly, lz)
        Z = z - mid  # .reshape(lx, ly, lz)

        U, V, W = -X, -Y, -Z

        # Normalize vectors to unit length
        magnitude = np.sqrt(U**2 + V**2 + W**2) + 1e-5  # Avoid division by zero

        # Normalize vectors (unit length)
        U /= magnitude
        V /= magnitude
        W /= magnitude

        # Define circular mask (radius = 25)
        radius = L * 0.4
        diff = 0.2

        # # outer mask
        # mask = X**2 + Y**2 + Z**2 <= (radius + diff) ** 2
        # U[~mask] = 0
        # V[~mask] = 0
        # W[~mask] = 0

        # # inner mask
        # mask = X**2 + Y**2 + Z**2 >= (radius - diff) ** 2
        # U[~mask] = 0
        # V[~mask] = 0
        # W[~mask] = 0

        # anything outside the horizontal middle plane = 0
        # mask = np.abs(Y) < 5  # & (X**2 + Z**2 >= (radius - diff) ** 2)
        mask = np.abs(Y) < 5 & (X**2 + Z**2 <= radius**2)
        U[~mask] = 0
        V[~mask] = 0
        W[~mask] = 0

        U *= 0.001
        V *= 0.001
        W *= 0.001

        return U, V, W

    def bx(x, y, z):
        return b(x, y, z, cells[0] + 5, cells[0] + 4, cells[0] + 4)[0]

    def by(x, y, z):
        return 0
        # return b(x, y, z, cells[0] + 4, cells[0] + 5, cells[0] + 4)[1]

    def bz(x, y, z):
        return b(x, y, z, cells[0] + 4, cells[0] + 4, cells[0] + 5)[2]

    def vxyz(x, y, z):
        return 0.0

    def vthxyz(x, y, z):
        return 0.00001

    C = "xyz"
    vvv = {
        **{f"vbulk{c}": vxyz for c in C},
        **{f"vth{c}": vthxyz for c in C},
        "nbr_part_per_cell": 100,
    }
    protons = {
        "charge": 1,
        "mass": 1,
        "density": density,
        **vvv,
        "init": {"seed": 12334},
    }
    ph.MaxwellianFluidModel(bx=bx, by=by, bz=bz, protons=protons)
    ph.ElectronModel(closure="isothermal", Te=0.0)

    # for quantity in ["density", "bulkVelocity"]:
    #     ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)
    # ph.FluidDiagnostics(
    #     quantity="density", write_timestamps=timestamps, population_name="protons"
    # )
    # for quantity in ["E", "B"]:
    #     ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    # ph.InfoDiagnostics(quantity="particle_count")  # defaults all coarse time steps

    return sim


def plot_file_for_qty(qty, time):
    return f"{plot_dir}/{name}_{qty}_t{time}.png"


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
            run.GetB(time, all_primal=False).plot(
                filename=plot_file_for_qty(f"b{c}", time),
                qty=f"B{c}",
                plot_patches=True,
            )
        for c in ["x", "y", "z"]:
            run.GetE(time, all_primal=False).plot(
                filename=plot_file_for_qty(f"e{c}", time),
                qty=f"E{c}",
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
    # Simulator(config()).run()
    Simulator(config()).setup(layout=3).initialize().run()

    if cpp.mpi_rank() == 0:
        plot(diag_outputs)
    cpp.mpi_barrier()


if __name__ == "__main__":
    startMPI()
    main()
