#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.run import Run

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
from pyphare.cpp import cpp_lib

cpp = cpp_lib()


def config(interp_order):
    """Configure the simulation

    This function defines the Simulation object,
    user initialization model and diagnostics.
    """
    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=20,
        time_step=0.005,  # number of time steps (not specified if time_step and final_time provided)
        final_time=30,  # simulation final time (not specified if time_step and time_step_nbr provided)
        boundary_types="periodic",  # boundary condition, string or tuple, length == len(cell) == len(dl)
        cells=2500,  # integer or tuple length == dimension
        dl=0.2,  # mesh size of the root level, float or tuple
        # max_nbr_levels=1,          # (default=1) max nbr of levels in the AMR hierarchy
        hyper_resistivity=0.01,
        nesting_buffer=0,
        interp_order=interp_order,
        # refinement_boxes = {"L0":{"B0":[(125,), (750,)]}},
        diag_options={
            "format": "phareh5",
            "options": {"dir": "shock_{}".format(interp_order), "mode": "overwrite"},
        },
    )

    def density(x):
        L = sim.simulation_domain()[0]
        v1 = 1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(x, L * 0.2, 1) - S(x, L * 0.8, 1))

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def bx(x):
        return 0.0

    def by(x):
        L = sim.simulation_domain()[0]
        v1 = 0.125
        v2 = 4.0
        return v1 + (v2 - v1) * (S(x, L * 0.2, 1) - S(x, L * 0.8, 1))

    def bz(x):
        return 0.0

    def T(x):
        return 0.1

    def vx(x):
        L = sim.simulation_domain()[0]
        v1 = 0.0
        v2 = 0.0
        return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vthx(x):
        return T(x)

    def vthy(x):
        return T(x)

    def vthz(x):
        return T(x)

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": 1, "density": density, **vvv}
    )

    ph.ElectronModel(closure="isothermal", Te=0.12)

    dt = 10 * sim.time_step
    nt = sim.final_time / dt + 1
    timestamps = dt * np.arange(nt)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["charge_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def main():
    import os
    import subprocess
    import glob
    import shlex

    for interp_order in (1, 2, 3):
        sim = config(interp_order)
        Simulator(sim).run()

        if cpp.mpi_rank() == 0:
            dt = 10 * sim.time_step
            nt = sim.final_time / dt + 1
            times = dt * np.arange(nt)
            r = Run("shock_{}".format(interp_order))
            for it, t in enumerate(times):
                fig, ax = plt.subplots()
                B = r.GetB(t, merged=True)
                title = "interp order {} - t = {:06.3f}".format(interp_order, t)
                x = B["By"][1][0]
                By = B["By"][0]
                ax.plot(x, By(x), color="k")
                ax.set_title(title)
                ax.set_ylim((-0.2, 5))
                ax.set_xlim((0, 250))
                fig.savefig(
                    "shock_{}/shock_By_{}_{:04d}.png".format(
                        interp_order, interp_order, it
                    )
                )
                plt.close(fig)
            cmd = shlex.split(
                "ffmpeg -r 10 -y -pattern_type glob -i 'shock_{}/shock_By_{}_*.png' -c:v libx264 -crf 0 shock_interp{}.mp4".format(
                    interp_order, interp_order, interp_order
                )
            )
            subprocess.call(cmd)

        ph.global_vars.sim = None
    pngs = glob.glob("shock*/*.png")
    for png in pngs:
        os.remove(png)

    if cpp.mpi_rank() == 0:
        t = 30
        runs = [Run(f"shock_{i+1}") for i in range(3)]
        fig, ax = plt.subplots()
        colors = ["k", "r", "b"]
        for r, color, interp_order in zip(runs, colors, (1, 2, 3)):
            print(r.path)
            B = r.GetB(t, merged=True)
            x = B["By"][1][0]
            By = B["By"][0]
            ax.plot(x, By(x), color=color, label=f"interp order {interp_order}")
        title = "interp order {} - t = {:06.3f}".format(interp_order, t)
        ax.set_title(title)
        ax.set_ylim((-0.2, 5))
        ax.set_xlim((0, 250))
        ax.legend()
        fig.savefig("shock_By.png")


if __name__ == "__main__":
    startMPI()
    main()
