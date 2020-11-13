#!/usr/bin/env python3

from pybindlibs import cpp

import pyphare.pharein as ph
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics,FluidDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
import numpy as np

from pyphare.pharesee.hierarchy import hierarchy_from


import matplotlib.pyplot as plt


import matplotlib as mpl
import shutil
import os
import time
mpl.use('Agg')




def config():
    """ Configure the simulation

    This function defines the Simulation object,
    user initialization model and diagnostics.
    """
    Simulation(
        smallest_patch_size=10,
        largest_patch_size=20,
        time_step_nbr=2000,        # number of time steps (not specified if time_step and final_time provided)
        final_time=25.,             # simulation final time (not specified if time_step and time_step_nbr provided)
        boundary_types="periodic", # boundary condition, string or tuple, length == len(cell) == len(dl)
        cells=80,                # integer or tuple length == dimension
        dl=0.3,                  # mesh size of the root level, float or tuple
        max_nbr_levels=3,          # (default=1) max nbr of levels in the AMR hierarchy
        refinement_boxes={"L0": {"B0": [(10, ), (40, )]},
                          "L1":{"B0":[(30,),(60,)]}},
        diag_options={"format": "phareh5", "options": {"dir": "phare_outputs","mode":"overwrite"}}
    )


    def density(x):
        return 1.


    def S(x,x0,l):
        return 0.5*(1+np.tanh((x-x0)/l))


    def bx(x):
        return 0.


    def by(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()[0]
        v1=-1
        v2=1.
        return v1 + (v2-v1)*(S(x,L*0.25,1) -S(x, L*0.75, 1))


    def bz(x):
        return 0.5


    def b2(x):
        return bx(x)**2 + by(x)**2 + bz(x)**2


    def T(x):
        K = 1
        return 1/density(x)*(K - b2(x)*0.5)


    def vx(x):
        return 0.


    def vy(x):
        return 0.


    def vz(x):
        return 0.


    def vthx(x):
        return T(x)


    def vthy(x):
        return T(x)


    def vthz(x):
        return T(x)


    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz
    }

    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 1337}}
    )

    ElectronModel(closure="isothermal", Te=0.12)



    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time +sim.time_step, 5*sim.time_step)



    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )


    for quantity in ["density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            )




def plot(bhier):
    times = np.sort(np.asarray(list(bhier.time_hier.keys())))

    components  =("B_y", "B_z")
    ylims = ((-1.25, 1.25),(0.,1.0))

    for component,ylim in zip(components,ylims):
        for it,t in enumerate(times):
            fig,ax = plt.subplots(figsize=(10,6))
            for il,level in bhier.levels(t).items():
                patches = level.patches
                if il == 0:
                    marker="+"
                    alpha=1
                    ls='-'
                else:
                    marker='o'
                    alpha=0.4
                    ls='none'

                for ip, patch in enumerate(patches):
                    val   = patch.patch_datas["EM_"+component].dataset[:]
                    x_val = patch.patch_datas["EM_"+component].x
                    label="${}$ level {} patch {}".format(component,il,ip)
                    ax.plot(x_val, val, label=label,
                            marker=marker, alpha=alpha, ls=ls)
                    ax.set_ylim(ylim)

            ax.legend(ncol=4)
            ax.set_title("t = {:05.2f}".format(t))
            fig.savefig("{}_{:04d}.png".format(component,it))
            plt.close(fig)



def main():
    config()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()


    if cpp.mpi_rank() == 0:
        b = hierarchy_from(h5_filename="phare_outputs/EM_B.h5")
        plot(b)

if __name__=="__main__":
    main()
