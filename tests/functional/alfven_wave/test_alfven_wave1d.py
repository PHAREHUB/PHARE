#!/usr/bin/env python3

from pybindlibs import cpp
import alfven_wave1d as aw1d
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharein import global_vars as gv
import numpy as np
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
import matplotlib.pyplot as plt


import matplotlib as mpl
import shutil
import os
import time
mpl.use('Agg')

def plot(bhier):
    times = np.sort(np.asarray(list(bhier.time_hier.keys())))

    components  =("B_x", "B_y", "B_z")
    ylims = ((0,1.5),(-0.25,0.25),(-0.25,0.25))

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

    simulator = Simulator(gv.sim)
    tic = time.time()
    simulator.initialize()
    toc = time.time()
    print("initialization took {}".format(float(toc-tic)))


    times =  np.arange(0, gv.sim.final_time, simulator.timeStep())
    perf = np.zeros_like(times)
    for it,t in enumerate(times):
        tic = time.time()
        simulator.advance()
        toc = time.time()
        perf[it] = toc-tic
        print("t = {}  -  {:6.5f}sec  - total {:7.4}sec".format(t, perf[it], np.sum(perf)))


    print("mean advance time = {}".format(np.mean(perf)))
    print("total advance time = {}".format(np.sum(perf)))

    print(perf)

    if cpp.mpi_rank() == 0:
        b = hierarchy_from(h5_filename="phare_outputs/EM_B.h5")
        plot(b)






if __name__ == "__main__":
    main()
