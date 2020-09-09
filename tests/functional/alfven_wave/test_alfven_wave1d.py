#!/usr/bin/env python3


import alfven_wave1d as aw1d
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
import numpy as np
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
import matplotlib.pyplot as plt


def main():
    simulator = Simulator(gv.sim)
    simulator.initialize()

    times =  np.arange(0,gv.sim.final_time, simulator.timeStep())
    for it,t in enumerate(times):
        simulator.diagnostics().dump(timestamp=simulator.currentTime(),
                                     timestep=simulator.timeStep())
        simulator.advance()


    b = hierarchy_from(h5_filename="phare_outputs/EM_B.h5")

    times = np.sort(np.asarray(list(b.time_hier.keys())))
    fig,ax = plt.subplots(figsize=(10,6))
    t=0
    for il,level in b.levels(t).items():
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
            val   = patch.patch_datas["EM_B_x"].dataset[:]
            x_val = patch.patch_datas["EM_B_x"].x
            label="$B_x$ level {} patch {}".format(il,ip)
            ax.plot(x_val, val, label=label,
                    marker=marker, alpha=alpha, ls=ls)
            ax.set_ylim((0, 1.5))

    ax.legend(ncol=4)
    ax.set_title("t = {:05.2f}".format(t))





if __name__ == "__main__":
    main()
