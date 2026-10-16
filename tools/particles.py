import sys
import numpy as np

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare
from matplotlib import pyplot as plt


def main():
    run = Run(sys.argv[1])
    # dim = run.simulation().ndim

    def stat(ndim, arrs, key):
        vavg = np.average(arrs)
        vmin = np.min(arrs)
        vmax = np.max(arrs)

        print(key)
        print("avg", vavg)
        print("min", vmin)
        print("max", vmax)
        print("")

    def lvl0(hier, time):
        lvl = hier.level(0)
        # print(lvl, list(lvl.values()))
        for patch in lvl:
            for pd in patch:
                particles = pd.dataset
                stat(hier.sim.ndim, particles.v, "v")
                stat(hier.sim.ndim, particles.weights, "weights")
                stat(hier.sim.ndim, particles.charges, "charges")
                # stat(hier.sim.ndim, particles.deltas, "delta")
                print(particles.weights)

    for time in run.all_times()["B"]:
        pops = run.all_pops()
        for pop in pops:
            lvl0(run.GetParticles(time, pop), time)


if __name__ == "__main__":
    main()
