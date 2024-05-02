# various plotting tools for PHARE development, NOT physics!
#

import sys
import logging
import argparse
import matplotlib.pyplot as plt

from pathlib import Path

from pyphare.pharesee.run import Run
import tools.python3.phloping as phloping

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

logging.getLogger("PIL").setLevel(logging.INFO)  # noise
logging.getLogger("h5py").setLevel(logging.INFO)  # noise
plt.set_loglevel(level="warning")  # noise


def plot_run_timer_data(diag_dir=None, rank=0):
    if diag_dir is None:  # assume cli
        parser = argparse.ArgumentParser()
        parser.add_argument("-d", "--dir", default=".", help="Diagnostics directory")
        diag_dir = parser.parse_args().dir
    run = Run(diag_dir, single_hier_for_all_quantities=True)
    res = phloping.file_parser(run, rank, Path(f".phare_times.{rank}.txt"))
    fig, ax = plt.subplots()
    L0X = res.time_steps_for_L(0)
    ax.plot(L0X, res.normalised_times_for_L(0), ":", label="L0 times", color="black")

    if res.n_levels() > 1:
        ax.plot(
            res.time_steps_for_L(1),
            res.normalised_times_for_L(1),
            "g:",
            label="L1 times",
        )
    if res.n_levels() > 2:
        ax.plot(
            res.time_steps_for_L(2),
            res.normalised_times_for_L(2),
            "b:",
            label="L2 times",
        )
    ax.legend()
    plt.ylabel("time in ns")
    plt.xlabel(f"timestep {res.sim.time_step}")
    ax.set_xticks([L0X[0], L0X[-1]])
    fig.savefig(f"run_timer.{rank}.png")


if __name__ == "__main__":
    # use first cli argument to pick plotting function
    if len(sys.argv) > 1:
        fn = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        globals()[fn]()
    else:
        print("available functions:")
        fns = [k for k in list(globals().keys()) if k.startswith("plot_")]
        for k in fns:
            print(" ", k)
