# various plotting tools for PHARE development, NOT physics!
#

import sys
import logging
import argparse
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

from pyphare.pharesee.run import Run
import tools.python3.run_timer as rt

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
    rtf = rt.file_parser(run, rank, Path(f".phare_times.{rank}.bin"))
    fig, ax = plt.subplots()
    L0X = rtf.time_steps_for_L(0)
    ax.plot(L0X, rtf.normalised_times_for_L(0), "b--", label="L0 times")
    ax.plot(
        rtf.time_steps_for_L(1), rtf.normalised_times_for_L(1), "g:", label="L1 times"
    )
    ax.legend()
    plt.ylabel("time in ns")
    plt.xlabel(f"timestep {rtf.sim.time_step}")
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
