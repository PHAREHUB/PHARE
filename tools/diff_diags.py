#!/usr/bin/env python3

import os
import sys
import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy.hierarchy_utils import overlap_diff_hierarchy


from tests.simulator.test_advance import AdvanceTestBase

ph.NO_GUI()

DO_PLOTS = False  # global plot skip
DOMAIN_ONLY = True
cpp = cpp_lib()

diag_dir = sys.argv[1] if len(sys.argv) > 1 else "phare_outputs/harris_3d"
print("diag_dir", os.getcwd(), diag_dir)


def plot_file_for_qty(plot_dir, qty, time, extra=""):
    return f"{plot_dir}/harris_t{"{:.10f}".format(time)}_{qty}_{extra}.png"


def diff_ranks(run, plot_dir, new_time):

    if DO_PLOTS:
        ranks = run.GetRanks(new_time)
        for ilvl in range(ranks.levelNbr()):
            ranks.plot(
                filename=plot_file_for_qty(plot_dir, f"ranks", new_time, f"L{ilvl}"),
                plot_patches=True,
                levels=(ilvl,),
                dpi=2000,
            )


def diff_E(run, plot_dir, new_time):
    differ = overlap_diff_hierarchy(
        run.GetE(new_time, all_primal=False), domain_only=DOMAIN_ONLY
    )
    print("E max: ", differ.max())
    print("E patch shapes: ", differ.min_max_patch_shape())
    print("Ex max: ", differ.max("Ex"))
    print("Ey max: ", differ.max("Ey"))
    print("Ez max: ", differ.max("Ez"))

    if DO_PLOTS and differ.has_non_zero():
        for c in ["x", "y", "z"]:
            for ilvl in range(differ.levelNbr()):
                differ.plot(
                    filename=plot_file_for_qty(
                        plot_dir, f"diffE{c}", new_time, f"L{ilvl}"
                    ),
                    plot_patches=True,
                    vmin=0,
                    vmax=+1e-16,
                    qty=f"E{c}",
                    levels=(ilvl,),
                    dpi=2000,
                )


def diff_B(run, plot_dir, new_time):
    differ = overlap_diff_hierarchy(
        run.GetB(new_time, all_primal=False), domain_only=DOMAIN_ONLY
    )
    print("B max: ", differ.max())
    print("Bx max: ", differ.max("Bx"))
    print("By max: ", differ.max("By"))
    print("Bz max: ", differ.max("Bz"))
    if DO_PLOTS and differ.has_non_zero():
        for c in ["x", "y", "z"]:
            for ilvl in range(differ.levelNbr()):
                differ.plot(
                    filename=plot_file_for_qty(
                        plot_dir, f"diffB{c}", new_time, f"L{ilvl}"
                    ),
                    plot_patches=True,
                    vmin=0,
                    vmax=+1e-16,
                    qty=f"B{c}",
                    levels=(ilvl,),
                    dpi=2000,
                )


def diff_current_density(run, plot_dir, new_time):
    differ = overlap_diff_hierarchy(run.GetNi(new_time), domain_only=DOMAIN_ONLY)
    print("ion charge rho max: ", differ.max())

    if differ.has_non_zero():
        for ilvl in range(differ.levelNbr()):
            differ.plot(
                filename=plot_file_for_qty(
                    plot_dir, f"ionCharge", new_time, f"L{ilvl}"
                ),
                plot_patches=True,
                vmin=0,
                vmax=+1e-16,
                levels=(ilvl,),
                dpi=2000,
            )


def diff_mass_density(run, plot_dir, new_time):
    differ = overlap_diff_hierarchy(run.GetMassDensity(new_time), domain_only=DOMAIN_ONLY)
    print("ion mass rho max: ", differ.max())

    if DO_PLOTS and differ.has_non_zero():
        for ilvl in range(differ.levelNbr()):
            differ.plot(
                filename=plot_file_for_qty(plot_dir, f"ionMass", new_time, f"L{ilvl}"),
                plot_patches=True,
                vmin=0,
                vmax=+1e-16,
                levels=(ilvl,),
                dpi=1000,
            )


def check_diag(run, plot_dir, new_time, fn):
    try:
        fn(run, plot_dir, new_time)

    except FileNotFoundError as e:
        print("File not found for", fn.__name__)
    except Exception as e:
        import traceback

        print(f"Unknown Exception caught: \n{e}")
        print(traceback.format_exc())


def check_time(run, plot_dir, new_time):
    if cpp.mpi_rank() == 0:
        check_diag(run, plot_dir, new_time, diff_E)
        check_diag(run, plot_dir, new_time, diff_B)
        check_diag(run, plot_dir, new_time, diff_current_density)
        check_diag(run, plot_dir, new_time, diff_mass_density)
    cpp.mpi_barrier()


def check_diags():

    plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
    plot_dir.mkdir(parents=True, exist_ok=True)
    run = Run(diag_dir)
    for time in run.all_times()["B"]:
        check_time(run, plot_dir, time)


if __name__ == "__main__":
    startMPI()
    check_diags()
