#!/usr/bin/env python3

import os
import sys
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy import hierarchy_utils as hootils


from tests.simulator.test_advance import AdvanceTestBase


if len(sys.argv) == 1:
    print("diag dir expected as first argument")
    sys.exit(1)

diag_dir = sys.argv[1]

ph.NO_GUI()
DO_PLOTS = True  # global plot skip
cpp = cpp_lib()
dpi=300


def plot_file_for_qty(plot_dir, qty, time, extra=""):
    return f"{plot_dir}/harris_t{"{:.10f}".format(time)}_{qty}_{extra}.png"


def diff_ranks(run, plot_dir, time):
    if DO_PLOTS:
        ranks = run.GetRanks(time)
        for ilvl in ranks.levels(time).keys():
            ranks.plot(
                filename=plot_file_for_qty(plot_dir, f"ranks", time, f"L{ilvl}"),
                plot_patches=True,
                levels=(ilvl,),
                dpi=dpi,
            )


def diff_B(run, plot_dir, time):
    differ = hootils.overlap_diff_hierarchy(run.GetB(time, all_primal=False), time)
    print("B max: ", hootils.max_from(differ, time))
    print("Bx max: ", hootils.max_from(differ, time, qty="Bx"))
    print("By max: ", hootils.max_from(differ, time, qty="By"))
    print("Bz max: ", hootils.max_from(differ, time, qty="Bz"))
    print("B patch shapes: ", hootils.min_max_patch_shape(differ, time))
    if DO_PLOTS and hootils.has_non_zero(differ, time):
        for c in ["x", "y", "z"]:
            for ilvl in differ.levels(time).keys():
                differ.plot(
                    filename=plot_file_for_qty(plot_dir, f"diffB{c}", time, f"L{ilvl}"),
                    plot_patches=True,
                    vmin=0,
                    vmax=+1e-16,
                    qty=f"B{c}",
                    levels=(ilvl,),
                    dpi=dpi,
                )


def diff_E(run, plot_dir, time):
    differ = hootils.overlap_diff_hierarchy(run.GetE(time, all_primal=False), time)

    print("E max: ", hootils.max_from(differ, time))
    print("Ex max: ", hootils.max_from(differ, time, qty="Ex"))
    print("Ey max: ", hootils.max_from(differ, time, qty="Ey"))
    print("Ez max: ", hootils.max_from(differ, time, qty="Ez"))
    print("E patch shapes: ", hootils.min_max_patch_shape(differ, time))

    if DO_PLOTS and hootils.has_non_zero(differ, time):
        for c in ["x", "y", "z"]:
            for ilvl in differ.levels(time).keys():
                differ.plot(
                    filename=plot_file_for_qty(plot_dir, f"diffE{c}", time, f"L{ilvl}"),
                    plot_patches=True,
                    vmin=0,
                    vmax=+1e-16,
                    qty=f"E{c}",
                    levels=(ilvl,),
                    dpi=dpi,
                )


def diff_current_density(run, plot_dir, time):
    differ = hootils.overlap_diff_hierarchy(run.GetNi(time), time)
    print("ion charge rho max: ", hootils.max_from(differ, time))

    if hootils.has_non_zero(differ, time):
        for ilvl in differ.levels(time).keys():
            differ.plot(
                filename=plot_file_for_qty(plot_dir, f"ionCharge", time, f"L{ilvl}"),
                plot_patches=True,
                vmin=0,
                vmax=+1e-16,
                levels=(ilvl,),
                dpi=dpi,
            )


def diff_mass_density(run, plot_dir, time):
    differ = hootils.overlap_diff_hierarchy(run.GetMassDensity(time), time)
    print("ion mass rho max: ", hootils.max_from(differ, time))

    if DO_PLOTS and hootils.has_non_zero(differ, time):
        for ilvl in differ.levels(time).keys():
            differ.plot(
                filename=plot_file_for_qty(plot_dir, f"ionMass", time, f"L{ilvl}"),
                plot_patches=True,
                vmin=0,
                vmax=+1e-16,
                levels=(ilvl,),
                dpi=dpi
            )


def check_diag(run, plot_dir, time, fn):
    try:
        fn(run, plot_dir, time)

    except FileNotFoundError as e:
        print("File not found for", fn.__name__)
    except Exception as e:
        import traceback

        print(f"Unknown Exception caught: \n{e}")
        print(traceback.format_exc())


def check_time(run, plot_dir, time):
    if cpp.mpi_rank() == 0:
        print("checking time", time)
        check_diag(run, plot_dir, time, diff_ranks)
        check_diag(run, plot_dir, time, diff_E)
        check_diag(run, plot_dir, time, diff_B)
        check_diag(run, plot_dir, time, diff_current_density)
        check_diag(run, plot_dir, time, diff_mass_density)
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
