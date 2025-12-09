#!/usr/bin/env python3

import os
import sys
import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy.hierarchy_utils import diff_hierarchy


from tests.simulator.test_advance import AdvanceTestBase

ph.NO_GUI()

cpp = cpp_lib()

diag_dir = sys.argv[1] if len(sys.argv) > 1 else "phare_outputs/advance/test_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts_6/2/2/3"
print("diag_dir", os.getcwd(), diag_dir)


def plot_file_for_qty(plot_dir, qty, time, extra=""):
    return f"{plot_dir}/harris_t{"{:.10f}".format(time)}_{qty}_{extra}.png"





def diff(new_time):
    print("diff", new_time)

    plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
    plot_dir.mkdir(parents=True, exist_ok=True)
    run = Run(diag_dir)
    ranks = run.GetRanks(new_time)

    for ilvl in range(ranks.levelNbr()):
        ranks.plot(
            filename=plot_file_for_qty(plot_dir, f"ranks", new_time, f"L{ilvl}"),
            plot_patches=True,
            levels=(ilvl,),
            dpi=2000,
        )

    differ = diff_hierarchy(run.GetE(new_time, all_primal=False))
    print("E max: ", differ.max())
    print("Ex max: ", differ.max("Ex"))
    print("Ey max: ", differ.max("Ey"))
    print("Ez max: ", differ.max("Ez"))

    if differ.has_non_zero():
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

    differ = diff_hierarchy(run.GetB(new_time, all_primal=False))
    print("B max: ", differ.max())
    print("Bx max: ", differ.max("Bx"))
    print("By max: ", differ.max("By"))
    print("Bz max: ", differ.max("Bz"))
    if differ.has_non_zero():
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

    differ = diff_hierarchy(run.GetNi(new_time))
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

    differ = diff_hierarchy(run.GetMassDensity(new_time))
    print("ion mass rho max: ", differ.max())

    if differ.has_non_zero():
        for ilvl in range(differ.levelNbr()):
            differ.plot(
                filename=plot_file_for_qty(plot_dir, f"ionMass", new_time, f"L{ilvl}"),
                plot_patches=True,
                vmin=0,
                vmax=+1e-16,
                levels=(ilvl,),
                dpi=1000,
            )

def get_time(path, time=None, datahier=None):
    if time is not None:
        time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    datahier = hierarchy_from(h5_filename=path + "/EM_E.h5", times=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path + "/EM_B.h5", times=time, hier=datahier)
    return datahier


test = AdvanceTestBase(rethrow=True)  # change to False for debugging images

def check_time(new_time):
    if cpp.mpi_rank() == 0:
        try:
            diff(new_time)
            hier = get_time(diag_dir, new_time)
            errors = test.base_test_overlaped_fields_are_equal(hier, new_time)
            if isinstance(errors, list):
                print("\n\n!!ERROR AT TIME: ", new_time)
        except FileNotFoundError as e:
            print("File not found for diag!\n", e)
        except KeyError as e:
            err = str(e)

            if not "Unable to synchronously open object" in err:  # no diag for time
                import traceback

                print(f"Exception caught: \n{e}")
                print(traceback.format_exc())

    cpp.mpi_barrier()


def check_diags():
    run = Run(diag_dir)
    for time in run.all_times()["B"]:
        check_time(time)



if __name__ == "__main__":
    startMPI()
    check_diags()

