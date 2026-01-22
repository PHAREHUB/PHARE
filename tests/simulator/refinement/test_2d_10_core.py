#!/usr/bin/env python3

"""
  Complex test to force Yee Lattice centering discrepancies due to SAMRAI Schedule/Algorithm
    grouping of quantity components with distinct geometries.

  Should be run with 10 cores
    mpirun -n 10 tests/simulator/refinement/test_2d_10_core.py
"""


import numpy as np
import pyphare.pharein as ph
import pyphare.core.box as boxm
from pyphare.simulator.simulator import Simulator, startMPI

from pyphare.cpp import cpp_lib
from tests.simulator.test_advance import AdvanceTestBase

cpp = cpp_lib()
test = AdvanceTestBase(rethrow=True)  # change to False for debugging images
L0_diags = "phare_outputs/test_x_homo_0"
L0L1_diags = "phare_outputs/test_x_homo_1"

ph.NO_GUI()


def config(diag_outputs, model_init={}, refinement_boxes=None):
    ph.global_vars.sim = None

    sim = ph.Simulation(
        smallest_patch_size=10,
        largest_patch_size=20,
        time_step_nbr=1,
        final_time=0.001,
        # boundary_types="periodic",
        cells=(50, 100),
        dl=(0.40, 0.40),
        # refinement="tagging",
        # max_nbr_levels = 3,
        refinement_boxes=refinement_boxes,
        hyper_resistivity=0.0050,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {
                "dir": diag_outputs,
                "mode": "overwrite",
                "fine_dump_lvl_max": 10,
            },
        },
    )

    def density(x, y):
        return 1.0

    def bx(x, y):
        return 0.1

    def by(x, y):
        return 0.2

    def bz(x, y):
        return 1.0

    def T(x, y):
        return 1.0

    def vx(x, y):
        return 1.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": 100,
        "init": model_init,
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": 1, "density": density, **vvv}
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    dt = 1 * sim.time_step
    nt = sim.final_time / dt + 1
    timestamps = dt * np.arange(nt)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["charge_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def get_time(path, time=None, datahier=None):
    if time is not None:
        time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    datahier = hierarchy_from(h5_filename=path + "/EM_E.h5", times=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path + "/EM_B.h5", times=time, hier=datahier)
    return datahier


def get_hier(path):
    return get_time(path)


def make_fig(hier, fig_name, ilvl, extra_collections=[]):
    if cpp.mpi_rank() == 0:
        l0_in_l1 = [boxm.refine(p.box, 2) for p in hier.level(0).patches]
        collections = [
            {
                "boxes": l0_in_l1,
                "facecolor": "grey",
            }
        ]
        if 1 in hier.levels():
            l1_over_l0 = [p.box for p in hier.level(1).patches]
            collections += [
                {
                    "boxes": l1_over_l0,
                    "facecolor": "yellow",
                }
            ]
        hier.plot_2d_patches(
            1,
            collections=collections + extra_collections,
        ).savefig(fig_name + ".png")


def check_hier(hier):
    for ilvl, lvl in hier.levels().items():
        for patch in lvl.patches:
            assert all(patch.box.shape > 5)
    return hier


def post_advance(new_time):
    if cpp.mpi_rank() == 0:
        L0_datahier = check_hier(get_hier(L0_diags))
        L0L1_datahier = check_hier(get_hier(L0L1_diags))
        extra_collections = []
        errors = test.base_test_overlaped_fields_are_equal(L0L1_datahier, new_time)

        # The test commented out, see warning in test_advance.py
        # errors = (
        #     test.base_test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
        #         L0_datahier, L0L1_datahier
        #     )
        # )
        print(f"errors {type(errors)}")
        if isinstance(errors, list):
            extra_collections += [
                {
                    "boxes": errors,
                    "facecolor": "black",
                }
            ]
        make_fig(L0L1_datahier, L0L1_diags.split("/")[-1], 1, extra_collections)


def main():
    import random

    startMPI()
    rando = random.randint(0, int(1e10))
    Simulator(config(L0_diags, {"seed": rando})).run().reset()
    refinement_boxes = {"L0": {"B0": [(5, 40), (20, 59)]}}
    sim = config(L0L1_diags, {"seed": rando}, refinement_boxes)
    Simulator(sim, post_advance=post_advance).run()


if __name__ == "__main__":
    main()
