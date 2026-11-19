#!/usr/bin/env python3
"""
BORDER-RANKS: a 2D fine level over several MPI ranks with a rank that holds no level border node
(levelBorderNodesOf empty on that rank): either no fine patch, or fine patches all away from the
level border. The step of a refined level must make the same collective calls on every rank (the
charge deposit's border sums and ghost fills before the second Gauss face set, solver_fk.hpp): a
rank skipping them because its own border is empty hangs the others. Under ctest the TIMEOUT
stops a hang; the run must end on all ranks.

SAMRAI gives every rank a piece of the fine level: here 9 ranks over a box of 56 x 56 fine cells
put rank 0 at its centre, fine cells 22..41 (3, 4, 6 and 8 ranks do not: every rank touches the
border; 5 ranks crash in SAMRAI's root partition). Before the per-level guard it hung in its
first step. The precondition is checked after the run on rank 0, from the patch boxes and ranks
of the E dump at t = 0: without such a rank the test proves nothing and fails.

Run: mpirun -n 9 python3 border_ranks.py (~5 s)
"""

import shutil
import sys

import h5py
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

CELLS, DL, DT, STEPS = 32, 0.04, 0.002, 4
BOX = (2, 29)  # root cells, both included, both directions
REACH = 3  # levelBorderNodesOf: border nodes within field ghost width - 1 (4 - 1, FK order 1)
DIAG = "phare_outputs/fully_kinetic_border_ranks"


def uniform(value):
    return lambda x, y: value + 0.0 * x


def wave(amplitude, direction):
    return lambda x, y: amplitude * np.sin(2 * np.pi * (x, y)[direction] / (CELLS * DL))


def simulate():
    lo, hi = BOX
    sim = ph.Simulation(
        time_step=DT,
        time_step_nbr=STEPS,
        cells=[CELLS, CELLS],
        dl=[DL, DL],
        interp_order=1,
        model_options=["FullyKineticModel"],
        c_over_va=10.0,
        current_filter_passes=1,
        divergence_cleaning=0.5,
        refinement_boxes={"L0": {"B0": [(lo, lo), (hi, hi)]}},
        smallest_patch_size=8,
        largest_patch_size=8,
        diag_options={"format": "phareh5", "options": {"dir": DIAG, "mode": "overwrite"}},
        strict=True,
    )
    common = {"density": uniform(1.0), "nbr_part_per_cell": 10}
    ph.MaxwellianFluidModel(
        bx=uniform(1.0),
        by=wave(0.1, 0),
        bz=wave(0.1, 1),
        protons={
            "charge": 1,
            "mass": 1.0,
            **{f"vth{k}": uniform(np.sqrt(0.03)) for k in "xyz"},
            "init": {"seed": 2026},
            **common,
        },
        electrons={
            "charge": -1,
            "mass": 0.04,
            "vbulkx": wave(0.05, 0),
            "vbulky": wave(0.05, 1),
            "vbulkz": uniform(0.0),
            **{f"vth{k}": uniform(np.sqrt(0.03 / 0.04)) for k in "xyz"},
            "init": {"seed": 2027},
            **common,
        },
    )
    ph.ElectromagDiagnostics(quantity="E", write_timestamps=[0.0])
    Simulator(sim).run().reset()


def borderless_ranks(size):
    """the ranks with no fine patch within REACH cells of the level border, and the fine patch
    count per rank"""
    lo, hi = 2 * BOX[0], 2 * BOX[1] + 1  # first and last fine cell
    on_border = dict.fromkeys(range(size), False)
    counts = dict.fromkeys(range(size), 0)
    with h5py.File(f"{DIAG}/EM_E.h5", "r") as f:
        for patch in f["t"][next(iter(f["t"]))]["pl1"].values():
            rank = int(patch.attrs["mpi_rank"])
            lower, upper = (np.array(patch.attrs[k]) for k in ("lower", "upper"))
            counts[rank] += 1
            on_border[rank] |= bool(any(lower - REACH <= lo) or any(upper + REACH >= hi))
    return [r for r, border in on_border.items() if not border], counts


def main():
    simulate()
    cpp.mpi_barrier()
    if cpp.mpi_rank() != 0:
        return
    borderless, counts = borderless_ranks(cpp.mpi_size())
    shutil.rmtree(DIAG, ignore_errors=True)
    print(
        f"border ranks: fine patches per rank {counts}, ranks without level border nodes"
        f" {borderless}; the run ended on every rank",
        flush=True,
    )
    if not borderless:
        sys.exit("border ranks: every rank holds level border nodes, nothing tested")


if __name__ == "__main__":
    main()
