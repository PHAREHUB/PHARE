#!/usr/bin/env python3
"""
BORDER-NOISE: on the edge of a fine level, the tangential E must not hold a grid-scale mode the
inside does not have. fine_gauss_2d.py's 2-level T2 plasma (box over root cells 8..23, 2 x 2 fine
patches), filter 1 pass, cleaning 0.5 with the boundary mask (run146's options).

The filter removes the grid-scale part of J, so the particles do not hold back a grid-scale
longitudinal E. Inside, Gauss's law (charge conservation, the cleaning) pins it. On the edge
nodes nothing did: the faces outside take any residual there (gaussOnLevelBorder), and the
cleaning skips them. The tangential E on the edge then random walked (t^0.4 over 3000 steps,
x10 the inside), and through Gauss the faces outside with it (run146, t 15: 1.7 rms). The fix is
one Marder step along the edge before the outside faces are set (marderAlongLevelBorder).

Measure, per E dump at coarse steps 100..STEPS every 25: on the four sides (3 fine nodes off the
corners) the cell-to-cell part of the tangential E, rms of (a_{j-1} - 2 a_j + a_{j+1}) / 4 along
the side, on the edge nodes and 3 nodes inside. Gate: mean edge / inside <= RATIO_MAX.
Measured (Release, ppc 20): 2.93 on 412e334b4, 1.13 with marderAlongLevelBorder.
"""

import argparse
import shutil
import sys

import numpy as np

import fine_gauss_2d as fg2
from fine_gauss import C_OVER_VA, DL, DT, MASS_E, RATIO, T

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

STEPS, EVERY, FROM = 200, 25, 100
RATIO_MAX = 1.5
diag_dir = "phare_outputs/fully_kinetic_border_noise"


def simulate(ppc):
    lo, hi = fg2.BOX
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=DT,
        time_step_nbr=STEPS,
        cells=[fg2.CELLS] * 2,
        dl=[DL] * 2,
        interp_order=1,
        model_options=["FullyKineticModel"],
        c_over_va=C_OVER_VA,
        current_filter_passes=1,
        divergence_cleaning=0.5,
        boundary_cleaning_mask=True,
        refinement_boxes={"L0": {"B0": [(lo, lo), (hi, hi)]}},
        smallest_patch_size=8,
        largest_patch_size=16,
        diag_options={"format": "phareh5", "options": {"dir": diag_dir, "mode": "overwrite"}},
        strict=True,
    )
    u = fg2.uniform
    common = {"density": u(1.0), "nbr_part_per_cell": ppc}
    ph.MaxwellianFluidModel(
        bx=u(1.0),
        by=fg2.wave(0.1, 0),
        bz=fg2.wave(0.1, 1),
        protons={
            "charge": 1,
            "mass": 1.0,
            **{f"vbulk{k}": u(0.0) for k in "xyz"},
            **{f"vth{k}": u(np.sqrt(T)) for k in "xyz"},
            "init": {"seed": 2026},
            **common,
        },
        electrons={
            "charge": -1,
            "mass": MASS_E,
            "vbulkx": fg2.wave(0.05, 0),
            "vbulky": fg2.wave(0.05, 1),
            "vbulkz": u(0.0),
            **{f"vth{k}": u(np.sqrt(T / MASS_E)) for k in "xyz"},
            "init": {"seed": 2027},
            **common,
        },
    )
    times = DT * np.arange(FROM, STEPS + 1, EVERY)
    ph.ElectromagDiagnostics(quantity="E", write_timestamps=times)
    Simulator(sim).run().reset()
    return times


def fine_e(t):
    """{(2 x / dx_f, 2 y / dx_f): value} of the fine Ex, Ey, ghosts included"""
    dl = DL / RATIO
    hier = hierarchy_from(h5_filename=f"{diag_dir}/EM_E.h5", times=f"{t:.10f}")
    E = {"x": {}, "y": {}}
    for patch in hier.level(1, t).patches:
        for c in "xy":
            pdata = patch.patch_datas[f"E{c}"]
            x, y = pdata.meshgrid()
            for px, py, v in zip(x.ravel(), y.ravel(), pdata.dataset[:].ravel()):
                E[c].setdefault((round(2 * px / dl), round(2 * py / dl)), float(v))
    return E


def cell_to_cell(a):
    a = np.asarray(a)[3:-3]
    return np.sqrt(np.mean(((a[:-2] - 2 * a[1:-1] + a[2:]) / 4) ** 2))


def edge_and_inside(E, depth=3):
    """mean over the four sides of the cell-to-cell part of the tangential E on the edge nodes
    and `depth` nodes inside"""
    lo, n = fg2.fine_nodes()
    along = np.arange(lo, lo + n - 1)
    out = {0: [], depth: []}
    for edge, inward in ((lo, 1), (lo + n - 1, -1)):
        for normal in (0, 1):
            tangential = "y" if normal == 0 else "x"
            for k in out:
                i = 2 * (edge + inward * k)
                keys = [(i, 2 * j + 1) if normal == 0 else (2 * j + 1, i) for j in along]
                out[k].append(cell_to_cell([E[tangential][key] for key in keys]))
    return np.mean(out[0]), np.mean(out[depth])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--ppc", type=int, default=20)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    times = simulate(args.ppc)
    rows = [edge_and_inside(fine_e(t)) for t in times]
    ratios = [e / i for e, i in rows]
    print(
        "border noise: cell-to-cell tangential E, edge / 3 nodes inside, steps "
        + " ".join(
            f"{round(t / DT)}: {e:.2e}/{i:.2e}" for t, (e, i) in zip(times, rows)
        )
        + f"; mean ratio {np.mean(ratios):.2f} (gate {RATIO_MAX})",
        flush=True,
    )
    shutil.rmtree(diag_dir, ignore_errors=True)
    if args.check and np.mean(ratios) > RATIO_MAX:
        sys.exit(f"border noise: edge / inside {np.mean(ratios):.2f} > {RATIO_MAX}")


if __name__ == "__main__":
    main()
