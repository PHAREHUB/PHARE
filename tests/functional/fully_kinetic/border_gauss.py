#!/usr/bin/env python3
"""
BORDER-GAUSS (docs/superpowers/specs/p4/border_gauss.md): on the edge nodes of a fine level, the
Gauss residual F = div E - c^2 rho must not follow the coarse E the ghost fills put on the faces
outside the level, nor accumulate the band refreshes.

Run: fine_gauss_2d.py's 2D two-level T2 plasma (box over root cells 8..23, 2 x 2 fine patches),
no cleaning, so that F changes only through the E update and the particles. Its electron velocity
wave peaks on the box edges (x = 8 and 24 root cells, y alike): a mean current crosses every edge.

The solver sets the faces outside the level after every E fill so that F = 0 on the edge nodes,
with the charge at E's time (fk_electric.hpp, gaussOnLevelBorder): refilled from the coarse E
instead, F moves by O(J dt / dx) per substep and drifts. The residual of a node reads its owner
patch's E and rho (fine_gauss_2d.py, fine_residual).

Gate 1 (--check): max |F| / scale over the edge nodes (corners and patch junctions included) at
the end of every substep that does not end a coarse step (here the mid step dumps), below
GATE = 1e-6: round-off and the single precision of the dumped E (4e-8 to 6e-8 of the scale).
Interior: the change of F over a substep, round-off too, printed. With N > 1 filter passes the
band nodes are not filtered again (fk_electric.hpp, depositFKCurrent): not exact, not gated.

Report, at the coarse times: the messenger replaces the band by the coarser particles refined at
the end of each coarse step (lastStep), after the last E fill of the step, which moves the charge
of the edge nodes until the next fill sets the faces again: the edge F of a coarse time dump is
that charge jump. Printed: the rms over the edge nodes of F_k - F_1 at coarse step k and its growth
exponent p (rms ~ (k - 1)^p, least squares): ~0 for jumps that do not add up; 0.54 when the faces
kept their values (ad46f273b), 0.88 when refilled from the coarse E.
"""

import argparse
import shutil
import sys
from argparse import Namespace

import numpy as np

import fine_gauss_2d as fg2
from fine_gauss import DT, RATIO

GATE = 1e-6


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--filter-passes", type=int, default=1)
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--ppc", type=int, default=20)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    fg2.diag_root = "phare_outputs/fully_kinetic_border_gauss"
    case = Namespace(
        order=1,
        filter_passes=args.filter_passes,
        cleaning=0.0,
        ppc=args.ppc,
        steps=args.steps,
        largest_patch_size=16,
    )
    diag_dir, times, seconds = fg2.simulate(case, 0)
    s0 = fg2.scale(diag_dir, case.order)
    cuts, _ = fg2.junctions(diag_dir, times[-1])
    edge = fg2.classify(cuts) >= 0

    def residual(t):
        return fg2.fine_residual(diag_dir, t, case.order, case.filter_passes)[0] / s0

    substep = [residual(t + DT / RATIO) for t in times[:-1]]
    edge_value = max(np.abs(f[edge]).max() for f in substep)
    inner_change = max(
        np.abs(residual(t + DT / RATIO) - residual(t))[~edge].max() for t in times[1:-1]
    )
    coarse = [residual(t)[edge] for t in times[1:]]
    walk = np.array([np.sqrt(np.mean((f - coarse[0]) ** 2)) for f in coarse[1:]])
    k = np.arange(2, len(coarse) + 1)
    exponent = np.polyfit(np.log(k - 1), np.log(np.maximum(walk, 1e-30)), 1)[0]
    print(
        f"border gauss N={case.filter_passes}: max |F| / scale on the edge at the mid steps"
        f" {edge_value:.2e} (gate {GATE:.0e}), max |dF| / scale over a substep inside"
        f" {inner_change:.2e}; at coarse steps 1..{args.steps}: max |F| / scale on the edge"
        f" {max(np.abs(f).max() for f in coarse):.2e}, rms edge (F_k - F_1) / scale, k = 2.."
        f"{args.steps}: {' '.join(f'{w:.2e}' for w in walk)}, ~ (k - 1)^{exponent:.2f};"
        f" run {seconds:.1f} s",
        flush=True,
    )
    shutil.rmtree(diag_dir, ignore_errors=True)
    if args.check and edge_value > GATE:
        sys.exit(f"border gauss: edge |F| / scale {edge_value:.2e} > {GATE:.0e} at a mid step")


if __name__ == "__main__":
    main()
