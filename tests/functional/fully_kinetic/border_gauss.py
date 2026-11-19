#!/usr/bin/env python3
"""
BORDER-GAUSS (docs/superpowers/specs/p4/border_gauss.md): on the edge nodes of a fine level, the
Gauss residual F = div E - c^2 rho must change only as the deposit does, not as the coarse E does.

Run: fine_gauss_2d.py's 2D two-level T2 plasma (box over root cells 8..23, 2 x 2 fine patches),
no cleaning, so that F changes only through the E update and the particles: with cleaning 0 the
Esirkepov deposit conserves charge on every node the fine E is advanced around. Its electron
velocity wave peaks on the box edges (x = 8 and 24 root cells, y alike): a mean current crosses
every edge.

Measured: per fine node of the box, the change of F over the first fine substep of each coarse
step from step 1 on (coarse time -> mid step, fine dumps on every substep): the band is renewed at
the coarse time, before it, so in that substep only the fine particles (domain, band, the leavers
handed to the bands) change rho. Inside the box E is the fine Ampere's, the change is round-off.
On an edge node, E reads its first outside normal face: refilled from the coarse E (before the
fix) it changes F by O(J dt / dx); advanced by the fine Ampere (the fix) the change is round-off.

Gate (--check): max |dF| / scale over the edge nodes (corners and patch junctions included)
below GATE = 1e-6. E is dumped in single precision: the interior, and the edge with the fix, sit
at 4e-8 to 6e-8 of the scale (order 1, N = 0 and 1); before the fix the edge moves by 1e-2 to
5e-2 of the scale per substep. With N > 1 filter passes the band nodes are not filtered again
(fk_electric.hpp, depositFKCurrent): not exact, not gated. A level border face is a ghost of
every patch along it, each patch advancing its own from its own band: the residual of a node
reads its owner patch's faces (fine_gauss_2d.py, fine_e).
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
    parser.add_argument("--steps", type=int, default=4)
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

    changes = np.array([residual(t + DT / RATIO) - residual(t) for t in times[1:-1]])
    edge_change = np.abs(changes[:, edge]).max()
    inner_change = np.abs(changes[:, ~edge]).max()
    drift = [residual(t)[edge].mean() for t in times[1:]]
    print(
        f"border gauss N={case.filter_passes}: max |dF| / scale over a fine substep, edge "
        f"{edge_change:.2e}, interior {inner_change:.2e} (gate {GATE:.0e} on the edge); mean edge"
        f" F / scale at coarse steps 1..{args.steps}: {' '.join(f'{d:+.3f}' for d in drift)};"
        f" run {seconds:.1f} s",
        flush=True,
    )
    shutil.rmtree(diag_dir, ignore_errors=True)
    if args.check and edge_change > GATE:
        sys.exit(f"border gauss: edge |dF| / scale {edge_change:.2e} > {GATE:.0e}")


if __name__ == "__main__":
    main()
