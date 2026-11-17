#!/usr/bin/env python3
"""
The gate of the sawtooth fix (docs/superpowers/specs/p4/sawtooth_fix_spec.md): fine_gauss.py's
1D run and measurement at order 1, no filter, s = 0.5, 1600 particles per cell, 100 steps, seed 0.

Before the fix, the domain particles leaving into the level ghost band were erased each fine step
while the band was renewed from the coarser level each coarse step: the jump (coarse time mean -
mid step mean) of F / scale at the boundary nodes, left / right, was +0.518 / +0.511, and the
two-phase mean -0.078 / -0.180.
With the fix (keeping them in the band until its renewal): +0.008 / +0.008 and +0.017 / -0.075.

Gate:
- |jump| <= JUMP_GATE on both sides. The no-sawtooth floor has sd 0.02 over 8 seeds and reaches
  0.040 at seed 0: 0.1 is 5 sd;
- the two-phase mean does not grow: |after| <= |before| + TWO_PHASE_SD (that floor's sd). Removing
  the mid step deficit moves it, by up to ~0.15 (UNVERIFIED): its shift is printed.
"""

import shutil
import time
from argparse import Namespace

import numpy as np

import fine_gauss as fg

CASE = Namespace(
    order=1,
    filter_passes=0,
    cleaning=0.5,
    ppc=1600,
    steps=100,
    vth_factor_protons=1.0,
    vth_factor_electrons=1.0,
    no_strict=False,
)
SEED = 0
JUMP_GATE, TWO_PHASE_SD = 0.1, 0.035
TWO_PHASE_BEFORE = (-0.078, -0.180)  # left / right boundary node
DUMPS_GB = 1.9  # every particle at every step (fine_gauss.py's diagnostics)


def main():
    fg.diag_root = "phare_outputs/fully_kinetic_sawtooth_gate"  # not fine_gauss.py's own runs
    # a full disk (a RAM-backed /tmp) cuts the dumps short: refuse to start rather than read them
    free_gb = shutil.disk_usage(".").free / 1e9
    if free_gb < 2 * DUMPS_GB:
        raise RuntimeError(f"sawtooth gate: {free_gb:.1f} GB free here, {2 * DUMPS_GB} needed")
    start = time.perf_counter()
    diag_dir, times, seconds = fg.simulate(CASE, SEED)
    s0 = fg.scale(diag_dir, CASE.order)

    def mean_from_jump_step(at):
        rows = [fg.fine_residual(diag_dir, t, CASE.order, CASE.filter_passes) for t in at]
        return np.array(rows)[fg.JUMP_FROM :].mean(axis=0) / s0

    coarse_mean = mean_from_jump_step(times)
    mid_mean = mean_from_jump_step(times[:-1] + fg.DT / fg.RATIO)
    jump, two_phase = coarse_mean - mid_mean, 0.5 * (coarse_mean + mid_mean)

    ends = (0, -1)
    bound = [abs(before) + TWO_PHASE_SD for before in TWO_PHASE_BEFORE]
    shift = [two_phase[e] - before for e, before in zip(ends, TWO_PHASE_BEFORE)]
    print(
        f"sawtooth gate: jump at the boundary nodes {jump[0]:+.3f} / {jump[-1]:+.3f} (|jump| <="
        f" {JUMP_GATE}); two-phase mean {two_phase[0]:+.3f} / {two_phase[-1]:+.3f} (|mean| <="
        f" {bound[0]:.3f} / {bound[1]:.3f}), shift from before the fix {shift[0]:+.3f} /"
        f" {shift[1]:+.3f}; run {seconds:.1f} s, total {time.perf_counter() - start:.1f} s",
        flush=True,
    )
    failed = [f"jump {jump[e]:+.3f}" for e in ends if abs(jump[e]) > JUMP_GATE] + [
        f"two-phase mean {two_phase[e]:+.3f}"
        for e, b in zip(ends, bound)
        if abs(two_phase[e]) > b
    ]
    if failed:
        raise RuntimeError("sawtooth gate: " + ", ".join(failed))


if __name__ == "__main__":
    main()
