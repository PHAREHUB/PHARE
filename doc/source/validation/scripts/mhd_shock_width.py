#!/usr/bin/env python3
"""Widths, in cells, of two waves of the shock-tube runs R2 (validation/data/mhd_shock_*.txt).
The two plateau values of each wave are taken once, from the Linear + VanLeer run (the sharpest), and used for every run:
10 % and 90 % of the step between them; count the cells whose density lies between them.
 - contact discontinuity near x = 452: plateaus = median density over x = 430..440 and 465..475, counted over x = 400..500
 - slow shock near x = 505: plateaus = median density over x = 465..475 and 545..555, counted over x = 490..540
A documentation measure, not a standard definition; it only serves to compare the runs with each other."""
import glob, os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
data = {}
for f in sorted(glob.glob(os.path.join(here, "..", "data", "mhd_shock_*.txt"))):
    data[os.path.basename(f)[len("mhd_shock_"):-4]] = np.loadtxt(f, comments="#").T
ref = data["Linear_VanLeer_Rusanov"]


def med(x, rho, a, b):
    return float(np.median(rho[(x >= a) & (x <= b)]))


waves = {"contact": ((430, 440), (465, 475), (400, 500)), "slow shock": ((465, 475), (545, 555), (490, 540))}
for wave, (wb, wf, wc) in waves.items():
    hi, lo = med(ref[0], ref[1], *wb), med(ref[0], ref[1], *wf)
    print(f"{wave}: plateaus {hi:.3f} -> {lo:.3f}; By just before/after (ref run): ", end="")
    x, rho, by = ref
    print(f"{med(x, by, *wb):.3f} / {med(x, by, *wf):.3f}")
    for name, (x, rho, by) in data.items():
        sel = (x >= wc[0]) & (x <= wc[1])
        n = int(np.sum((rho[sel] < lo + 0.9 * (hi - lo)) & (rho[sel] > lo + 0.1 * (hi - lo))))
        print(f"   {name:26s} {n:3d} cells")
