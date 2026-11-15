#!/usr/bin/env python3
"""Error ratios for each halving of the time step (documentation run R3), from validation/data/mhd_time_order_*.txt.
For TVDRK3 the space error is estimated by the SSPRK4_5 error at the same time step and subtracted."""
import glob, os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
data = {}
for f in glob.glob(os.path.join(here, "..", "data", "mhd_time_order_*.txt")):
    name = os.path.basename(f)[len("mhd_time_order_"):-4]
    st, dt = name.rsplit("_", 1)
    row = [l.split() for l in open(f) if l[0] != "#" and l.split()[0] == "64"][0]
    data.setdefault(st, {})[float(dt)] = float(row[2])
for st in ("Euler", "TVDRK2", "TVDRK3", "SSPRK4_5"):
    dts = sorted(data[st])
    print(st, {dt: data[st][dt] for dt in dts})
    ok = [dt for dt in dts if np.isfinite(data[st][dt])]
    print("  ratio per halving of dt:", [round(data[st][a] / data[st][b], 2) for a, b in zip(ok[1:], ok[:-1])])
space = data["SSPRK4_5"]
ok = [dt for dt in sorted(data["TVDRK3"]) if np.isfinite(data["TVDRK3"][dt])]
sub = {dt: data["TVDRK3"][dt] - space[dt] for dt in ok}
print("TVDRK3 minus SSPRK4_5 error, ratio per halving:", [round(sub[a] / sub[b], 2) for a, b in zip(ok[1:], ok[:-1])])
