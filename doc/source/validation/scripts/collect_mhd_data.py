#!/usr/bin/env python3
"""Collects the outputs of the documentation runs into data/*.txt, each with a header
(commit, build, command, date). Nothing is typed by hand.

usage: python3 collect_mhd_data.py RUNS_DIR
RUNS_DIR holds r1d_<REC>_<LIM>_<T>/conv.txt (mhd_convergence_scan.py), r2/<REC>_<LIM>_<RIEMANN>/shock.txt
(mhd_shock_solvers.py) and r3/<STEPPER>_<DT>/conv.txt (mhd_convergence_scan.py MP5 ... 64 10).
"""
import datetime
import pathlib
import sys

runs = pathlib.Path(sys.argv[1])
out = pathlib.Path(__file__).resolve().parent.parent / "data"
out.mkdir(exist_ok=True)
head = ("# PHARE commit dc62fdd4e (documentation branch), build build/rel, serial, 1 rank\n"
        "# collected {date} by doc/source/validation/scripts/collect_mhd_data.py\n")
date = datetime.date.today().isoformat()


def write(name, src, command):
    lines = [l for l in pathlib.Path(src).read_text().splitlines() if l.strip()]
    (out / name).write_text(head.format(date=date) + f"# command: {command}\n" + "\n".join(lines) + "\n")


for d in sorted(runs.glob("r1d_*")):
    rec, lim, t = d.name[4:].split("_")
    write(f"mhd_convergence_{rec}_{lim}_t{t}.txt", d / "conv.txt",
          f"python3 mhd_convergence_scan.py {rec} {lim} conv.txt SSPRK4_5 5e-4 1e-6 <nx list> {t}")
for d in sorted((runs / "r3").glob("*")):
    st, dt = d.name.rsplit("_", 1)
    write(f"mhd_time_order_{st}_{dt}.txt", d / "conv.txt",
          f"python3 mhd_convergence_scan.py MP5 None conv.txt {st} {dt} 1e-6 64 10")
for d in sorted((runs / "r2").glob("*")):
    rec, lim, riem = d.name.split("_")
    src = d / "shock.txt"
    lines = [l for l in src.read_text().splitlines() if l.startswith("#") or (l.strip() and 0 <= float(l.split()[0]) < 800)]
    (out / f"mhd_shock_{rec}_{lim}_{riem}.txt").write_text(
        head.format(date=date) + f"# command: python3 mhd_shock_solvers.py {rec} {lim} {riem} shock.txt (ghost rows removed)\n"
        + "\n".join(lines) + "\n")
