#!/usr/bin/env python3
"""Runs every documentation run of the MHD and hybrid pages and rewrites validation/data/*.txt.

    python3 run_doc_runs.py --build <PHARE build dir> [--quick] [--jobs N] [--only PATTERN]
    python3 run_doc_runs.py --build <PHARE build dir> --check [--quick]
    python3 run_doc_runs.py --list-modules

The runs (all serial, one MPI rank, 1D unless said):
  R1  convergence of each reconstruction on a small Alfven wave      mhd_convergence_scan.py   -> mhd_convergence_*.txt
  R2  the shock tube with different reconstructions / Riemann solvers  mhd_shock_solvers.py      -> mhd_shock_*.txt
  R3  order in time, MP5, 64 cells, ten periods                       mhd_convergence_scan.py   -> mhd_time_order_*.txt
  R4  hybrid energy budget, 2D uniform plasma (numerics/hybrid_ohm)  hybrid_energy.py          -> hybrid_energy_*.txt

Every data file starts with a header: PHARE commit (git rev-parse HEAD at run time, and whether src/ and pyphare/ match it), build directory, command, date.
--check does not touch the committed files: it reruns into a temporary folder and compares with them (relative
tolerance per quantity, see TOLERANCES; errors at the single-precision output floor are not compared, so a build that
writes double precision, or another compiler, can be checked); the exit status is non-zero on any drift or failure.
--quick runs only the first grids of each convergence scan (for a smoke test); with --check it compares those rows.
Works from any directory; needs numpy and the PHARE python modules listed by --list-modules.
"""
import argparse
import concurrent.futures
import datetime
import fnmatch
import math
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import time

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
REPO = HERE.parents[3]  # scripts -> validation -> source -> doc -> repo root

# ---- the runs -----------------------------------------------------------------------------------------------
CONV = "mhd_convergence_scan.py"
SHOCK = "mhd_shock_solvers.py"
RUNS = []  # (data file, script, argv, build permutation)
DIVB = "mhd_divb_runs.py"
BUILD_ARG = ["<build>"]  # the --build string as given (no absolute path in the data files)


def conv(rec, lim, t, nx):
    RUNS.append((f"mhd_convergence_{rec}_{lim}_t{t}.txt", CONV,
                 [rec, lim, "OUT", "SSPRK4_5", "5e-4", "1e-6", nx, str(t)], f"1_{rec}_{lim}_Rusanov_false"))


conv("Constant", "None", 1, "50:800")
conv("Linear", "VanLeer", 1, "50:1600")
conv("Linear", "MinMod", 1, "50:1600")
conv("WENO3", "None", 1, "50:1600")
for rec in ("WENO3", "WENOZ", "MP5"):
    conv(rec, "None", 10, "16,24,32,48,64")
for rec, lim, riem in [("Constant", "None", "Rusanov"), ("Constant", "None", "HLL"), ("Constant", "None", "HLLD"),
                       ("Linear", "VanLeer", "Rusanov"), ("Linear", "MinMod", "Rusanov")]:
    RUNS.append((f"mhd_shock_{rec}_{lim}_{riem}.txt", SHOCK, [rec, lim, riem, "OUT"], f"1_{rec}_{lim}_{riem}_false"))
for st, dts in [("Euler", ["0.0001953125", "0.000390625", "0.00078125", "0.0015625", "0.003125", "0.00625", "0.0125"]),
                ("TVDRK2", ["0.0015625", "0.003125", "0.00625", "0.0125"]),
                ("TVDRK3", ["0.0015625", "0.003125", "0.00625", "0.0125"]),
                ("SSPRK4_5", ["0.0015625", "0.003125", "0.00625", "0.0125"])]:
    for dt in dts:
        RUNS.append((f"mhd_time_order_{st}_{dt}.txt", CONV, ["MP5", "None", "OUT", st, dt, "1e-6", "64", "10"],
                     "1_MP5_None_Rusanov_false"))

RUNS.append(("mhd_divb_init.txt", DIVB, ["init", "OUT"], "2_Linear_VanLeer_Rusanov_false"))
RUNS.append(("mhd_divb_amr.txt", DIVB, ["amr", "OUT"], "2_Linear_VanLeer_Rusanov_false"))
ENERGY = "hybrid_energy.py"
RUNS.append(("hybrid_energy_z_o1.txt", ENERGY, ["z", "1", "0.01", "20", "OUT"], "2_1_4"))
RUNS.append(("hybrid_energy_z_o3.txt", ENERGY, ["z", "3", "0.01", "20", "OUT"], "2_3_4"))
RUNS.append(("hybrid_energy_xy_o1_dt0.01.txt", ENERGY, ["xy", "1", "0.01", "2", "OUT"], "2_1_4"))
RUNS.append(("hybrid_energy_xy_o1_dt0.005.txt", ENERGY, ["xy", "1", "0.005", "2", "OUT"], "2_1_4"))

# --check tolerances. The committed files were written with single-precision diagnostics (float32 output, ~1e-13
# resolution on B ~ 1e-6); a CI build may write double precision and use another compiler (gcc, clang, -O3). So:
#  * an error is compared if the new OR the committed value is >= ERR_FLOOR (when both are below it the float32
#    rounding dominates; a regression that raises a floor-level error is caught);
#  * a local order is compared only if the errors it is made from are >= ORDER_FLOOR in the new or the committed run;
#  * no bitwise equality is assumed: relative tolerances, per quantity.
ERR_FLOOR = 1e-12
ORDER_FLOOR = 1e-11
TOLERANCES = {"L1": 0.10, "local_order": 0.15, "rho": 1e-3, "By": 1e-3}   # relative (local_order: absolute)
DIVB_REL = 1e-2          # divergences of order 1e-3..1e-1 (single-precision output noise is ~1e-4 relative or less)
DIVB_ABS_BOUNDS = {3: 5e-5, 4: 1e-6}   # mhd_divb_amr: |div fine - div parent| and flux mismatch must stay below these
SHOCK_ABS = 1e-5  # absolute floor for the shock-tube profiles (single vs double output, compiler)
# hybrid energy: particle noise makes the runs chaotic, so another compiler gives other bits. Checked: the
# physics bounds of the page (ions + B gain below ENERGY_GAIN_MAX % of E0 at every sample, ions cool at
# the last time), and dK at the last time within ENERGY_DK_REL of the committed value.
ENERGY_GAIN_MAX = 1e-4
ENERGY_DK_REL = 0.25


def quick_args(name, argv):
    """--quick: only the first three grids of a convergence scan, nothing else changes."""
    if name.startswith("mhd_convergence_"):
        nx = argv[6]
        if ":" in nx:
            a, b = (int(v) for v in nx.split(":"))
            argv = argv[:6] + [",".join(str(a * 2 ** k) for k in range(3) if a * 2 ** k <= b)] + argv[7:]
        else:
            argv = argv[:6] + [",".join(nx.split(",")[:3])] + argv[7:]
    return argv


# ---- helpers --------------------------------------------------------------------------------------------------
def git_head():
    try:
        return subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def code_state():
    """The code the runs are made on: HEAD, and whether src/ and pyphare/ in the working tree are identical to it."""
    try:
        same = subprocess.run(["git", "-C", str(REPO), "diff", "--quiet", "HEAD", "--", "src", "pyphare"]).returncode == 0
    except Exception:
        same = False
    return f"code {git_head()[:9]} (src/ and pyphare/ identical to it: {'yes' if same else 'NO'})"


def header(build, command):
    return (f"# PHARE {code_state()}, build {build}, serial (1 MPI rank)\n"
            f"# written {datetime.date.today().isoformat()} by doc/source/validation/scripts/run_doc_runs.py\n"
            f"# command: {command}\n")


def one_run(item, build, workdir, quick):
    name, script, argv, perm = item
    d = pathlib.Path(workdir) / name[:-4]
    d.mkdir(parents=True, exist_ok=True)
    out = d / "out.txt"
    argv = quick_args(name, argv) if quick else argv
    argv = [str(out) if a == "OUT" else a for a in argv]
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(build), str(REPO / "pyphare"), env.get("PYTHONPATH", "")])
    cmd = [sys.executable, str(HERE / script)] + argv
    t0 = time.time()
    r = subprocess.run(cmd, cwd=d, env=env, stdout=open(d / "log.txt", "w"), stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if r.returncode != 0 or not out.exists():
        return name, dt, None, f"FAILED (see {d / 'log.txt'})", cmd
    lines = [l for l in out.read_text().splitlines() if l.strip()]
    if name.startswith("mhd_shock_"):  # drop the ghost rows
        lines = [l for l in lines if l.startswith("#") or 0 <= float(l.split()[0]) < 800]
    command = "python3 " + script + " " + " ".join(a if a != "OUT" else "out.txt" for a in argv_display(item, quick))
    return name, dt, header(BUILD_ARG[0], command) + "\n".join(lines) + "\n", None, cmd


def argv_display(item, quick):
    name, script, argv, perm = item
    return quick_args(name, argv) if quick else argv


def numbers(text):
    """rows of a data file: list of (label columns..., float columns...) as lists of tokens."""
    return [l.split() for l in text.splitlines() if l.strip() and not l.startswith("#")]


def compare(name, new, old):
    """Returns a list of drift messages comparing two data files (see the tolerance notes above)."""
    problems = []
    rows_new, rows_old = numbers(new), numbers(old)
    # an energy run must reach every committed sample: a run that stops early (e.g. only t=0)
    # would otherwise pass, since the last row is the one checked for dK
    if name.startswith("hybrid_energy_") and [r[0] for r in rows_new] != [r[0] for r in rows_old]:
        return [f"{name}: sample times {[r[0] for r in rows_new]} do not match the committed "
                f"{[r[0] for r in rows_old]} (incomplete history)"]
    old_by = {r[0]: r for r in rows_old}
    for i, r in enumerate(rows_new):
        o = old_by.get(r[0])
        if o is None:
            problems.append(f"{name}: row {r[0]} not in the committed file")
            continue
        if name.startswith("mhd_divb_"):
            cols = (2, 3) if name == "mhd_divb_init.txt" else (1, 2)
            for c in cols:
                x, y = float(r[c]), float(o[c])
                if abs(x - y) > DIVB_REL * abs(y):
                    problems.append(f"{name} row {r[0]} column {c}: {x:.6g} vs committed {y:.6g}")
            if name == "mhd_divb_amr.txt":
                for c, bound in DIVB_ABS_BOUNDS.items():
                    if float(r[c]) > bound:
                        problems.append(f"{name} row {r[0]} column {c}: {float(r[c]):.3g} above {bound:g}")
            continue
        if name.startswith("hybrid_energy_"):
            values = [float(x) for x in r[1:4]]
            if not all(math.isfinite(x) for x in values):
                problems.append(f"{name} t={r[0]}: non-finite energy {r[1:4]}")
                continue
            gain, dK = values[2], values[0]
            # the bound holds at every sample: the page plots the whole history
            if abs(gain) > ENERGY_GAIN_MAX:
                problems.append(f"{name} t={r[0]}: |gain| {abs(gain):.3g} % above {ENERGY_GAIN_MAX:g}")
            if r is rows_new[-1]:
                if dK > 0:
                    problems.append(f"{name} t={r[0]}: the ions heat (dK = {dK:.3g} %)")
                if abs(dK - float(o[1])) > ENERGY_DK_REL * abs(float(o[1])):
                    problems.append(f"{name} t={r[0]}: dK {dK:.4g} vs committed {float(o[1]):.4g}")
            continue
        if name.startswith("mhd_shock_"):
            for label, c in (("rho", 1), ("By", 2)):
                x, y = float(r[c]), float(o[c])
                if abs(x - y) > max(TOLERANCES[label] * abs(y), SHOCK_ABS):
                    problems.append(f"{name} x={r[0]} {label}: {x:.6g} vs committed {y:.6g}")
            continue
        x, y = float(r[2]), float(o[2])           # L1 error
        if math.isnan(y) or math.isnan(x):
            if math.isnan(x) != math.isnan(y):
                problems.append(f"{name} row {r[0]}: L1 {x} vs committed {y} (NaN differs)")
            continue
        if max(x, y) >= ERR_FLOOR and abs(x - y) > TOLERANCES["L1"] * max(abs(y), 1e-300):
            problems.append(f"{name} row {r[0]}: L1 {x:.4g} vs committed {y:.4g}")
        if len(r) > 3 and len(o) > 3 and i > 0:  # local order, between this row and the previous one
            prev = old_by.get(rows_new[i - 1][0])
            newp = float(rows_new[i - 1][2])
            if prev is not None and max(x, y) >= ORDER_FLOOR and max(newp, float(prev[2])) >= ORDER_FLOOR \
                    and abs(float(r[3]) - float(o[3])) > TOLERANCES["local_order"]:
                problems.append(f"{name} row {r[0]}: local order {r[3]} vs committed {o[3]}")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", help="PHARE build directory (the one holding pybindlibs)")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4, help="runs in parallel (each one is a serial process)")
    ap.add_argument("--only", default="*", help="glob on the data file names, e.g. 'mhd_shock_*'")
    ap.add_argument("--list-modules", action="store_true")
    a = ap.parse_args()
    items = [i for i in RUNS if fnmatch.fnmatch(i[0], a.only)]
    if a.list_modules:
        for m in sorted({i[3] for i in items}):
            print("pybindlibs.cpp_" + m)
        return 0
    if not a.build:
        ap.error("--build is required")
    BUILD_ARG[0] = a.build
    build = pathlib.Path(a.build).resolve()
    if not (build / "pybindlibs").is_dir():
        print(f"error: {build}/pybindlibs not found", file=sys.stderr)
        return 2
    missing = [m for m in sorted({i[3] for i in items})
               if not list((build / "pybindlibs").glob(f"cpp_{m}.*"))]
    if missing:
        print("error: missing PHARE modules in the build:", *missing, sep="\n  ", file=sys.stderr)
        return 2
    workdir = tempfile.mkdtemp(prefix="phare_doc_runs_")
    t0 = time.time()
    results, failures = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
        futs = [ex.submit(one_run, i, build, workdir, a.quick) for i in items]
        for f in concurrent.futures.as_completed(futs):
            name, dt, text, err, cmd = f.result()
            print(f"{dt:7.1f} s  {name}" + (f"  {err}" if err else ""), flush=True)
            if err:
                failures.append(f"{name}: {err}")
            else:
                results.append((name, text))
    problems = []
    for name, text in sorted(results):
        committed = DATA / name
        if a.check:
            if not committed.exists():
                problems.append(f"{name}: no committed file")
            else:
                problems += compare(name, text, committed.read_text())
        else:
            DATA.mkdir(exist_ok=True)
            committed.write_text(text)
    print(f"total wall time {time.time() - t0:.0f} s for {len(items)} runs, {a.jobs} in parallel "
          f"({'check' if a.check else 'files written to ' + str(DATA)})")
    for p in failures + problems:
        print("DRIFT/FAIL:", p)
    if not failures and not problems:
        print("OK" + (" (no drift)" if a.check else ""))
    shutil.rmtree(workdir, ignore_errors=True) if not (failures or problems) else print("work dir kept:", workdir)
    return 1 if (failures or problems) else 0


if __name__ == "__main__":
    sys.exit(main())
