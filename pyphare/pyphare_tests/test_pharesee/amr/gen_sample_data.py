#!/usr/bin/env python3
"""Write the sample runs of the pharesee.amr tests: each simulation in phareh5 AND vtkhdf.

usage:
  gen_sample_data.py [--root DIR] [--force] all            every single-process case, one
                                                           subprocess per (case, format)
  gen_sample_data.py [--root DIR] [--force] CASE FORMAT    one run in this process; under
                                                           ``mpirun -n 3`` with --mpi3
  gen_sample_data.py [--root DIR] [--force] --mpi3 FORMAT  the 3-rank run (2d_tagged), to be
                                                           launched with ``mpirun -n 3``
CASE: 1d, 2d_static, 2d_tagged, 3d; FORMAT: h5, vtk.
Output: ``<root>/<case>/<format>`` (single process) and ``<root>_mpi3/2d_tagged/<format>``
(3 ranks). root: --root, else ``$PHARESEE_SAMPLE_DATA``, else ``<repo>/build/sample_data``
(see sample_data.py). PHARE's per-rank logs go to ``<case dir>/.log``.
A run whose stamp ``<case dir>/.<format>.done`` matches this file's content is skipped;
--force reruns it. Needs PYTHONPATH=<build dir>:<repo>/pyphare (as every python test).

The expected values of the tests depend on these simulations: change a parameter and the
data, and the tests, change.
"""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time

CASES = ("2d_static", "2d_tagged", "1d", "3d")
FORMATS = {"h5": "phareh5", "vtk": "pharevtkhdf"}
MPI3_CASE = "2d_tagged"


def _default_root():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from sample_data import sample_root

    return sample_root()


def _fingerprint():
    with open(os.path.abspath(__file__), "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def harris_funcs(L=0.5):
    import numpy as np

    def mk(sim):
        Ly = sim.simulation_domain()[1]
        Lx = sim.simulation_domain()[0]

        def density(x, y):
            return (
                0.4
                + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
                + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
            )

        def S(y, y0, l):
            return 0.5 * (1.0 + np.tanh((y - y0) / l))

        def by(x, y):
            sigma, dB = 1.0, 0.1
            x0 = x - 0.5 * Lx
            y1 = y - 0.3 * Ly
            y2 = y - 0.7 * Ly
            return 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / sigma**2) - 2 * dB * x0 * np.exp(
                -(x0**2 + y2**2) / sigma**2
            )

        def bx(x, y):
            sigma, dB = 1.0, 0.1
            x0 = x - 0.5 * Lx
            y1 = y - 0.3 * Ly
            y2 = y - 0.7 * Ly
            dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / sigma**2)
            dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / sigma**2)
            return -1 + 2.0 * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

        def bz(x, y):
            return 0.0

        def T(x, y):
            b2 = bx(x, y) ** 2 + by(x, y) ** 2
            return 1.0 / density(x, y) * (0.7 - b2 * 0.5)

        return density, bx, by, bz, T

    return mk


def build(case, fmt, out_dir):
    import numpy as np

    import pyphare.pharein as ph

    ph.NO_GUI()
    diag = {"format": FORMATS[fmt], "options": {"dir": out_dir, "mode": "overwrite"}}
    common = dict(interp_order=1, diag_options=diag, strict=True, nesting_buffer=1)

    def zero(*a):
        return 0.0 * a[0] if hasattr(a[0], "shape") else 0.0

    vx = zero  # bulk velocity along x: zero except in 2d_tagged
    if case == "2d_static":
        dt, nsteps = 0.001, 20
        sim = ph.Simulation(
            time_step=dt, time_step_nbr=nsteps, cells=(100, 80), dl=(0.4, 0.5),
            refinement_boxes={
                "L0": {"B0": [(8, 8), (30, 28)], "B1": [(55, 25), (85, 52)]},
                "L1": {"B0": [(24, 24), (50, 50)]},
            },
            **common,
        )
        ts = dt * np.array([0, 7, 14, 20])
        density, bx, by, bz, T = harris_funcs()(sim)
        ppc = 25
    elif case == "2d_tagged":
        dt, nsteps = 0.005, 40
        sim = ph.Simulation(
            time_step=dt, time_step_nbr=nsteps, cells=(80, 80), dl=(0.5, 0.5),
            refinement="tagging", max_nbr_levels=3, tagging_threshold=0.1,
            hyper_resistivity=0.002, resistivity=0.001,
            **common,
        )
        ts = dt * np.linspace(0, nsteps, 5).astype(int)
        density0, bx, by, bz, T = harris_funcs()(sim)
        Lx, Ly = sim.simulation_domain()
        VD = 8.0

        def blob(x, y):  # localized blob, drifting in +x
            return np.exp(-((x - 0.25 * Lx) ** 2 + (y - 0.3 * Ly) ** 2) / 2.0)

        def density(x, y):
            return density0(x, y) + 1.5 * blob(x, y)

        def vx(x, y):
            return VD * blob(x, y)

        ppc = 10
    elif case == "1d":
        dt, nsteps = 0.001, 30
        sim = ph.Simulation(
            time_step=dt, time_step_nbr=nsteps, cells=200, dl=0.2,
            refinement_boxes={"L0": {"B0": [(40,), (110,)]}},
            **common,
        )
        ts = dt * np.array([0, 15, 30])
        L = sim.simulation_domain()[0]

        def density(x):
            return 0.6 + 0.3 * np.cos(2 * np.pi * x / L)

        def bx(x):
            return 1.0

        def by(x):
            return 0.3 * np.sin(2 * np.pi * x / L)

        def bz(x):
            return 0.2 * np.cos(4 * np.pi * x / L)

        def T(x):
            return 0.3 + 0 * x

        ppc = 100
    elif case == "3d":
        dt, nsteps = 0.001, 6
        sim = ph.Simulation(
            time_step=dt, time_step_nbr=nsteps, cells=(20, 20, 20), dl=(0.4, 0.4, 0.4),
            refinement_boxes={"L0": {"B0": [(5, 5, 5), (14, 14, 14)]}},
            **common,
        )
        ts = dt * np.array([0, 6])
        L = sim.simulation_domain()

        def density(x, y, z):
            return 0.6 + 0.2 * np.cos(2 * np.pi * x / L[0]) * np.cos(
                2 * np.pi * y / L[1]) * np.cos(2 * np.pi * z / L[2])

        def bx(x, y, z):
            return 1.0

        def by(x, y, z):
            return 0.3 * np.sin(2 * np.pi * x / L[0])

        def bz(x, y, z):
            return 0.3 * np.cos(2 * np.pi * y / L[1])

        def T(x, y, z):
            return 0.3 + 0 * x

        ppc = 8
    else:
        raise ValueError(f"unknown case {case!r}; known: {CASES}")

    def vth(*a):
        return np.sqrt(T(*a))

    vvv = {"vbulkx": vx, "vbulky": zero, "vbulkz": zero,
           "vthx": vth, "vthy": vth, "vthz": vth}
    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"mass": 1, "charge": 1, "density": density, **vvv,
                 "nbr_part_per_cell": ppc, "init": {"seed": 1337}},
        alpha={"mass": 4, "charge": 2, "density": density, **vvv,
               "nbr_part_per_cell": ppc, "init": {"seed": 2334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)

    ts = [float(t) for t in ts]
    for q in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=q, write_timestamps=ts)
    for q in ["charge_density", "mass_density", "bulkVelocity", "charge_flux", "momentum_tensor"]:
        ph.FluidDiagnostics(quantity=q, write_timestamps=ts)
    for pop in ["protons", "alpha"]:
        for q in ["density", "flux", "momentum_tensor"]:
            ph.FluidDiagnostics(quantity=q, write_timestamps=ts, population_name=pop)
    ph.MetaDiagnostics(quantity="tags", write_timestamps=ts)
    return sim


def _stamp(case_dir, fmt):
    return os.path.join(case_dir, f".{fmt}.done")


def _is_done(case_dir, fmt):
    try:
        with open(_stamp(case_dir, fmt)) as f:
            return f.read().strip() == _fingerprint()
    except OSError:
        return False


def run_one(root, case, fmt, force):
    """One simulation in this process (every rank of an mpirun calls it)."""
    from pyphare import cpp
    from pyphare.simulator.simulator import Simulator, startMPI

    case_dir = os.path.join(root, case)
    out = os.path.join(case_dir, fmt)
    if not force and _is_done(case_dir, fmt):
        print(f"up to date: {out}", flush=True)
        return
    startMPI()
    if cpp.mpi_rank() == 0:
        os.makedirs(case_dir, exist_ok=True)
        if os.path.exists(_stamp(case_dir, fmt)):
            os.remove(_stamp(case_dir, fmt))
        shutil.rmtree(out, ignore_errors=True)  # no stale file of an earlier run
    cpp.mpi_barrier()
    os.chdir(case_dir)  # PHARE's .log directory lands in the case directory
    sim = build(case, fmt, out)
    t0 = time.time()
    Simulator(sim).run().reset()
    cpp.mpi_barrier()
    if cpp.mpi_rank() == 0:
        with open(_stamp(case_dir, fmt), "w") as f:
            f.write(_fingerprint() + "\n")
        print(f"DONE {out} {time.time() - t0:.1f}s", flush=True)


def run_all(root, force):
    """Every single-process case, each run in its own process (one Simulator per process)."""
    failed = []
    for case in CASES:
        for fmt in FORMATS:
            cmd = [sys.executable, "-u", os.path.abspath(__file__), "--root", root, case, fmt]
            if force:
                cmd.append("--force")
            if subprocess.run(cmd, check=False).returncode != 0:
                failed.append(f"{case}/{fmt}")
    if failed:
        sys.exit(f"gen_sample_data: failed runs: {', '.join(failed)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("what", nargs="+", help="'all' | CASE FORMAT | FORMAT with --mpi3")
    parser.add_argument("--root", default=None, help="output root (default: see module doc)")
    parser.add_argument("--mpi3", action="store_true",
                        help="the 3-rank 2d_tagged run, under mpirun -n 3, into <root>_mpi3")
    parser.add_argument("--force", action="store_true", help="rerun up-to-date runs")
    args = parser.parse_args()
    root = os.path.abspath(args.root or _default_root())

    if args.mpi3:
        if len(args.what) != 1 or args.what[0] not in FORMATS:
            parser.error("--mpi3 takes one FORMAT (h5 or vtk)")
        run_one(root.rstrip(os.sep) + "_mpi3", MPI3_CASE, args.what[0], args.force)
    elif args.what == ["all"]:
        run_all(root, args.force)
    elif len(args.what) == 2 and args.what[0] in CASES and args.what[1] in FORMATS:
        run_one(root, args.what[0], args.what[1], args.force)
    else:
        parser.error(f"expected 'all' or CASE FORMAT with CASE in {CASES}, FORMAT in "
                     f"{tuple(FORMATS)}")


if __name__ == "__main__":
    main()
