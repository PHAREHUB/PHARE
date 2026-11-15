#!/usr/bin/env python3
"""Documentation run R1 (not part of the test suite).

Same set-up as tests/functional/mhd_convergence/convergence.py (small Alfven wave,
Rusanov, SSPRK4_5, dt = 5e-4, t = 1), but: one reconstruction per call, errors are
printed and written to a file, no assertion.

usage: python3 mhd_convergence_scan.py RECONSTRUCTION LIMITER OUTFILE TIMESTEPPER DT AMPLITUDE NX_LIST FINAL_TIME
NX_LIST: comma-separated cell numbers, or FIRST:LAST (nx doubles). The test itself uses 1e-6, 50..800, t=1.
Diagnostics are written in single precision (float32): errors below ~1e-13 cannot be resolved.
"""
import sys
import numpy as np
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy.hierarchy_utils import single_patch_for_LO

ph.NO_GUI()
reconstruction = sys.argv[1]
limiter = sys.argv[2] if len(sys.argv) > 2 else "None"
outfile = sys.argv[3] if len(sys.argv) > 3 else f"conv_{reconstruction}.txt"
timestepper = sys.argv[4] if len(sys.argv) > 4 else "SSPRK4_5"
time_step = float(sys.argv[5]) if len(sys.argv) > 5 else 5e-4
amp = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-6
nx_arg = sys.argv[7] if len(sys.argv) > 7 else "50:1600"
if ":" in nx_arg:
    a, b = (int(v) for v in nx_arg.split(":"))
    nx_list = [a * 2**k for k in range(64) if a * 2**k <= b]
else:
    nx_list = [int(v) for v in nx_arg.split(",")]
final_time = float(sys.argv[8]) if len(sys.argv) > 8 else 1.0
diag_dir = f"phare_outputs/conv_{reconstruction}_{limiter}_{timestepper}_{time_step}"


def config(nx, dx):
    sim = ph.Simulation(
        smallest_patch_size=15, time_step=time_step, final_time=final_time,
        cells=(nx,), dl=(dx,), refinement="tagging", max_mhd_level=1, max_nbr_levels=1,
        hyper_resistivity=0.0, resistivity=0.0,
        diag_options={"format": "phareh5", "options": {"dir": diag_dir, "mode": "overwrite"}},
        strict=True, nesting_buffer=1, eta=0.0, nu=0.0, gamma=5.0 / 3.0,
        reconstruction=reconstruction, limiter=limiter, riemann="Rusanov",
        mhd_timestepper=timestepper, model_options=["MHDModel"],
    )
    ph.MHDModel(
        density=lambda x: 1.0, vx=lambda x: 0.0,
        vy=lambda x: -amp * np.cos(2 * np.pi * x), vz=lambda x: 0.0,
        bx=lambda x: 1.0, by=lambda x: amp * np.cos(2 * np.pi * x), bz=lambda x: 0.0,
        p=lambda x: 0.1,
    )
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=[0.0, final_time])
    return sim


def by_at(run, t):
    return (single_patch_for_LO(run.GetB(t, all_primal=False).By)
            .levels()[0].patches[0].patch_datas["By"].dataset[:])


rows = []
for nx in nx_list:
    ph.global_vars.sim = None
    sim = config(nx, 1.0 / nx)
    Simulator(sim).run().reset()
    run = Run(diag_dir)
    a, b = by_at(run, final_time), by_at(run, 0.0)
    err = float(np.sum(np.abs(a - b)) / len(a))     # same L1 as the test
    print(f"nx={nx} dx={1.0/nx:.6g} L1={err:.6e} nan={bool(np.isnan(a).any())}", flush=True)
    rows.append((nx, 1.0 / nx, err))
    if not np.isfinite(err):
        break

# local order between consecutive grids (never a fit through points near the float32 floor)
rates = [None]
for (n0, d0, e0), (n1, d1, e1) in zip(rows[:-1], rows[1:]):
    rates.append(np.log(e0 / e1) / np.log(d0 / d1) if e0 > 0 and e1 > 0 else None)
with open(outfile, "w") as f:
    f.write(f"# {reconstruction} {limiter} {timestepper} dt={time_step} amp={amp} final_time={final_time}\n")
    f.write("# nx dx L1 local_order\n")
    for r, q in zip(rows, rates):
        f.write(f"{r[0]} {r[1]:.8g} {r[2]:.8e} {'' if q is None else f'{q:.3f}'}\n")
print(open(outfile).read())
