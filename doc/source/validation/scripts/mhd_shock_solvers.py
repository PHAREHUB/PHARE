#!/usr/bin/env python3
"""Documentation run R2 (not part of the test suite).

The shock tube of tests/functional/mhd_shock (same initial state, 800 cells of size 1,
time step 0.2, t = 80, periodic) with the reconstruction, limiter and Riemann solver given
on the command line. Writes x, density, By at t = 80 to a text file.

usage: python3 mhd_shock_solvers.py RECONSTRUCTION LIMITER RIEMANN OUTFILE
"""
import sys
import numpy as np
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy.hierarchy_utils import single_patch_for_LO

ph.NO_GUI()
rec, lim, riem, outfile = sys.argv[1:5]
gamma, BX, ncells, dx, final_time = 5.0 / 3.0, 0.75, 800, 1.0, 80
x0 = ncells * dx / 2
diag_dir = f"phare_outputs/shock_{rec}_{lim}_{riem}"

sim = ph.Simulation(
    smallest_patch_size=15, time_step=0.2, final_time=final_time, cells=(ncells,), dl=(dx,),
    refinement="tagging", max_mhd_level=1, max_nbr_levels=1,
    hyper_resistivity=0.0, resistivity=0.0,
    diag_options={"format": "phareh5", "options": {"dir": diag_dir, "mode": "overwrite"}},
    strict=True, eta=0.0, nu=0.0, gamma=gamma,
    reconstruction=rec, limiter=lim, riemann=riem, mhd_timestepper="TVDRK2", hall=False,
    model_options=["MHDModel"],
)
ph.MHDModel(
    density=lambda x: np.where(x < x0, 1.0, 0.125), vx=lambda x: 0.0, vy=lambda x: 0.0,
    vz=lambda x: 0.0, bx=lambda x: BX, by=lambda x: np.where(x < x0, 1.0, -1.0), bz=lambda x: 0.0,
    p=lambda x: np.where(x < x0, 1.0, 0.1),
)
ph.ElectromagDiagnostics(quantity="B", write_timestamps=[final_time])
ph.MHDDiagnostics(quantity="rho", write_timestamps=[final_time])
Simulator(sim).run().reset()

run = Run(diag_dir)
rho = single_patch_for_LO(run.GetMHDrho(final_time)).levels()[0].patches[0].patch_datas["value"]
by = single_patch_for_LO(run.GetB(final_time, all_primal=False).By).levels()[0].patches[0].patch_datas["By"]
with open(outfile, "w") as f:
    f.write(f"# shock tube, {rec} {lim} {riem}, TVDRK2, dt=0.2, t={final_time}, {ncells} cells\n# x rho By\n")
    xr, xb = np.asarray(rho.x), np.asarray(by.x)
    r, b = np.asarray(rho.dataset[:]), np.asarray(by.dataset[:])
    n = min(len(xr), len(xb))
    for i in range(n):
        f.write(f"{xr[i]:.6g} {r[i]:.8g} {b[i]:.8g}\n")
