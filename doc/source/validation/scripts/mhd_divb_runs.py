#!/usr/bin/env python3
"""Documentation runs R5 and R6 (not part of the test suite): the divergence of B in a 2D MHD run.

    python3 mhd_divb_runs.py init OUTFILE      one level, B sampled from a smooth non-solenoidal-on-the-grid field
    python3 mhd_divb_runs.py amr  OUTFILE      two levels (static refined box), same field, fixed time step

The field is the curl of a Gaussian potential A = exp(-r^2 / 2 s^2), s = 0.1, in the unit square, sampled POINTWISE at the
position of each face by PHARE's initialiser: B_x = -0.1 (y - 1/2) A / s^2, B_y = 0.1 (x - 1/2) A / s^2. It is divergence-free
analytically; its discrete divergence on the grid is not zero. Diagnostics are written in single precision (float32): B
values near 0.6 are resolved to about 3e-8, so a divergence of order 1e-2 is far above that floor, and the differences
reported as 'mismatch' are to be read against it.
"""
import sys
import numpy as np
import pyphare.pharein as ph
from pyphare.core.box import Box2D
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()
mode, outfile = sys.argv[1], sys.argv[2]
S = 0.1
AMP = 0.1


def A(x, y):
    return np.exp(-((x - 0.5) ** 2 + (y - 0.5) ** 2) / (2 * S ** 2))


def bx(x, y):
    return -AMP * (y - 0.5) / S ** 2 * A(x, y)


def by(x, y):
    return AMP * (x - 0.5) / S ** 2 * A(x, y)


def make(n, dt, nsteps, diag, times, boxes):
    kw = dict(refinement_boxes={"L0": {"B0": boxes}}, max_mhd_level=2) if boxes is not None else \
        dict(refinement="tagging", max_mhd_level=1, max_nbr_levels=1)
    sim = ph.Simulation(
        time_step=dt, time_step_nbr=nsteps, cells=(n, n), dl=(1.0 / n, 1.0 / n), smallest_patch_size=6,
        nesting_buffer=1, strict=True, eta=0.0, nu=0.0, gamma=5.0 / 3.0, reconstruction="Linear", limiter="VanLeer",
        riemann="Rusanov", mhd_timestepper="TVDRK2", hall=False, model_options=["MHDModel"],
        diag_options={"format": "phareh5", "options": {"dir": diag, "mode": "overwrite"}}, **kw)
    ph.MHDModel(density=lambda x, y: 1.0 + 0 * x, vx=lambda x, y: 0 * x, vy=lambda x, y: 0 * x, vz=lambda x, y: 0 * x,
                bx=bx, by=by, bz=lambda x, y: 0 * x, p=lambda x, y: 1.0 + 0 * x)
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=times)
    return sim


def level_arrays(run, lvl, t):
    B = run.GetB(t, all_primal=False)
    patches = B.levels(t)[lvl].patches
    assert len(patches) == 1, "these runs have one patch per level"
    p = patches[0]
    bxp, byp = p.patch_datas["Bx"], p.patch_datas["By"]
    return (bxp.dataset[:].astype(float), byp.dataset[:].astype(float), int(bxp.ghosts_nbr[0]),
            tuple(p.box.lower), tuple(p.box.upper), bxp.dl)


def divergence(run, lvl, t):
    bxa, bya, ng, lo, up, dl = level_arrays(run, lvl, t)
    nx, ny = up[0] - lo[0] + 1, up[1] - lo[1] + 1
    Bx, By = bxa[ng:ng + nx + 1, ng:ng + ny], bya[ng:ng + nx, ng:ng + ny + 1]
    return (Bx[1:] - Bx[:-1]) / dl[0] + (By[:, 1:] - By[:, :-1]) / dl[1], lo


if mode == "init":
    rows = []
    for n in (16, 32, 64, 128):
        dt, nsteps = 0.2 / n, 10
        ph.global_vars.sim = None
        sim = make(n, dt, nsteps, f"o_init_{n}", [0.0, nsteps * dt], None)
        Simulator(sim).run().reset()
        run = Run(f"o_init_{n}")
        d0, _ = divergence(run, 0, 0.0)
        d1, _ = divergence(run, 0, nsteps * dt)
        rows.append((n, 1.0 / n, abs(d0).max(), abs(d1).max()))
        print(rows[-1], flush=True)
    with open(outfile, "w") as f:
        f.write("# R5: one level, 2D MHD, Gaussian-curl field sampled at the faces; max |div B| over the interior cells\n")
        f.write("# n dx maxdiv_t0 maxdiv_after_10_steps\n")
        for r in rows:
            f.write(f"{r[0]} {r[1]:.8g} {r[2]:.6e} {r[3]:.6e}\n")
else:
    n, dt, nsteps = 32, 2e-3, 10
    times = [0.0] + [2 * dt * k for k in range(1, 6)]
    sim = make(n, dt, nsteps, "o_amr", times, Box2D(10, 21))
    Simulator(sim).run().reset()
    run = Run("o_amr")
    rows = []
    for t in times:
        d0, l0 = divergence(run, 0, t)
        d1, l1 = divergence(run, 1, t)
        # fine cell versus its parent coarse cell (inside the refined box)
        diff = max(abs(d1[fi, fj] - d0[(fi + l1[0]) // 2 - l0[0], (fj + l1[1]) // 2 - l0[1]])
                   for fi in range(d1.shape[0]) for fj in range(d1.shape[1]))
        bx0, by0, ng0, lo0, _, _ = level_arrays(run, 0, t)
        bx1, by1, ng1, lo1, _, _ = level_arrays(run, 1, t)
        mmx = max(abs(bx0[I - lo0[0] + ng0, J - lo0[1] + ng0]
                      - 0.5 * (bx1[2 * I - lo1[0] + ng1, 2 * J - lo1[1] + ng1] + bx1[2 * I - lo1[0] + ng1, 2 * J + 1 - lo1[1] + ng1]))
                  for I in range(10, 23) for J in range(10, 22))
        mmy = max(abs(by0[I - lo0[0] + ng0, J - lo0[1] + ng0]
                      - 0.5 * (by1[2 * I - lo1[0] + ng1, 2 * J - lo1[1] + ng1] + by1[2 * I + 1 - lo1[0] + ng1, 2 * J - lo1[1] + ng1]))
                  for I in range(10, 22) for J in range(10, 23))
        rows.append((t, abs(d0).max(), abs(d1).max(), diff, max(mmx, mmy)))
        print(rows[-1], flush=True)
    bxa = level_arrays(run, 0, times[0])[0]
    bxb = level_arrays(run, 0, times[-1])[0]
    change = abs(bxb - bxa).max()
    with open(outfile, "w") as f:
        f.write("# R6: two levels (static refined box, coarse cells 10..21), 2D MHD, same field, fixed time step 2e-3 (the coarse step)\n")
        f.write("# every time listed is the end of a coarse step (levels synchronised)\n")
        f.write("# t maxdiv_level0 maxdiv_level1 max|div_fine_cell - div_parent_cell| max|flux_coarse - flux_fine_sum_per_face|\n")
        for r in rows:
            f.write(f"{r[0]:.6g} {r[1]:.6e} {r[2]:.6e} {r[3]:.6e} {r[4]:.6e}\n")
        f.write(f"# max change of the coarse B_x over the run: {change:.6e}\n")
