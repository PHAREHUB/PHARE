#!/usr/bin/env python3
"""
The dumped FK_rho (the solver's charge density, FullyKineticState::rho) pairs with the E of the
fine step before: on the fine level of T2 (fine_gauss.py's 1D setup, N = 0, s = 0.5), FK_rho
dumped at fine time t(m+1) equals, to the precision of the dump, the rho fine_gauss.py rebuilds
from the domain and level ghost particles dumped at t(m), at x - (dt_f / 2) v. So the Gauss
residual dE/dx - c^2 rho on every fine node, the level's boundary nodes included, is the same
with either rho and the E dumped at t(m).

Why: step 1 of a fine step from t(m) deposits rho from the particles at t(m), at their position
half a step back, and E at t(m) is E^{t(m) - dt_f/2} (solver_fk.hpp:36-39, 46-47); nothing changes
the fine rho until the next deposit, after the dump at t(m+1).

Scope: the finest level only. On a level with a finer one, the dumped FK_rho is the step-6
replay's deposit, from the particles at the dump time t(m+1) (solver_fk.hpp:381-384), not the one
of step 1. With divergence_cleaning = 0, rho is never deposited (fk_electric.hpp:69, 89-90, 106-107)
and FK_rho stays at its initial value: hence s = 0.5 here.

Tolerance: the fields are dumped in single precision (PHARE_DIAG_DOUBLES off), the particles in
double. Per node |FK_rho - rho_particles| <= 4 eps_float |FK_rho| + 1e-12 x the gross charge
(sum of |q| w S over both species, the scale of a missing particle being w / ppc).
"""

import argparse
import os
import sys
from itertools import pairwise

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import fine_gauss as fg  # pV's T5 script: setup, particle charge, E residual

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator

EPS_FLOAT = np.finfo(np.float32).eps


def fk_rho(diag_dir, t, nodes):
    """the dumped FK_rho on the given global primal nodes of level 1, physical values only"""
    import h5py

    values = {}
    with h5py.File(f"{diag_dir}/ions_FK_rho.h5", "r") as f:
        key = next(k for k in f["t"] if abs(float(k) - t) < 1e-9)
        for patch in f["t"][key]["pl1"].values():
            lower, upper = int(patch.attrs["lower"][0]), int(patch.attrs["upper"][0])
            dataset = patch["FK_rho"]
            ghosts = int(dataset.attrs["ghosts"])
            data = dataset[:]
            # the primal nodes of the cells lower..upper
            for node in range(lower, upper + 2):
                values.setdefault(node, float(data[node - lower + ghosts]))
    return np.array([values[j] for j in nodes])


def gross_charge(diag_dir, t, nodes, order):
    """sum over both species of |q| w S at x - (dt_f / 2) v: the scale of the node's charge"""
    dl, dt = fg.DL / fg.RATIO, fg.DT / fg.RATIO
    gross = np.zeros(len(nodes))
    reach = (order + 1) // 2 + 1
    for pop in ("protons", "electrons"):
        for quantity in ("domain", "levelGhost"):
            x, v, w = fg.particles(diag_dir, pop, quantity, t, 1)
            pos = x - 0.5 * dt * v / dl
            for i, node in enumerate(nodes):
                near = np.abs(pos - node) < reach
                gross[i] += np.sum(w[near] * fg.bspline(pos[near] - node, order))
    return gross


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--steps", type=int, default=6)
    parser.add_argument("--ppc", type=int, default=50)
    parser.add_argument("--out", default="phare_outputs/fully_kinetic_fk_rho_pairing")
    args = parser.parse_args()

    run = argparse.Namespace(
        order=args.order,
        filter_passes=0,
        cleaning=0.5,
        steps=args.steps,
        ppc=args.ppc,
        vth_factor_protons=1.0,
        vth_factor_electrons=1.0,
        no_strict=False,
    )
    fg.diag_root = args.out
    times = fg.DT * np.arange(args.steps + 1)

    # fine_gauss.simulate declares its diagnostics, then runs: add FK_rho at the same times
    def simulator_with_fk_rho(sim):
        ph.FluidDiagnostics(quantity="FK_rho", write_timestamps=times)
        return Simulator(sim)

    fg.Simulator = simulator_with_fk_rho
    diag_dir, _, seconds = fg.simulate(run, seed=0)

    lo, hi = fg.RATIO * fg.BOX[0], fg.RATIO * (fg.BOX[1] + 1)
    nodes = np.arange(lo, hi + 1)
    fine_times = np.sort(np.concatenate([times, times[:-1] + fg.DT / fg.RATIO]))
    s0 = fg.scale(diag_dir, args.order)

    worst, worst_boundary, failures = 0.0, 0.0, []
    for t_m, t_next in pairwise(fine_times):
        rho_p = fg.charge(diag_dir, t_m, 1, nodes, args.order, ("domain", "levelGhost"))
        rho_s = fk_rho(diag_dir, t_next, nodes)
        gross = gross_charge(diag_dir, t_m, nodes, args.order)
        diff = np.abs(rho_s - rho_p)
        tolerance = 4 * EPS_FLOAT * np.abs(rho_s) + 1e-12 * gross
        # residuals differ by c^2 (rho_p - rho_s), the same E on both sides
        dF = fg.C_OVER_VA**2 * diff / s0
        worst = max(worst, dF.max())
        worst_boundary = max(worst_boundary, dF[0], dF[-1])
        if np.any(diff > tolerance):
            bad = np.flatnonzero(diff > tolerance)
            failures.append(
                f"t(m) {t_m:.4f}: {bad.size} nodes, first {nodes[bad[0]]}: FK_rho "
                f"{rho_s[bad[0]]:+.9e} particles {rho_p[bad[0]]:+.9e}"
            )

    print(
        f"FK_rho(t(m+1)) vs particle rho(t(m)), order {args.order}, {fine_times.size - 1} pairs, "
        f"{nodes.size} fine nodes: max |dF| / scale {worst:.2e}, at the boundary nodes "
        f"{worst_boundary:.2e}; run {seconds:.1f} s",
        flush=True,
    )
    # the residual itself at the last pair, both ways, for the eye (fine_gauss.py's profile)
    F_p = fg.fine_residual(diag_dir, fine_times[-2], args.order, 0)
    rho_p = fg.charge(
        diag_dir, fine_times[-2], 1, nodes, args.order, ("domain", "levelGhost")
    )
    F_s = F_p + fg.C_OVER_VA**2 * (rho_p - fk_rho(diag_dir, fine_times[-1], nodes))
    print(
        f"boundary node F / scale at t(m) = {fine_times[-2]:.4f}, left / right: particles "
        f"{F_p[0] / s0:+.6f} / {F_p[-1] / s0:+.6f}, FK_rho {F_s[0] / s0:+.6f} / {F_s[-1] / s0:+.6f}",
        flush=True,
    )
    for failure in failures:
        print("MISMATCH", failure, flush=True)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
