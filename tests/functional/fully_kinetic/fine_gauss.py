#!/usr/bin/env python3
"""
T5 (P4): the Gauss residual F = div E - c^2 F^N rho on a fine level next to its boundary, 1D, the
plasma of T2 (tests/simulator/test_fully_kinetic_amr.py: 64 root cells, uniform thermal ions and
electrons T = 0.03, m_i / m_e = 25, a small B and electron velocity wave, box over root cells 24 to
39), measured as a profile F(x, t) against the distance x to the nearer level boundary, every
coarse step. Model and hypotheses: docs/superpowers/specs/p4/t5_fine_gauss.md (t5_model.py).

F on fine primal node j at a coarse time t: the dumped fine E (E^{t - dt_f/2}, ghosts included, so
the boundary node uses the ghost the solver filled from the coarse E), and rho at x - (dt_f/2) v
from the domain AND the level ghost particles of the fine level (the band deposits since D5),
filtered N times as the solver does. Relative to the scale of T2: c^2 max |rho_root| at t = 0.

Prints, per run, F / scale at distances 0, 3, 6, 10 fine cells from each boundary at the last
step and at step 20, and saves the whole profile history (fine_gauss_<options>_seed<k>.npz).

The dump phase: a coarse time comes right after the band renewal, while the fine level is also
dumped one fine step later (mid step, t + dt_f). The boundary node then differs by the charge of
the domain particles erased in that fine step. The npz also holds the mid step profiles and, from
step JUMP_FROM on, the jump (coarse time mean - mid step mean) and the two-phase mean.
"""

import argparse
import time
from math import comb, factorial

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

CELLS, DL, DT = 64, 0.04, 0.002
C_OVER_VA, MASS_E, T = 10.0, 0.04, 0.03
BOX = (24, 39)  # root cells, both included
RATIO = 2
PROBES = (0, 3, 6, 10)
JUMP_FROM = 10  # coarse step, past the start-up
diag_root = "phare_outputs/fully_kinetic_fine_gauss"


def bspline(u, order):
    h = (order + 1) / 2
    return sum(
        (-1) ** k * comb(order + 1, k) * np.maximum(u + h - k, 0) ** order
        for k in range(order + 2)
    ) / factorial(order)


def uniform(value):
    return lambda x: value + 0.0 * x


def wave(amplitude):
    return lambda x: amplitude * np.sin(2 * np.pi * x / (CELLS * DL))


def simulate(args, seed):
    tag = f"o{args.order}_N{args.filter_passes}_s{args.cleaning:g}_ppc{args.ppc}_steps{args.steps}"
    if (args.vth_factor_protons, args.vth_factor_electrons) != (1.0, 1.0):
        tag += f"_vi{args.vth_factor_protons:g}_ve{args.vth_factor_electrons:g}"
    diag_dir = f"{diag_root}/{tag}_seed{seed}"
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=DT,
        time_step_nbr=args.steps,
        cells=[CELLS],
        dl=[DL],
        interp_order=args.order,
        model_options=["FullyKineticModel"],
        c_over_va=C_OVER_VA,
        current_filter_passes=args.filter_passes,
        divergence_cleaning=args.cleaning,
        refinement_boxes={"L0": {"B0": [(BOX[0],), (BOX[1],)]}},
        smallest_patch_size=8,
        largest_patch_size=16,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite", "fine_dump_lvl_max": 1},
        },
        strict=not args.no_strict,
    )
    common = {"density": uniform(1.0), "nbr_part_per_cell": args.ppc}
    vth = {
        "protons": args.vth_factor_protons * np.sqrt(T),
        "electrons": args.vth_factor_electrons * np.sqrt(T / MASS_E),
    }
    ph.MaxwellianFluidModel(
        bx=uniform(1.0),
        by=wave(0.1),
        bz=uniform(0.0),
        protons={
            "charge": 1,
            "mass": 1.0,
            **{f"vbulk{k}": uniform(0.0) for k in "xyz"},
            **{f"vth{k}": uniform(vth["protons"]) for k in "xyz"},
            "init": {"seed": 2026 + 10 * seed},
            **common,
        },
        electrons={
            "charge": -1,
            "mass": MASS_E,
            "vbulkx": wave(0.05),
            "vbulky": uniform(0.0),
            "vbulkz": uniform(0.0),
            **{f"vth{k}": uniform(vth["electrons"]) for k in "xyz"},
            "init": {"seed": 2027 + 10 * seed},
            **common,
        },
    )
    times = DT * np.arange(args.steps + 1)
    ph.ElectromagDiagnostics(quantity="E", write_timestamps=times)
    for pop in ("protons", "electrons"):
        for quantity in ("domain", "levelGhost"):
            ph.ParticleDiagnostics(
                quantity=quantity, write_timestamps=times, population_name=pop
            )
    start = time.perf_counter()
    Simulator(sim).run().reset()
    return diag_dir, times, time.perf_counter() - start


def particles(diag_dir, pop, quantity, t, level):
    """x in cells of the level (global index + delta), v_x, weight"""
    import h5py

    xs, vs, ws = [], [], []
    with h5py.File(f"{diag_dir}/ions_pop_{pop}_{quantity}.h5", "r") as f:
        key = next(k for k in f["t"] if abs(float(k) - t) < 1e-9)
        group = f["t"][key].get(f"pl{level}")
        for patch in group.values() if group is not None else []:
            if "iCell" not in patch:
                continue
            xs.append(patch["iCell"][:, 0] + patch["delta"][:, 0])
            vs.append(patch["v"][:, 0])
            ws.append(patch["weight"][:].ravel())
    if not xs:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    return np.concatenate(xs), np.concatenate(vs), np.concatenate(ws)


def charge(diag_dir, t, level, nodes, order, quantities):
    """rho at x - (dt_l / 2) v on the given global primal nodes of a level, unfiltered"""
    dl, dt = DL / RATIO**level, DT / RATIO**level
    rho = np.zeros(len(nodes))
    reach = (order + 1) // 2 + 1
    for pop, q in (("protons", 1.0), ("electrons", -1.0)):
        for quantity in quantities:
            x, v, w = particles(diag_dir, pop, quantity, t, level)
            pos = x - 0.5 * dt * v / dl
            for i, node in enumerate(nodes):
                near = np.abs(pos - node) < reach
                rho[i] += q * np.sum(w[near] * bspline(pos[near] - node, order))
    return rho


def scale(diag_dir, order):
    """c^2 max |rho_root| at t = 0 (T2's scale), the root periodic"""
    nodes = np.arange(CELLS)
    rho = np.zeros(CELLS)
    for pop, q in (("protons", 1.0), ("electrons", -1.0)):
        x, v, w = particles(diag_dir, pop, "domain", 0.0, 0)
        pos = x - 0.5 * DT * v / DL
        for i in nodes:
            u = (pos - i + CELLS / 2) % CELLS - CELLS / 2
            near = np.abs(u) < (order + 1) // 2 + 1
            rho[i] += q * np.sum(w[near] * bspline(u[near], order))
    return C_OVER_VA**2 * np.abs(rho).max()


def fine_residual(diag_dir, t, order, passes):
    """F on the fine primal nodes of the box, both boundary nodes included, in global indices"""
    dl = DL / RATIO
    lo, hi = RATIO * BOX[0], RATIO * (BOX[1] + 1)
    nodes = np.arange(lo - passes, hi + passes + 1)
    rho = charge(diag_dir, t, 1, nodes, order, ("domain", "levelGhost"))
    for _ in range(passes):
        rho[1:-1] = 0.25 * rho[:-2] + 0.5 * rho[1:-1] + 0.25 * rho[2:]
    rho = rho[passes : len(rho) - passes]
    Ex = {}
    hier = hierarchy_from(h5_filename=f"{diag_dir}/EM_E.h5", times=f"{t:.10f}")
    for patch in hier.level(1, t).patches:
        pdata = patch.patch_datas["Ex"]
        for x, value in zip(pdata.x, pdata.dataset[:]):
            key = round(x / dl - 0.5)
            inside = lo <= key < hi
            if inside or key not in Ex:  # a ghost never overwrites a physical value
                Ex[key] = float(value)
    E = np.array([Ex[j] for j in range(lo - 1, hi + 1)])
    return np.diff(E) / dl - C_OVER_VA**2 * rho


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--filter-passes", type=int, default=0)
    parser.add_argument("--cleaning", type=float, default=0.5)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--ppc", type=int, default=50)
    parser.add_argument("--seeds", type=int, default=1)
    parser.add_argument("--first-seed", type=int, default=0)
    parser.add_argument("--vth-factor-protons", type=float, default=1.0)
    parser.add_argument("--vth-factor-electrons", type=float, default=1.0)
    # setup checks as warnings, e.g. electrons x2 past the non relativistic pusher's 0.3 c
    parser.add_argument("--no-strict", action="store_true")
    args = parser.parse_args()

    def history(diag_dir, times, s0):
        return (
            np.array(
                [
                    fine_residual(diag_dir, t, args.order, args.filter_passes)
                    for t in times
                ]
            )
            / s0
        )

    for seed in range(args.first_seed, args.first_seed + args.seeds):
        diag_dir, times, seconds = simulate(args, seed)
        s0 = scale(diag_dir, args.order)
        profiles = history(diag_dir, times, s0)
        mid_times = times[:-1] + DT / RATIO
        mid_profiles = history(diag_dir, mid_times, s0)
        coarse_mean = profiles[JUMP_FROM:].mean(axis=0)
        mid_mean = mid_profiles[JUMP_FROM:].mean(axis=0)
        jump, two_phase = coarse_mean - mid_mean, 0.5 * (coarse_mean + mid_mean)
        np.savez(
            f"{diag_dir}.npz",
            times=times,
            profiles=profiles,
            mid_times=mid_times,
            mid_profiles=mid_profiles,
            jump=jump,
            two_phase=two_phase,
            scale=s0,
            args=str(vars(args)),
        )
        print(
            f"seed {seed} steps {JUMP_FROM}-{args.steps}: jump (coarse - mid) at the boundary "
            f"node {jump[0]:+.3f} / {jump[-1]:+.3f}, one cell in {jump[1]:+.3f} / {jump[-2]:+.3f}; "
            f"two-phase mean {two_phase[0]:+.3f} / {two_phase[-1]:+.3f}, one cell in "
            f"{two_phase[1]:+.3f} / {two_phase[-2]:+.3f}",
            flush=True,
        )
        for step in sorted({min(20, args.steps), args.steps}):
            p = profiles[step]
            row = ", ".join(f"{x}: {p[x]:+.3f} / {p[-1 - x]:+.3f}" for x in PROBES)
            print(
                f"seed {seed} step {step}: F / scale at distance (left / right) {row}; change "
                f"since t = 0 at 3 cells {p[3] - profiles[0][3]:+.2e} / "
                f"{p[-4] - profiles[0][-4]:+.2e}; run {seconds:.1f} s",
                flush=True,
            )


if __name__ == "__main__":
    main()
