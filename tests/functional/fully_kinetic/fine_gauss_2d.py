#!/usr/bin/env python3
"""
T5 in 2D (P4): the Gauss residual F = div E - c^2 F^N rho on a fine level along its boundary, the
plasma of T2 2D (tests/simulator/test_fully_kinetic_amr.py 2 1: uniform thermal ions and
electrons T = 0.03, m_i / m_e = 25, small B and electron velocity waves along x and y), on 32^2
root cells with a box over root cells 8 to 23 in x and y (32 fine cells per side), for the
sawtooth fix (docs/superpowers/specs/p4/sawtooth_fix_spec.md, point 3). 1D: fine_gauss.py.

F on fine primal node (i, j) at a coarse time t, as in 1D: the dumped fine E (ghosts included, so
the boundary nodes use the ghosts the solver filled from the coarse E), rho at x - (dt_f / 2) v
from the domain particles of the fine level and the level ghost particles of the node's owner
patch (smallest GlobalId) only, filtered N times as the solver does. The solver deposits each
patch's level ghost band on that patch alone, and near a junction the bands of two patches
overlap and hold the same particles: summing every band counts them twice.
At coarse times and mid steps (fine dumps on every substep): per boundary node the jump (coarse
time mean - mid step mean) and the two-phase mean from step JUMP_FROM on.

Boundary nodes, by class: "corner" within GUARD fine cells of a box corner along its edge,
"junction" within GUARD fine cells of a fine patch junction (a fine patch boundary meeting the
edge), "clean" elsewhere. The fix keeps today's sawtooth near the corners and junctions only.
--largest-patch-size 16 (default) tiles the fine level 2 x 2, one junction per edge (the run
prints the patch boxes); 32 gives one fine patch, no junction.

Units: T2's scale c^2 max |rho_root| at t = 0, as fine_gauss.py, which comes from the particle
noise and so changes with ppc; and J0 = c^2 n (u_e - u_i) dt_f / (dx_f sqrt(2 pi)), the charge
density the fine domain particles leaving through the boundary in one fine step put on it (one
way Maxwellian flux, order 1 weight ~1), which does not: the 1D jump at ppc 1600 (o1 N0 s0.5,
seed 0) is 1.19 / 1.18 J0. --compare <npz> prints, per class, the change of |jump| against
another run of the same options (before / after the fix). Report only, no gate.
"""

import argparse
import itertools
import time

import numpy as np
from fine_gauss import (
    C_OVER_VA,
    DL,
    DT,
    JUMP_FROM,
    MASS_E,
    RATIO,
    T,
    bspline,
)

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

CELLS = 32
BOX = (8, 23)  # root cells, both included, both directions
GUARD = 3  # fine cells
CLASSES = ("clean", "corner", "junction")
diag_root = "phare_outputs/fully_kinetic_fine_gauss_2d"


def uniform(value):
    return lambda x, y: value + 0.0 * x


def wave(amplitude, direction):
    return lambda x, y: amplitude * np.sin(2 * np.pi * (x, y)[direction] / (CELLS * DL))


def simulate(args, seed):
    tag = (
        f"o{args.order}_N{args.filter_passes}_s{args.cleaning:g}_ppc{args.ppc}"
        f"_steps{args.steps}_lps{args.largest_patch_size}"
    )
    diag_dir = f"{diag_root}/{tag}_seed{seed}"
    lo, hi = BOX
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=DT,
        time_step_nbr=args.steps,
        cells=[CELLS, CELLS],
        dl=[DL, DL],
        interp_order=args.order,
        model_options=["FullyKineticModel"],
        c_over_va=C_OVER_VA,
        current_filter_passes=args.filter_passes,
        divergence_cleaning=args.cleaning,
        refinement_boxes={"L0": {"B0": [(lo, lo), (hi, hi)]}},
        smallest_patch_size=8,
        largest_patch_size=args.largest_patch_size,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite", "fine_dump_lvl_max": 1},
        },
        strict=True,
    )
    common = {"density": uniform(1.0), "nbr_part_per_cell": args.ppc}
    ph.MaxwellianFluidModel(
        bx=uniform(1.0),
        by=wave(0.1, 0),
        bz=wave(0.1, 1),
        protons={
            "charge": 1,
            "mass": 1.0,
            **{f"vbulk{k}": uniform(0.0) for k in "xyz"},
            **{f"vth{k}": uniform(np.sqrt(T)) for k in "xyz"},
            "init": {"seed": 2026 + 10 * seed},
            **common,
        },
        electrons={
            "charge": -1,
            "mass": MASS_E,
            "vbulkx": wave(0.05, 0),
            "vbulky": wave(0.05, 1),
            "vbulkz": uniform(0.0),
            **{f"vth{k}": uniform(np.sqrt(T / MASS_E)) for k in "xyz"},
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


def particles(diag_dir, pop, quantity, t, level, patch_name=None):
    """positions in cells of the level (global index + delta), (v_x, v_y), weights; of the
    patch patch_name only if given"""
    import h5py

    xs, vs, ws = [], [], []
    with h5py.File(f"{diag_dir}/ions_pop_{pop}_{quantity}.h5", "r") as f:
        key = next(k for k in f["t"] if abs(float(k) - t) < 1e-9)
        group = f["t"][key].get(f"pl{level}")
        for name, patch in group.items() if group is not None else []:
            if "iCell" not in patch or patch_name not in (None, name):
                continue
            xs.append(patch["iCell"][:, :2] + patch["delta"][:, :2])
            vs.append(patch["v"][:, :2])
            ws.append(patch["weight"][:].ravel())
    if not xs:
        return np.zeros((0, 2)), np.zeros((0, 2)), np.zeros(0)
    return np.concatenate(xs), np.concatenate(vs), np.concatenate(ws)


def deposit(rho, pos, weights, lower, order, periodic=False):
    """adds weights x the product b-spline to rho, whose node (0, 0) is the global node lower"""
    reach = (order + 1) // 2 + 1
    base = np.floor(pos).astype(int)
    shape = np.asarray(rho.shape)
    for offsets in itertools.product(range(1 - reach, reach + 1), repeat=2):
        node = base + np.asarray(offsets)
        w = weights * np.prod(bspline(pos - node, order), axis=1)
        index = node - np.asarray(lower)
        if periodic:
            index %= shape
        else:
            keep = np.all((index >= 0) & (index < shape), axis=1)
            index, w = index[keep], w[keep]
        np.add.at(rho, tuple(index.T), w)


def charge(
    diag_dir, t, level, lower, n, order, quantities, periodic=False, patch_name=None
):
    """rho at x - (dt_l / 2) v on the n x n global primal nodes from lower, unfiltered"""
    dl, dt = DL / RATIO**level, DT / RATIO**level
    rho = np.zeros((n, n))
    for pop, q in (("protons", 1.0), ("electrons", -1.0)):
        for quantity in quantities:
            x, v, w = particles(diag_dir, pop, quantity, t, level, patch_name)
            deposit(rho, x - 0.5 * dt * v / dl, q * w, lower, order, periodic)
    return rho


def scale(diag_dir, order):
    """c^2 max |rho_root| at t = 0 (T2's scale), the root periodic"""
    rho = charge(diag_dir, 0.0, 0, (0, 0), CELLS, order, ("domain",), periodic=True)
    return C_OVER_VA**2 * np.abs(rho).max()


def outflux_charge():
    """J0: c^2 x the charge density of one fine step of leaving particles on a boundary node"""
    u_e, u_i = np.sqrt(T / MASS_E), np.sqrt(T)
    dt, dl = DT / RATIO, DL / RATIO
    return C_OVER_VA**2 * (u_e - u_i) * dt / (dl * np.sqrt(2 * np.pi))


def fine_nodes():
    """the global fine primal nodes of the box, both boundary nodes included: lower, count"""
    lo, hi = RATIO * BOX[0], RATIO * (BOX[1] + 1)
    return lo, hi - lo + 1


def binomial(rho, passes):
    """the current filter (1/4, 1/2, 1/4) along x then y, interior nodes only"""
    for _ in range(passes):
        for axis in (0, 1):
            r = np.moveaxis(rho, axis, 0).copy()
            r[1:-1] = 0.25 * r[:-2] + 0.5 * r[1:-1] + 0.25 * r[2:]
            rho = np.moveaxis(r, 0, axis)
    return rho


def fine_e(diag_dir, t):
    """per fine patch (its name in the dumps), {(2 x / dx_f, 2 y / dx_f): value} of its Ex and Ey,
    ghosts included: a level border face is a ghost of every patch along it, and each patch
    advances its own (its band deposit), so the residual of a node reads its owner's. And the
    largest disagreement between two patches on a node of the box of one of them"""
    dl = DL / RATIO
    hier = hierarchy_from(h5_filename=f"{diag_dir}/EM_E.h5", times=f"{t:.10f}")
    E, own, mismatch = {}, {"x": {}, "y": {}}, 0.0
    for patch in hier.level(1, t).patches:
        lower, upper = patch.box.lower * dl, (patch.box.upper + 1) * dl
        E[patch.id] = {}
        for c in "xy":
            pdata = patch.patch_datas[f"E{c}"]
            x, y = pdata.meshgrid()
            keys = [(round(2 * px / dl), round(2 * py / dl)) for px, py in zip(x.ravel(), y.ravel())]
            values = pdata.dataset[:].ravel().astype(float)
            E[patch.id][c] = dict(zip(keys, values))
            mine = (
                (x > lower[0] - 1e-9)
                & (x < upper[0] + 1e-9)
                & (y > lower[1] - 1e-9)
                & (y < upper[1] + 1e-9)
            ).ravel()
            for key, v, m in zip(keys, values, mine):
                if m:
                    if key in own[c]:
                        mismatch = max(mismatch, abs(own[c][key] - v))
                    own[c].setdefault(key, v)
    return E, mismatch


def node_owners(diag_dir, t):
    """per fine primal node of the box (n x n), the name of its owner patch: the patch of
    smallest GlobalId (rank, local index) among those holding the node, as the solver"""
    import h5py

    lo, n = fine_nodes()
    with h5py.File(f"{diag_dir}/ions_pop_protons_levelGhost.h5", "r") as f:
        key = next(k for k in f["t"] if abs(float(k) - t) < 1e-9)
        patches = f["t"][key]["pl1"]
        boxes = {
            name: (p.attrs["lower"][:2], p.attrs["upper"][:2] + 1)
            for name, p in patches.items()
        }
    owners = np.empty((n, n), dtype=object)
    for name in sorted(boxes, key=lambda s: tuple(map(int, s[1:].split("#"))))[::-1]:
        (xl, yl), (xu, yu) = (np.asarray(c) - lo for c in boxes[name])
        owners[xl : xu + 1, yl : yu + 1] = name
    return owners


def fine_residual(diag_dir, t, order, passes):
    """F on the fine primal node (global lo + (0..n-1) per direction) of the box, from the rho
    of its owner patch, and the largest patch disagreement on E"""
    dl = DL / RATIO
    lo, n = fine_nodes()
    nodes = ((lo - passes,) * 2, n + 2 * passes, order)
    domain = charge(diag_dir, t, 1, *nodes, ("domain",))
    owners = node_owners(diag_dir, t)
    rho = np.zeros((n, n))
    for name in set(owners.ravel()):
        band = charge(diag_dir, t, 1, *nodes, ("levelGhost",), patch_name=name)
        own = owners == name
        rho[own] = binomial(domain + band, passes)[passes : n + passes, passes : n + passes][own]
    E, mismatch = fine_e(diag_dir, t)
    div = np.zeros((n, n))
    for a, b in itertools.product(range(n), repeat=2):
        i, j = 2 * (lo + a), 2 * (lo + b)
        Ex, Ey = E[owners[a, b]]["x"], E[owners[a, b]]["y"]
        div[a, b] = (Ex[(i + 1, j)] - Ex[(i - 1, j)] + Ey[(i, j + 1)] - Ey[(i, j - 1)]) / dl
    return div - C_OVER_VA**2 * rho, mismatch


def junctions(diag_dir, t):
    """the fine patch boundaries strictly inside the box, per direction (global fine indices),
    and the patch boxes"""
    lo, n = fine_nodes()
    hier = hierarchy_from(h5_filename=f"{diag_dir}/EM_E.h5", times=f"{t:.10f}")
    boxes = [patch.box for patch in hier.level(1, t).patches]
    cuts = [
        sorted(
            ({int(b.lower[d]) for b in boxes} | {int(b.upper[d]) + 1 for b in boxes})
            - {lo, lo + n - 1}
        )
        for d in (0, 1)
    ]
    return cuts, boxes


def classify(cuts):
    """per node of the n x n grid: -1 inside, else the index in CLASSES of a boundary node.
    Along an edge x = const the node moves in y and meets the junctions of y (cuts[1])"""
    lo, n = fine_nodes()
    out = -np.ones((n, n), dtype=int)
    for a, b in itertools.product(range(n), repeat=2):
        along = []  # (coordinate along the edge, junctions on it)
        if a in (0, n - 1):
            along.append((b, cuts[1]))
        if b in (0, n - 1):
            along.append((a, cuts[0]))
        if not along:
            continue
        corner = min(min(s, n - 1 - s) for s, _ in along)
        junction = min((abs(lo + s - c) for s, cs in along for c in cs), default=n)
        out[a, b] = 1 if corner <= GUARD else 2 if junction <= GUARD else 0
    return out


def by_class(values, classes):
    """{class: (mean, sd, max |value|, nodes)} over the boundary nodes"""
    out = {}
    for k, name in enumerate(CLASSES):
        v = values[classes == k]
        if v.size:
            sd = v.std(ddof=1) if v.size > 1 else 0.0
            out[name] = (v.mean(), sd, np.abs(v).max(), v.size)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--filter-passes", type=int, default=1)
    parser.add_argument("--cleaning", type=float, default=0.5)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--ppc", type=int, default=50)
    parser.add_argument("--largest-patch-size", type=int, default=16)
    parser.add_argument("--seeds", type=int, default=1)
    parser.add_argument("--first-seed", type=int, default=0)
    parser.add_argument(
        "--compare", help="npz of another run, same options (before / after)"
    )
    args = parser.parse_args()
    J0 = outflux_charge()

    def history(diag_dir, times, s0):
        out = [
            fine_residual(diag_dir, t, args.order, args.filter_passes) for t in times
        ]
        return np.array([F for F, _ in out]) / s0, max(m for _, m in out)

    for seed in range(args.first_seed, args.first_seed + args.seeds):
        diag_dir, times, seconds = simulate(args, seed)
        s0 = scale(diag_dir, args.order)
        cuts, boxes = junctions(diag_dir, times[-1])
        classes = classify(cuts)
        profiles, m0 = history(diag_dir, times, s0)
        mid_times = times[:-1] + DT / RATIO
        mid_profiles, m1 = history(diag_dir, mid_times, s0)
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
            classes=classes,
            scale=s0,
            J0=J0,
            args=str(vars(args)),
        )
        print(
            f"seed {seed}: run {seconds:.1f} s; fine patches "
            f"{[(tuple(b.lower), tuple(b.upper)) for b in boxes]}, junctions x {cuts[0]} "
            f"y {cuts[1]}; scale {s0:.3f}, J0 {J0:.3f} (J0 / scale {J0 / s0:.3f}); largest "
            f"patch disagreement on E {max(m0, m1):.1e}",
            flush=True,
        )
        to_J0 = s0 / J0
        for name, values in (("jump", jump), ("two-phase mean", two_phase)):
            for cls, (mean, sd, top, count) in by_class(values, classes).items():
                print(
                    f"seed {seed} steps {JUMP_FROM}-{args.steps} {name}, {cls} ({count} "
                    f"nodes): mean {mean:+.3f} sd {sd:.3f} max |.| {top:.3f} of the scale; "
                    f"mean {mean * to_J0:+.3f} max |.| {top * to_J0:.3f} J0",
                    flush=True,
                )
        ring = np.zeros_like(
            classes, dtype=bool
        )  # one node in, along the clean stretches
        n = len(classes)
        for a, b in zip(*np.nonzero(classes == 0)):
            ring[min(max(a, 1), n - 2), min(max(b, 1), n - 2)] = True
        print(
            f"seed {seed} jump one node inside the clean nodes: mean "
            f"{jump[ring].mean():+.3f} sd {jump[ring].std(ddof=1):.3f} of the scale",
            flush=True,
        )
        if args.compare:
            other = np.load(args.compare)
            before = np.abs(other["jump"]) * float(other["scale"]) / float(other["J0"])
            after = np.abs(jump) * to_J0
            for k, cls in enumerate(CLASSES):
                change = (after - before)[classes == k]
                if change.size:
                    print(
                        f"seed {seed} |jump| change against {args.compare}, {cls}: mean "
                        f"{change.mean():+.3f} max {change.max():+.3f} J0",
                        flush=True,
                    )


if __name__ == "__main__":
    main()
