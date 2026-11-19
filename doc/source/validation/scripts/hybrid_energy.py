#!/usr/bin/env python3
"""Energy budget of a uniform, periodic, undamped hybrid plasma (page numerics/hybrid_ohm).

    python3 hybrid_energy.py BDIR INTERP DT T OUT

2D, 24 x 24 cells, dx = 0.5, 50 particles per cell, vth = 0.3, n = 1, |B| = 1 along BDIR
(z: normal to the plane, xy: diagonal in the plane), Te = 0, no resistivity, no
hyper-resistivity. Nothing but the ions and B can hold energy, so their sum must stay
constant up to the time-integration error.

OUT gets one row per measurement: t, then the changes since t = 0 in % of the initial
energy E0 = K(0) + W_B(0): ion kinetic energy dK, magnetic energy dW_B, and their sum.
K = 1/2 sum_p w_p v_p^2, W_B = 1/2 sum B^2 per cell (level 0, serial run).
"""

import sys

import numpy as np

B_DIRS = {"z": (0.0, 0.0, 1.0), "xy": (1 / np.sqrt(2), 1 / np.sqrt(2), 0.0)}
CELLS, DL, PPC, VTH, SEED, MEASURES = 24, 0.5, 50, 0.3, 1, 20

_PRIMAL = {"Bx": (1, 0), "By": (0, 1), "Bz": (0, 0)}  # 2D: primal (1) or dual (0) along x, y


def _owned(patch_data, quantity):
    """a patch's values on its own cells plus its lower primal node: shared nodes counted once"""
    g = patch_data.nGhosts
    n = [int(u - l + 1) for l, u in zip(patch_data.lower, patch_data.upper)]
    shape = [n[d] + 2 * g + _PRIMAL[quantity][d] for d in range(2)]
    values = np.asarray(patch_data.data).reshape(shape)
    return values[tuple(slice(g, g + n[d]) for d in range(2))]


def energies(simulator):
    level = simulator.data_wrangler().getPatchLevel(0)
    W_B = 0.0
    for name, getter in (("Bx", level.getBx), ("By", level.getBy), ("Bz", level.getBz)):
        for patch in getter():
            W_B += 0.5 * (_owned(patch, name) ** 2).sum()
    K = 0.0
    for patch in level.getParticles()["protons"]["domain"]:
        w = np.asarray(patch.data.weight)
        v = np.asarray(patch.data.v).reshape(-1, 3)
        K += 0.5 * (w * (v**2).sum(axis=1)).sum()
    return K, W_B


def main(bdir, interp, dt, T, out):
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator

    steps = round(T / dt)
    ph.NO_GUI()
    sim = ph.Simulation(
        time_step=dt, time_step_nbr=steps, cells=[CELLS] * 2, dl=[DL] * 2, interp_order=interp,
        largest_patch_size=12, smallest_patch_size=12, resistivity=0.0, hyper_resistivity=0.0,
        strict=True,
    )

    def constant(v):
        return lambda x, y: v + 0.0 * x

    b = B_DIRS[bdir]
    ph.MaxwellianFluidModel(
        bx=constant(b[0]), by=constant(b[1]), bz=constant(b[2]),
        protons={"charge": 1, "mass": 1.0, "density": constant(1.0),
                 **{f"vbulk{c}": constant(0.0) for c in "xyz"},
                 **{f"vth{c}": constant(VTH) for c in "xyz"},
                 "nbr_part_per_cell": PPC, "init": {"seed": SEED}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    simulator = Simulator(sim, log_to_file=False).initialize()
    every = max(1, steps // MEASURES)
    rows = [(0.0, *energies(simulator))]
    for s in range(1, steps + 1):
        simulator.advance()
        if s % every == 0:
            rows.append((s * dt, *energies(simulator)))
    simulator.reset()

    K0, W0 = rows[0][1], rows[0][2]
    E0 = K0 + W0
    with open(out, "w") as f:
        f.write(f"# B along {bdir}, interp order {interp}, dt {dt}; changes in % of E0\n")
        f.write("# t dK dW_B dK+dW_B\n")
        for t, K, W in rows:
            dK, dW = (K - K0) / E0 * 100, (W - W0) / E0 * 100
            f.write(f"{t:.4f} {dK:.6e} {dW:.6e} {dK + dW:.6e}\n")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
