#!/usr/bin/env python3
"""
Self-force of one particle on ONE level, fully kinetic, 1D periodic, 64 cells, PHARE units: the
floor the P4 gate T10 compares against at a refinement boundary.

Frozen Gauss (lesson): the FK start has E^0 = 0 and the deposit conserves charge, so div E - c^2 rho
keeps its t = 0 value forever. A charge present at t = 0 never gets its Coulomb field (no light
crossing builds it, Gauss is a constraint): it carries a frozen opposite charge where it started.

Start. pyphare refuses a non-neutral start (fully_kinetic.py, "frozen background"). Here the
start is neutral: the test particle (charge
-1, mass m_e of waves.py) and its partner (charge +1, mass 1e6, at rest) share one cell, each
alone in its population (density in that one cell only, 1 particle per cell, v_th = 1e-12: the
C++ normal distribution needs v_th > 0). When the test particle drifts away the deposit builds
the exact discrete Coulomb field of the pair: no frozen charge left, except the dipole of the two
positions inside the cell. A particle at rest therefore has no field at all: its force is exactly
zero (measured below), and its static self-force is the one of the slow drift: in 1D the
longitudinal E is the Gauss field of the charge at the segment end x^{n+1/2}, so the force of
step n is the static self-force at x^{n+1/2}, for any drift speed.

Force. From the particle dumps (float64): E_felt = (m/q) (v^{n+1} - v^n) / dt, the field the
pusher applied at x^{n+1/2} = x^n + v^n dt / 2 (Boris, half drifts). The field of everything else
is the mean of the two plateaus of the field dump around the particle, just outside its deposit
and gather stencils, E_ext = (E_left + E_right) / 2, and E_self = E_felt - E_ext. Its scale is
the jump of the field across the particle, Q = E_right - E_left = c^2 q w dx: the Coulomb field
of the partner, a charge sheet, is Q/2. Cases:
- rest: test particle at rest, 100 steps, the largest |dv|
- x: drift v_x = 0.1 v_th,e = 0.1 (v_th,e = sqrt(Te/m_e) = 1 in waves.py, 0.01 c); once its
  stencil has left the partner's, 400 steps over 2 cells, E_self against the sub-cell position
- y: drift v_y = 0.1 at four sub-cell positions (cells chosen from the loader's positions):
  J_y is a direct deposit, B_z gathered like E_x: magnetic self-force F_x / q = v_y B_self,
  scale the jump of B_z, Q v_y / c^2. Steps within 8 steps of the light fronts the start of the
  current sends around the box (period L/c = 128 steps) are not used.
Prediction: the discrete Coulomb field of the deposit (B-spline of the order on the primal nodes,
Gauss on the Yee grid, plateaus -Q/2, +Q/2) gathered with the same B-spline on the dual nodes:
E_self / Q = g_p(delta), the same function for B_z / (Q v_y / c^2).

Measured, Debug, C++ at 1714f3d49, Q = 1e-4 (~9.5 s per order, all cases):
- rest: largest |dv| per step 2.6e-18, 8.4e-19, 2.2e-18 (orders 1, 2, 3): no field, no force
- x: E_self / Q matches the model to 1.4e-8, 1.3e-8, 1.3e-8. g_p is zero at delta = 0 and 1/2,
  odd around them, extreme at delta = 1/4 and 3/4: +-0.0625 (order 1, g_1 = delta (1/2 - delta)
  on [0, 1/2]), +-0.0342, +-0.0098, i.e. 12.5%, 6.8%, 2.0% of the partner's Coulomb field Q/2.
  Bounded: the same extremes over both cells of the sweep, mean over a cell ~1e-4 of Q.
- y: B_self / jump agrees with the same g_p within 1% at delta 0.25 and 0.75 (0.0625 / 0.0346 /
  0.0097 against 0.0625 / 0.0342 / 0.0098), equal over both halves of the run; per-step
  scatter 1e-3 to 3e-2 of the jump (the light fronts' dispersive tails cross the stencil).
--check (registered, level 9, label physics): orders 1 to 3 in one process, x drift only:
|E_self / Q - g_p(delta)| <= 1e-7 at every step of the sweep, and the peak |E_self / Q| equal
to 1/16, 0.0342, 0.0098 within 1%.
Reading: the self-force is the Yee staggering (charge on primal nodes, field on dual nodes,
gathered with the same B-spline), not the Esirkepov averaging: in 1D the field of the segment
deposit is the Gauss field of the shape at the segment end, and the direct J_y deposit gives the
same g_p for B_z.
"""

import argparse
import time
from math import comb, factorial

import numpy as np
import waves
from waves import c_over_va, dl, mass_e, time_step

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

CELLS = 64
DENSITY = 2.5e-5  # Q = c^2 n dx = 1e-4: the test particle slows by < 4% over a run
VELOCITY = 0.1  # 0.1 v_th,e of waves.py
PARTNER_MASS = 1e6
TARGETS = (0.0, 0.25, 0.5, 0.75)
SWEEP_STEPS = 400  # 2 cells at v = 0.1
Y_STEPS = 200
FRONT_PERIOD = round(CELLS * dl / c_over_va / time_step)  # 128 steps
diag_root = "phare_outputs/fully_kinetic_self_force"


def bspline(u, order):
    """centred cardinal B-spline of degree `order`, the particle shape, u in cells"""
    h = (order + 1) / 2
    return sum(
        (-1) ** k * comb(order + 1, k) * np.maximum(u + h - k, 0) ** order
        for k in range(order + 2)
    ) / factorial(order)


def predicted(delta, order):
    """E_self / Q of a unit charge at delta in cell 0: Gauss on the Yee grid (charge on the
    primal nodes, field on the dual nodes, plateaus -1/2 and +1/2), gathered on the dual nodes"""
    nodes = np.arange(-order - 2, order + 3)
    field = -0.5 + np.cumsum(bspline(nodes - delta, order))  # dual node j + 1/2
    return np.sum(bspline(nodes + 0.5 - delta, order) * field)


def single_cell(cell, density):
    def fn(x):
        return np.where(np.floor(x / dl + 1e-9) == cell, density, 0.0)

    return fn


def population(charge, mass, density, vx=0.0, vy=0.0):
    pop = waves.maxwellian(charge, mass, 1.0, 1e-24 * mass, ppc=1, seed=1)
    pop.update(density=density, vbulkx=waves.uniform(vx), vbulky=waves.uniform(vy))
    return pop


def simulate(name, interp_order, steps, populations, dumps, fields=("E",)):
    diag_dir = f"{diag_root}/{name}_{interp_order}"
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=steps,
        cells=[CELLS],
        dl=[dl],
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=False,  # a cold two-particle box, not a plasma
    )
    ph.MaxwellianFluidModel(
        bx=waves.uniform(0.0),
        by=waves.uniform(0.0),
        bz=waves.uniform(0.0),
        **populations,
    )
    timestamps = time_step * np.asarray(dumps)
    for quantity in fields:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    ph.ParticleDiagnostics(
        quantity="domain", population_name="test", write_timestamps=timestamps
    )
    start = time.perf_counter()
    Simulator(sim).run().reset()
    return diag_dir, time.perf_counter() - start


def times_of(path):
    import h5py

    with h5py.File(path, "r") as f:
        return sorted(f["t"].keys(), key=float)


def particles(diag_dir):
    """x (cells) and v per dump of the test particles, arrays [dump, particle, ...]"""
    import h5py

    path = f"{diag_dir}/ions_pop_test_domain.h5"
    xs, vs = [], []
    with h5py.File(path, "r") as f:
        for t in times_of(path):
            (patch,) = f["t"][t]["pl0"].values()
            xs.append(patch["iCell"][:, 0] + patch["delta"][:, 0])
            vs.append(patch["v"][:])
    return np.array(xs), np.array(vs)


def dual_field(diag_dir, quantity):
    """a dual component (E_x, B_z) on the physical dual nodes j + 1/2, per dump"""
    import h5py

    path = f"{diag_dir}/EM_{quantity[0]}.h5"
    values = []
    with h5py.File(path, "r") as f:
        for t in times_of(path):
            (patch,) = f["t"][t]["pl0"].values()
            data = patch[f"EM_{quantity[0]}_{quantity[1]}"][:].astype(np.float64)
            ghosts = (data.size - CELLS) // 2
            values.append(data[ghosts : ghosts + CELLS])
    return np.array(values)


def plateaus(field, x, order):
    """values of a dual field just outside the stencils of a particle at x (cells): left of
    its first deposit node, right of its last"""
    first = int(np.floor(x - (order + 1) / 2)) + 1
    last = int(np.ceil(x + (order + 1) / 2)) - 1
    return field[(first - 1) % CELLS], field[last % CELLS]


def self_field(felt, field, x, order):
    left, right = plateaus(field, x, order)
    jump = right - left
    return (felt - (left + right) / 2) / jump, jump


def delta_table(interp_order):
    """sub-cell position of the test population's particle in every cell (the loader's hash
    of seed, population and cell): one run without steps"""
    diag_dir, _ = simulate(
        "deltas",
        interp_order,
        1,
        {
            "test": population(-1, mass_e, waves.uniform(DENSITY)),
            "partner": population(1, mass_e, waves.uniform(DENSITY)),
        },
        dumps=[0],
    )
    x, _ = particles(diag_dir)
    x = np.sort(x[0])
    return x - np.floor(x)


def cell_for(target, deltas, avoid=()):
    distance = np.abs((deltas - target + 0.5) % 1 - 0.5)
    distance[list(avoid)] = np.inf
    cell = int(np.argmin(distance))
    return cell, deltas[cell]


def pair(cell, vx=0.0, vy=0.0):
    return {
        "test": population(-1, mass_e, single_cell(cell, DENSITY), vx=vx, vy=vy),
        "partner": population(1, PARTNER_MASS, single_cell(cell, DENSITY)),
    }


def rest(order, deltas):
    cell, _ = cell_for(0.25, deltas)
    diag_dir, seconds = simulate("rest", order, 100, pair(cell), range(101))
    _, v = particles(diag_dir)
    return np.abs(np.diff(v[:, 0, :], axis=0)).max(), seconds


def drift_x(order, deltas):
    """E_self / Q against delta over the sweep, the jump Q, and the run time"""
    # the stencils of the test particle and of the partner overlap over order + 2 cells
    separation = round((order + 3) / (VELOCITY * time_step / dl))
    steps = separation + SWEEP_STEPS
    cell, _ = cell_for(0.0, deltas, avoid=range(CELLS - 8, CELLS))
    dumps = range(separation, steps + 1)
    diag_dir, seconds = simulate("x", order, steps, pair(cell, vx=VELOCITY), dumps)
    x, v = particles(diag_dir)
    E = dual_field(diag_dir, "Ex")
    x, vx = x[:, 0], v[:, 0, 0]
    rows = []
    for n in range(len(x) - 1):
        felt = mass_e / -1 * (vx[n + 1] - vx[n]) / time_step
        half = x[n] + vx[n] * time_step / 2 / dl
        ratio, jump = self_field(felt, E[n + 1], half, order)
        rows.append((half % 1, ratio, jump))
    return np.array(rows), seconds


def drift_y(order, cell):
    """B_self / (jump of B_z) over the usable steps, and the run time"""
    steps = Y_STEPS
    diag_dir, seconds = simulate(
        f"y_{cell}",
        order,
        steps,
        pair(cell, vy=VELOCITY),
        range(steps + 1),
        fields=("B",),
    )
    x, v = particles(diag_dir)
    B = dual_field(diag_dir, "Bz")
    rows = []
    for n in range(steps):
        if min(n % FRONT_PERIOD, FRONT_PERIOD - n % FRONT_PERIOD) < 8:
            continue
        vy = (v[n, 0, 1] + v[n + 1, 0, 1]) / 2
        felt = mass_e / -1 * (v[n + 1, 0, 0] - v[n, 0, 0]) / time_step / vy
        ratio, jump = self_field(felt, (B[n] + B[n + 1]) / 2, x[n, 0], order)
        rows.append((ratio, jump))
    return x[0, 0] % 1, np.array(rows), seconds


PEAKS = {1: 1 / 16, 2: 0.0342, 3: 0.0098}  # g_p(1/4)


def check():
    for order in (1, 2, 3):
        rows, seconds = drift_x(order, delta_table(order))
        delta, ratio, _ = rows.T
        error = np.abs(ratio - [predicted(d, order) for d in delta]).max()
        peak = np.abs(ratio).max()
        print(
            f"order {order}: |E_self / Q - g_p| max {error:.1e}, "
            f"peak |E_self / Q| {peak:.5f}, run {seconds:.1f} s"
        )
        if error > 1e-7:
            raise RuntimeError(
                f"order {order}: E_self / Q off the model by {error:.1e}"
            )
        if abs(peak / PEAKS[order] - 1) > 0.01:
            raise RuntimeError(f"order {order}: peak {peak:.5f}, not {PEAKS[order]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--order", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument("--check", action="store_true", help="the registered test")
    args = parser.parse_args()
    if args.check:
        return check()
    order = args.order
    deltas = delta_table(order)

    dv, seconds = rest(order, deltas)
    print(f"order {order} rest: largest |dv| per step {dv:.1e}, run {seconds:.1f} s")

    rows, seconds = drift_x(order, deltas)
    delta, ratio, jump = rows.T
    model = np.array([predicted(d, order) for d in delta])
    print(
        f"order {order} x drift: Q = {np.mean(jump):.4e} +- {np.std(jump):.1e}, "
        f"E_self / Q max {np.abs(ratio).max():.4f} (model {np.abs(model).max():.4f}), "
        f"|measured - model| max {np.abs(ratio - model).max():.1e}, "
        f"cell 1 / cell 2: max {np.abs(ratio[:200]).max():.4f} / {np.abs(ratio[200:]).max():.4f}, "
        f"mean {ratio[:200].mean():+.1e} / {ratio[200:].mean():+.1e}, "
        f"run {seconds:.1f} s"
    )
    for target in TARGETS:
        n = np.argmin(np.abs((delta - target + 0.5) % 1 - 0.5))
        print(
            f"   delta {delta[n]:.3f}: E_self / Q {ratio[n]:+.5f} model {model[n]:+.5f}, "
            f"E_self {ratio[n] * jump[n]:+.3e}"
        )

    used = set()
    for target in TARGETS:
        cell, _ = cell_for(target, deltas, avoid=used)
        used.add(cell)
        delta, rows, seconds = drift_y(order, cell)
        ratio, jump = rows.T
        half = len(ratio) // 2
        print(
            f"order {order} y drift, delta {delta:.3f}: B_self / jump {ratio.mean():+.5f} "
            f"+- {ratio.std():.1e} (model {predicted(delta, order):+.5f}), "
            f"halves {ratio[:half].mean():+.5f} / {ratio[half:].mean():+.5f}, "
            f"jump {np.mean(jump):.3e}, run {seconds:.1f} s"
        )


if __name__ == "__main__":
    main()
