#!/usr/bin/env python3
"""
Self-force of one particle on ONE level, fully kinetic, 1D periodic, 64 cells, PHARE units: the
floor the P4 gate T10 compares against at a refinement boundary, with the options of the refined
runs (--filter-passes 1 --cleaning 0.5) or without.

Frozen Gauss (lesson): the FK start has E^0 = 0 and the deposit conserves charge, so div E - c^2 rho
keeps its t = 0 value forever. A charge present at t = 0 never gets its Coulomb field (no light
crossing builds it, Gauss is a constraint): it carries a frozen opposite charge where it started.

Start. pyphare refuses a non-neutral start (fully_kinetic.py, "frozen background"). Here the
start is neutral: the test particle (charge -1, mass m_e of waves.py) and its partner (charge +1,
mass 1e6, at rest) share one cell, each alone in its population (density in that one cell only,
1 particle per cell, v_th = 1e-12: the C++ normal distribution needs v_th > 0). When the test
particle drifts away the deposit builds the exact discrete Coulomb field of the pair: no frozen
charge left, except the dipole of the two positions inside the cell. Without cleaning, a particle
at rest therefore has no field at all, and its static self-force is the one of the slow drift: in
1D the longitudinal E is the Gauss field of the charge at the segment end x^{n+1/2}, so the force
of step n is the static self-force at x^{n+1/2}, for any drift speed.

Force. From the particle dumps (float64): E_felt = (m/q) (v^{n+1} - v^n) / dt, the field the
pusher applied at x^{n+1/2} = x^n + v^n dt / 2 (Boris, half drifts). The field of everything else
is the mean of the two plateaus of the field dump around the particle, just outside its deposit
(widened by one node per filter pass) and gather stencils, E_ext = (E_left + E_right) / 2, and
E_self = E_felt - E_ext. Its scale is the jump of the field across the particle,
Q = E_right - E_left = c^2 q w dx: the Coulomb field of the partner, a charge sheet, is Q/2.
Cases:
- rest: test particle at rest. Without cleaning, 100 steps, the largest |dv|. With cleaning, both
  charges heavy, 6000 steps: E / Q felt, against the Marder start-up field below and against
  the converged field of the pair (Gauss, zero mean)
- x: drift v_x = 0.1 v_th,e = 0.1 (v_th,e = sqrt(Te/m_e) = 1 in waves.py, 0.01 c); once its
  stencil has left the partner's, 400 steps over 2 cells, E_self against the sub-cell position
- y: drift v_y = 0.1 at four sub-cell positions (cells chosen from the loader's positions):
  J_y is a direct deposit, B_z gathered like E_x: magnetic self-force F_x / q - E_x = v_y B_self
  (E_x gathered from the dumps), scale the jump of B_z, Q v_y / c^2. Steps within 8 steps of the
  light fronts the start of the current sends around the box (period L/c = 128 steps) are not
  used.
Prediction: the discrete Coulomb field of the deposit (B-spline of the order on the primal nodes,
Gauss on the Yee grid, plateaus -Q/2, +Q/2) gathered with the same B-spline on the dual nodes:
E_self / Q = g_p(delta), the same function for B_z / (Q v_y / c^2).
Filter: in 1D, E_x only changes by -c^2 dt J_x (curl B has no x component), from the start
E^{-1/2} built with the filtered J_back too, so E_x is the filtered Gauss field: the (1/4, 1/2,
1/4) filter is linear, shift invariant and commutes with the Yee difference, the same as
depositing the smoothed shape F*S and gathering with S. The same for B_z through J_y. One pass:
g_1 = delta (1/2 - delta) / 4 on [0, 1/2] (a quarter of the unfiltered one, peak 1/64), peaks
0.008301 and 0.001636 at orders 2 and 3 (0.243 and 0.166 of the unfiltered ones).
Cleaning: the start's Gauss residual F(0) = -c^2 rho(x^0) is the in-cell dipole of the pair, and
the deposit keeps F between two cleanings, so the field the cleanings build from it does not
depend on the motion of the particles: marder_field computes it exactly, and it is taken out of
the felt and dumped fields in the x case. It decays like the slowest mode, 1 - s sin^2(pi / 64)
per step: -1.3e-2 of Q off the converged field at step 100, 9e-4 at step 1000, 7e-7 at 4000.

Measured, Debug, C++ at dd9d3aa37 plus uncommitted part 2 edits (one level only), Q = 1e-4:
- rest without cleaning: largest |dv| per step 2.6e-18, 8.4e-19, 2.2e-18 (orders 1, 2, 3)
- x: |E_self / Q - g_p| <= 1.4e-8 (N = 0 or 1, s = 0), <= 2.0e-8 (s = 0.5, N = 0 or 1) at
  every step of the sweep; peaks:
      N, s     order 1    order 2    order 3
      0, 0     0.062500   0.034179   0.009833   (12.5%, 6.8%, 2.0% of the partner's Q/2)
      1, 0     0.015625   0.008301   0.001636   (3.1%, 1.7%, 0.33%)
      0, 0.5   0.062500   0.034179   0.009833
      1, 0.5   0.015625   0.008301   0.001636
  Bounded: the same extremes over both cells of the sweep, mean over a cell ~1e-4 of Q.
- rest with cleaning, both heavy, test at delta 0.249, partner at 0.579: E / Q matches the Marder
  start-up model to 3e-12 to 2e-11 (step 100) up to 7e-8 (step 6000: the test charge drifts,
  4e-8 cells by estimate), and the converged pair field to 2e-8 to 6e-8 at step 6000 (orders 1 to 3, N = 0, 1)
- y: B_self / jump within 1e-3 of g_p at N = 0 (0.0626 / 0.0347 / 0.0097 at delta 0.25, model
  0.0625 / 0.0342 / 0.0098), the same with s = 0.5 (within 3e-4 of s = 0); at N = 1 0.0157 /
  0.0076 / 0.0032 against 0.0156 / 0.0083 / 0.0016: the halves of the run differ by up to 4e-3
  (per-step scatter 1e-3 to 3e-2 of the jump, the light fronts' dispersive tails cross the
  stencil), so the y case resolves g_p to a few 1e-3 of the jump only: report only.
--check: orders 1 to 3 in one process, x drift only: |E_self / Q - g_p(delta)| <= 1e-7 at every
step of the sweep, and the peak |E_self / Q| equal to PEAKS within 1%. Registered (level 9,
label physics) twice: N = 0, s = 0 (11.6 to 13.8 s), and --filter-passes 1 --cleaning 0.5, the
reference the refinement boundary gate subtracts (16.3 s).
Reading: the self-force is the Yee staggering (charge on primal nodes, field on dual nodes,
gathered with the same B-spline), not the Esirkepov averaging: in 1D the field of the segment
deposit is the Gauss field of the shape at the segment end, and the direct J_y deposit gives the
same g_p for B_z. The filter divides it by 4 to 6; the cleaning does not change it.

T10 at the refinement boundary (--refinement-box [lo,hi], default 16,48 = half the domain, hi
excluded; one filter pass and s = 0.5 unless given): the pair in one cell, --outside 4 cells
before the box edge (case in) or --inside 12 cells inside it (case out), the test charge drifting
+x at 0.01 c through the edge (3521 steps, every step dumped, both levels), and the same run on one
level. In the box the test charge is its split children (2 per particle in 1D, split_1d.hpp): per
coarse step the composite force is the children's in the box at both ends, weighted by their
charge fraction, plus the coarse particle's for the rest; a child entering (a level ghost copy
pushed in) or erased (pushed out) is a hand-off impulse. Extra force = composite force minus the
isolated self-force of its carrier (coarse g_p, or the fine cloud: `cloud`), minus the same
difference in the one-level run: the partner's field, the Marder start-up field (marder_field) and
the drift of the mean field cancel; in case out the steps where the pair's stencils still overlap
are left out (split in the box run, coarse in the other). Reported per cell crossed (one carrier),
and per crossing window (part of the cloud in the box): the net extra momentum (Q steps) and the
child impulses. Gate (--t10-check, orders 1 to 3, not pinned): extra <= the one-level peak; per
case, the per-step rows are saved (crossing_<case>_<order>_<tag>.npy).
Predicted by the two-grid model (docs/superpowers/specs/p4/t10_interface_model.py, the coupling as
the code does it, x drift: B is zero, the B ghost does not enter), max |extra| / Q, in / out:
      s      order 1          order 2          order 3       where
      0.5    0.38 / 0.49      0.33 / 0.43      0.29 / 0.39   cells -2 .. +1 of the edge
      0      0.034 / 0.030    0.025 / 0.022    0.022 / 0.016 cells -1, 0
(24-31, 40-52, 180-240 one-level peaks at s = 0.5): the gate will be red. At s = 0.5 the fine
Marder at the two boundary nodes dominates (their F reads the coarse E ghost; with F = 0 there:
0.10 / 0.09 at order 1); the coarse Marder under the box (plan D7(iv)) changes nothing; the E
dual child at the particle adds <= 3e-3 (the 3/4, 1/4 child 4.4e-2 instead of 3.4e-2). Crossing
windows: net +13 to +27 Q steps (mean 0.1 to 0.2 Q); out, the two child impulses are about
-125 and +150 Q steps (the children repel each other in the box; the coarse particle has the mean
velocity) and nearly cancel. Measured (11-17), max |extra| / Q, in / out:
      0.5    0.385 / 0.488    0.325 / 0.415    0.284 / 0.388
      0      0.034 / 0.030    0.026 / 0.022    0.022 / 0.016
crossing windows, s = 0.5: in +19 to +26, out -26 to -28 Q steps (model out -21 to -22); s = 0:
in +13 / +16 / +18, out +7 / +6 / +7 (as the model). The run total from the last velocities, box
run minus one-level run, matches the sum of the rows within 1 Q step.
"""

import argparse
import time
from dataclasses import dataclass
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
JUMP = -(c_over_va**2) * DENSITY * dl  # Q = c^2 q w dx of the test charge (-1)
REST_STEPS = 6000  # Marder: the slowest mode decays by 1 - s sin^2(pi / 64) per step
diag_root = "phare_outputs/fully_kinetic_self_force"


@dataclass(frozen=True)
class Options:
    passes: int = 0  # current_filter_passes
    cleaning: float = 0.0  # divergence_cleaning s

    @property
    def tag(self):
        return f"N{self.passes}_s{self.cleaning:g}"


BASE = Options()
AMR = Options(passes=1, cleaning=0.5)  # the refined runs' operating point (P4 D7)


def bspline(u, order):
    """centred cardinal B-spline of degree `order`, the particle shape, u in cells"""
    h = (order + 1) / 2
    return sum(
        (-1) ** k * comb(order + 1, k) * np.maximum(u + h - k, 0) ** order
        for k in range(order + 2)
    ) / factorial(order)


def filtered(field, passes):
    """(1/4, 1/2, 1/4) passes on a field with constant ends"""
    for _ in range(passes):
        f = np.concatenate([field[:1], field, field[-1:]])
        field = 0.25 * f[:-2] + 0.5 * f[1:-1] + 0.25 * f[2:]
    return field


def predicted(delta, order, passes=0):
    """E_self / Q of a unit charge at delta in cell 0: Gauss on the Yee grid (charge on the
    primal nodes, field on the dual nodes, plateaus -1/2 and +1/2), filtered as the current is
    (the same as filtering the charge: the filter commutes with the Yee difference), gathered on
    the dual nodes"""
    nodes = np.arange(-order - 2 - passes, order + 3 + passes)
    field = -0.5 + np.cumsum(bspline(nodes - delta, order))  # dual node j + 1/2
    return np.sum(bspline(nodes + 0.5 - delta, order) * filtered(field, passes))


def pair_charge(x_test, x_partner, order, passes):
    """charge of the test charge (x in cells) and of its opposite partner on the periodic primal
    nodes, filtered as the deposit, in units of Q / (c^2 dx), Q = c^2 q w dx of the test charge"""
    nodes = np.arange(CELLS)

    def charge(x):
        u = (nodes - x + CELLS / 2) % CELLS - CELLS / 2
        return bspline(u, order)

    rho = charge(x_test) - charge(x_partner)
    for _ in range(passes):
        rho = 0.25 * np.roll(rho, 1) + 0.5 * rho + 0.25 * np.roll(rho, -1)
    return rho


def predicted_pair(x_test, x_partner, order, passes):
    """E / Q at the test charge of the pair at rest once Marder has converged: E satisfies Gauss
    on the periodic grid with a zero mean (Marder leaves the mean of E unchanged, zero here)"""
    field = np.cumsum(
        pair_charge(x_test, x_partner, order, passes)
    )  # dual node j + 1/2
    return gather(field - field.mean(), x_test, order)


def marder_field(rho0, cleaning, cleanings):
    """E / Q (dual nodes) built by the first k Marder cleanings, k = 0 ... cleanings, from the
    Gauss residual of the start F(0) = -c^2 rho(x^0) (E^0 = 0), rho0 from pair_charge. The
    deposit conserves charge, so F only changes at the cleanings and this field does not depend
    on how the particles move. With f = F dx / Q and d dt = s dx^2 / 2, a cleaning is
    e += (s/2) (f_{j+1} - f_j), then f += (s/2) (f_{j+1} - 2 f_j + f_{j-1}). E^{n+1/2} has had
    n + 1 cleanings."""
    f, e = -rho0, np.zeros(CELLS)
    fields = [e]
    for _ in range(cleanings):
        e = e + 0.5 * cleaning * (np.roll(f, -1) - f)
        f = f + 0.5 * cleaning * (np.roll(f, -1) - 2 * f + np.roll(f, 1))
        fields.append(e)
    return np.array(fields)


def start_field(diag_dir, order, options, cleanings):
    """the Marder field of the start of a run (zero without cleaning), from the t = 0 dumps"""
    if not options.cleaning:
        return np.zeros((cleanings + 1, CELLS))
    x_test, _ = particles(diag_dir)
    x_partner, _ = particles(diag_dir, "partner")
    rho0 = pair_charge(x_test[0, 0], x_partner[0, 0], order, options.passes)
    return marder_field(rho0, options.cleaning, cleanings)


def gather(field, x, order):
    """a dual field (node j at j + 1/2) at x (cells), periodic, with the particle shape"""
    nodes = np.arange(int(np.floor(x)) - order - 1, int(np.floor(x)) + order + 2)
    return np.sum(bspline(nodes + 0.5 - x, order) * field[nodes % CELLS])


def single_cell(cell, density):
    def fn(x):
        return np.where(np.floor(x / dl + 1e-9) == cell, density, 0.0)

    return fn


def population(charge, mass, density, vx=0.0, vy=0.0):
    pop = waves.maxwellian(charge, mass, 1.0, 1e-24 * mass, ppc=1, seed=1)
    pop.update(density=density, vbulkx=waves.uniform(vx), vbulky=waves.uniform(vy))
    return pop


def simulate(
    name, interp_order, steps, populations, dumps, options, fields=("E",), box=None
):
    """box: (lo, hi) coarse cells, hi excluded, one static refined box"""
    diag_dir = f"{diag_root}/{name}_{interp_order}_{options.tag}"
    refinement = {}
    if box is not None:
        diag_dir += f"_box{box[0]}-{box[1]}"
        refinement = {"refinement_boxes": {"L0": {"B0": [(box[0],), (box[1] - 1,)]}}}
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=steps,
        cells=[CELLS],
        dl=[dl],
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=options.passes,
        divergence_cleaning=options.cleaning,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=False,  # a cold two-particle box, not a plasma
        **refinement,
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
    for pop in populations:
        ph.ParticleDiagnostics(
            quantity="domain", population_name=pop, write_timestamps=timestamps
        )
    start = time.perf_counter()
    Simulator(sim).run().reset()
    return diag_dir, time.perf_counter() - start


def times_of(path):
    import h5py

    with h5py.File(path, "r") as f:
        return sorted(f["t"].keys(), key=float)


def particles(diag_dir, pop="test"):
    """x (cells) and v per dump of the particles of a population, arrays [dump, particle, ...]"""
    import h5py

    path = f"{diag_dir}/ions_pop_{pop}_domain.h5"
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


def plateaus(field, x, order, passes):
    """values of a dual field just outside the stencils of a particle at x (cells): left of
    its first deposit node, right of its last, each filter pass widening the deposit by one"""
    first = int(np.floor(x - (order + 1) / 2)) + 1 - passes
    last = int(np.ceil(x + (order + 1) / 2)) - 1 + passes
    return field[(first - 1) % CELLS], field[last % CELLS]


def self_field(felt, field, x, order, passes):
    left, right = plateaus(field, x, order, passes)
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
        options=BASE,
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


def rest(order, deltas, options):
    """without cleaning: the largest |dv| over 100 steps. With cleaning, both charges heavy
    (they stay where they are): E / Q felt by the test charge at a few steps up to REST_STEPS,
    against the pair model"""
    cell, _ = cell_for(0.25, deltas)
    if not options.cleaning:
        diag_dir, seconds = simulate(
            "rest", order, 100, pair(cell), range(101), options
        )
        _, v = particles(diag_dir)
        return np.abs(np.diff(v[:, 0, :], axis=0)).max(), seconds

    populations = {
        "test": population(-1, PARTNER_MASS, single_cell(cell, DENSITY)),
        "partner": population(1, PARTNER_MASS, single_cell(cell, DENSITY)),
    }
    steps = [100, 500, 1000, 2000, 4000, REST_STEPS - 1]
    dumps = sorted({0, *(k for n in steps for k in (n, n + 1))})
    diag_dir, seconds = simulate("rest", order, REST_STEPS, populations, dumps, options)
    x_test, v = particles(diag_dir)
    x_partner, _ = particles(diag_dir, "partner")
    converged = predicted_pair(x_test[0, 0], x_partner[0, 0], order, options.passes)
    e = start_field(diag_dir, order, options, REST_STEPS)
    rows = []
    for n in steps:
        i = dumps.index(n)
        felt = PARTNER_MASS / -1 * (v[i + 1, 0, 0] - v[i, 0, 0]) / time_step
        rows.append((n, felt / JUMP, gather(e[n + 1], x_test[i, 0], order)))
    return (x_test[0, 0] % 1, x_partner[0, 0] % 1, converged, rows), seconds


def drift_x(order, deltas, options):
    """E_self / Q against delta over the sweep, the jump Q, and the run time"""
    # the stencils of the test particle and of the partner overlap over order + 2 cells, one
    # more per filter pass
    separation = round((order + 3 + options.passes) / (VELOCITY * time_step / dl))
    steps = separation + SWEEP_STEPS
    cell, _ = cell_for(0.0, deltas, avoid=range(CELLS - 8, CELLS))
    dumps = [0, *range(separation, steps + 1)]
    diag_dir, seconds = simulate(
        "x", order, steps, pair(cell, vx=VELOCITY), dumps, options
    )
    e = JUMP * start_field(diag_dir, order, options, steps + 1)
    x, v = particles(diag_dir)
    E = dual_field(diag_dir, "Ex")
    x, vx, E = x[1:, 0], v[1:, 0, 0], E[1:]
    rows = []
    for n in range(len(x) - 1):
        cleanings = separation + n + 1  # of E^{n+1/2}
        felt = mass_e / -1 * (vx[n + 1] - vx[n]) / time_step
        half = x[n] + vx[n] * time_step / 2 / dl
        felt -= gather(e[cleanings], half, order)
        ratio, jump = self_field(
            felt, E[n + 1] - e[cleanings], half, order, options.passes
        )
        rows.append((half % 1, ratio, jump))
    return np.array(rows), seconds


def drift_y(order, cell, options):
    """B_self / (jump of B_z) over the usable steps, and the run time. The E_x the test charge
    feels is gathered from the dumps and taken out: zero without cleaning; with cleaning, the
    field Marder builds from the in-cell dipole of the start pushes the test charge along x, and
    its own E_x then dwarfs the magnetic self-force (v^2 / c^2 = 1e-4 smaller)"""
    steps = Y_STEPS
    diag_dir, seconds = simulate(
        f"y_{cell}",
        order,
        steps,
        pair(cell, vy=VELOCITY),
        range(steps + 1),
        options,
        fields=("B", "E"),
    )
    x, v = particles(diag_dir)
    B = dual_field(diag_dir, "Bz")
    E = dual_field(diag_dir, "Ex")
    rows = []
    for n in range(steps):
        if min(n % FRONT_PERIOD, FRONT_PERIOD - n % FRONT_PERIOD) < 8:
            continue
        vy = (v[n, 0, 1] + v[n + 1, 0, 1]) / 2
        felt = mass_e / -1 * (v[n + 1, 0, 0] - v[n, 0, 0]) / time_step
        half = x[n, 0] + v[n, 0, 0] * time_step / 2 / dl
        felt = (felt - gather(E[n + 1], half, order)) / vy
        ratio, jump = self_field(
            felt, (B[n] + B[n + 1]) / 2, x[n, 0], order, options.passes
        )
        rows.append((ratio, jump))
    return x[0, 0] % 1, np.array(rows), seconds


# g_p(1/4) per order
PEAKS = {
    BASE: {1: 1 / 16, 2: 0.0342, 3: 0.0098},
    AMR: {1: 1 / 64, 2: 0.0083, 3: 0.00164},
}


def check(options):
    for order in (1, 2, 3):
        rows, seconds = drift_x(order, delta_table(order), options)
        delta, ratio, _ = rows.T
        model = [predicted(d, order, options.passes) for d in delta]
        error = np.abs(ratio - model).max()
        peak, expected = np.abs(ratio).max(), PEAKS[options][order]
        print(
            f"{options.tag} order {order}: |E_self / Q - g_p| max {error:.1e}, "
            f"peak |E_self / Q| {peak:.6f}, run {seconds:.1f} s"
        )
        if error > 1e-7:
            raise RuntimeError(
                f"{options.tag} order {order}: E_self / Q off the model by {error:.1e}"
            )
        if abs(peak / expected - 1) > 0.01:
            raise RuntimeError(
                f"{options.tag} order {order}: peak {peak:.6f}, not {expected}"
            )


BOX = (16, 48)  # coarse cells, hi excluded: half the domain
SPLIT = {
    1: 0.551569,
    2: 0.663959,
    3: 0.752399,
}  # 1D 2-particle split, fine cells (split_1d.hpp)


def cloud(xs, fractions, order, passes):
    """sum_i f_i E(x_i) / Q of a cloud of charges f_j at x_j (cells of one grid, unit total): the
    isolated self-force of a split particle, `predicted` for a single charge"""
    xs, fractions = np.atleast_1d(xs), np.atleast_1d(fractions)
    nodes = np.arange(
        int(np.floor(xs.min())) - order - 3 - passes,
        int(np.floor(xs.max())) + order + 4 + passes,
    )
    field = sum(
        f * (np.cumsum(bspline(nodes - x, order)) - 0.5) for x, f in zip(xs, fractions)
    )
    field = filtered(field, passes)
    return sum(
        f * np.sum(bspline(nodes + 0.5 - x, order) * field)
        for x, f in zip(xs, fractions)
    )


def level_particles(diag_dir, pop, level):
    """per dump, (x in the level's cells, v_x, weight) arrays of a population on a level"""
    import h5py

    path = f"{diag_dir}/ions_pop_{pop}_domain.h5"
    out = []
    with h5py.File(path, "r") as f:
        for t in times_of(path):
            x, v, w = [], [], []
            for patch in f["t"][t].get(f"pl{level}", {}).values():
                if "iCell" not in patch or patch["iCell"].shape[0] == 0:
                    continue
                x.append(patch["iCell"][:, 0] + patch["delta"][:, 0].astype(np.float64))
                v.append(patch["v"][:, 0].astype(np.float64))
                w.append(patch["weight"][:].ravel().astype(np.float64))
            out.append(
                tuple(np.concatenate(a) if a else np.zeros(0) for a in (x, v, w))
            )
    return out


def composite(diag_dir, order, passes, refined):
    """per coarse step: position of the composite particle (coarse cells), its force minus the
    isolated self-force of what carries it (E / Q), the hand-off impulse (Q steps) and the
    carrier: "coarse" (no child in the box at either end), "fine" (the whole cloud in the box at
    both ends), "mixed". Children in the box at both ends of the step are matched by position;
    the force is theirs, weighted by their charge fraction (child weight / (2 coarse weight)),
    plus the coarse particle's for the rest; a child entering (a band copy pushed into the box)
    or erased (pushed out) changes the composite momentum by its fraction times its velocity
    minus the coarse one: the hand-off impulse"""
    coarse = level_particles(diag_dir, "test", 0)
    fine = level_particles(diag_dir, "test", 1) if refined else None
    u = (
        VELOCITY * time_step / dl
    )  # cells per step on either level (fine: dl/2 per dt/2)
    to_force = mass_e / -1 / time_step / JUMP
    rows = []
    for n in range(len(coarse) - 1):
        (xc, vc, wc), (_, vc1, _) = coarse[n], coarse[n + 1]
        xc, vc, vc1, wc = xc[0], vc[0], vc1[0], wc[0]
        kids = list(zip(*fine[n])) if refined else []
        kids1 = list(zip(*fine[n + 1])) if refined else []
        kids = [(x, v, w / (2 * wc)) for x, v, w in kids]
        kids1 = [(x, v, w / (2 * wc)) for x, v, w in kids1]
        matched = []
        for x, v, f in kids:
            gaps = [abs(x1 - x - 2 * u) for x1, _, _ in kids1]
            if gaps and min(gaps) < 0.1:
                matched.append((x, v, f, kids1[int(np.argmin(gaps))][1]))
        share = sum(m[2] for m in matched)
        force = to_force * (
            sum(f * (v1 - v) for _, v, f, v1 in matched) + (1 - share) * (vc1 - vc)
        )
        momentum = [
            sum(f * v for _, v, f in k) + (1 - sum(f for *_, f in k)) * c
            for k, c in ((kids, vc), (kids1, vc1))
        ]
        impulse = to_force * (momentum[1] - momentum[0]) - force
        f0, f1 = sum(f for *_, f in kids), sum(f for *_, f in kids1)
        X = sum(f * x / 2 for x, _, f in kids) + (1 - f0) * xc
        half = xc + vc * time_step / dl / 2
        if f0 == 0 and f1 == 0:
            carrier, ref = "coarse", predicted(half % 1, order, passes)
        elif min(f0, f1) > 1 - 1e-6 and len(matched) == len(kids):
            carrier = "fine"
            xs = [x for x, *_ in matched]
            fs = [f for _, _, f, _ in matched]
            ref = np.mean(
                [cloud([x + k * u for x in xs], fs, order, passes) for k in (0.5, 1.5)]
            )
        else:  # the box share of a cloud split from the coarse particle, the rest coarse
            carrier = "mixed"
            d = SPLIT[order]
            ref = share * cloud([2 * half - d, 2 * half + d], [0.5, 0.5], order, passes)
            ref += (1 - share) * predicted(half % 1, order, passes)
        rows.append((X, force - ref, impulse, carrier))
    return rows


def crossing(order, options, box, case, outside, inside):
    """T10 at the refinement boundary: the pair in one cell, `outside` cells before the box edge
    (in) or `inside` cells inside it (out), the test charge drifting +x through the edge; the box
    run against the one-level run. Rows: (position from the edge, extra E / Q, hand-off impulse,
    carrier, partner overlap)"""
    lo, hi = box
    cell = lo - outside if case == "in" else hi - inside
    steps = int(np.ceil((outside + inside) / (VELOCITY * time_step / dl) * 1.1))
    runs = {}
    for refined in (True, False):
        diag_dir, seconds = simulate(
            f"crossing_{case}",
            order,
            steps,
            pair(cell, vx=VELOCITY),
            range(steps + 1),
            options,
            box=box if refined else None,
        )
        runs[refined] = (composite(diag_dir, order, options.passes, refined), seconds)
    x_partner = particles(diag_dir, "partner")[0][0, 0]
    edge = lo if case == "in" else hi
    # out: the pair starts in the box, split there and coarse in the one-level run: their near
    # fields differ while the stencils overlap; in: the partner is coarse in both runs
    overlap = order + 3 + options.passes
    rows = [
        (X - edge, d - d1, imp, carrier, case == "out" and abs(X - x_partner) < overlap)
        for (X, d, imp, carrier), (_, d1, _, _) in zip(runs[True][0], runs[False][0])
    ]
    return rows, runs[True][1] + runs[False][1]


def crossing_report(order, rows):
    """max |extra| per cell crossed (one carrier, no partner overlap), and per crossing window
    (consecutive mixed steps): length, net extra momentum (Q steps), child hand-off impulses"""
    cells = {}
    for rel, d, _, carrier, near in rows:
        if carrier != "mixed" and not near:
            cells.setdefault(int(np.floor(rel)), []).append(abs(d))
    per_cell = {c: max(v) for c, v in sorted(cells.items())}
    windows, i = [], 0
    while i < len(rows):
        j = i
        while j < len(rows) and rows[j][3] == "mixed":
            j += 1
        if j > i:
            kicks = [r[2] for r in rows[i : j + 1] if abs(r[2]) > 1e-6]
            windows.append((j - i, sum(r[1] + r[2] for r in rows[i:j]), kicks))
        i = j + 1
    return per_cell, windows


def interface(orders, options, box, cases, outside, inside, gate):
    failures = []
    for order in orders:
        # the one-level self-force peak at these options (PEAKS at the registered ones)
        peak = max(
            abs(predicted(d, order, options.passes)) for d in np.linspace(0, 1, 401)
        )
        for case in cases:
            rows, seconds = crossing(order, options, box, case, outside, inside)
            np.save(
                f"{diag_root}/crossing_{case}_{order}_{options.tag}.npy",
                np.array([r[:3] for r in rows]),
            )
            per_cell, windows = crossing_report(order, rows)
            worst = max(per_cell.values())
            print(
                f"T10 {options.tag} order {order} {case}: max |extra| / Q {worst:.2e} = "
                f"{worst / peak:.2f} one-level peak ({peak:.6f}), runs {seconds:.1f} s"
            )
            print(
                "   per cell from the edge: "
                + " ".join(f"{c:+d}:{v:.1e}" for c, v in per_cell.items())
            )
            for n, net, kicks in windows:
                print(
                    f"   crossing window {n} steps: net extra momentum {net:+.2f} Q steps "
                    f"(mean {net / n:+.1e}), child hand-off impulses "
                    + " ".join(f"{k:+.2f}" for k in kicks)
                )
            if worst > peak:
                failures.append(f"order {order} {case}: {worst / peak:.2f} x the peak")
    if failures:
        message = "T10 interface gate (extra <= one-level peak): " + "; ".join(failures)
        if gate:
            raise RuntimeError(message)
        print("report only, " + message)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--order", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument(
        "--filter-passes",
        type=int,
        default=None,
        help="0, or 1 with --refinement-box",
    )
    parser.add_argument(
        "--cleaning", type=float, default=None, help="0, or 0.5 with --refinement-box"
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="x drift, orders 1 to 3, against the model: registered with the default options",
    )
    parser.add_argument(
        "--refinement-box",
        nargs="?",
        const=f"{BOX[0]},{BOX[1]}",
        help=f"T10 at the box edge: lo,hi coarse cells, hi excluded (default {BOX[0]},{BOX[1]})",
    )
    parser.add_argument("--case", choices=["in", "out", "both"], default="both")
    parser.add_argument(
        "--outside", type=int, default=4, help="in: start cells before the edge"
    )
    parser.add_argument("--inside", type=int, default=12, help="out: start cells inside")
    parser.add_argument(
        "--t10-check",
        action="store_true",
        help="orders 1 to 3, assert extra <= the one-level peak (not pinned yet)",
    )
    args = parser.parse_args()
    refined = args.refinement_box is not None or args.t10_check
    options = Options(
        args.filter_passes if args.filter_passes is not None else int(refined),
        args.cleaning if args.cleaning is not None else 0.5 * refined,
    )
    if refined:
        box = tuple(
            int(c) for c in (args.refinement_box or f"{BOX[0]},{BOX[1]}").split(",")
        )
        cases = ["in", "out"] if args.case == "both" else [args.case]
        orders = [1, 2, 3] if args.t10_check else [args.order]
        return interface(
            orders, options, box, cases, args.outside, args.inside, args.t10_check
        )
    if args.check:
        return check(options)
    order = args.order
    passes = options.passes
    deltas = delta_table(order)

    result, seconds = rest(order, deltas, options)
    if options.cleaning:
        x_test, x_partner, model, rows = result
        print(
            f"{options.tag} order {order} rest: test at {x_test:.3f}, partner at "
            f"{x_partner:.3f} (in their cell), converged pair model E / Q {model:+.8f}, "
            f"run {seconds:.1f} s"
        )
        for n, value, transient in rows:
            print(
                f"   step {n}: E / Q {value:+.8f}, - converged {value - model:+.1e}, "
                f"- Marder transient model {value - transient:+.1e}"
            )
    else:
        print(
            f"{options.tag} order {order} rest: largest |dv| per step {result:.1e}, "
            f"run {seconds:.1f} s"
        )

    rows, seconds = drift_x(order, deltas, options)
    delta, ratio, jump = rows.T
    model = np.array([predicted(d, order, passes) for d in delta])
    print(
        f"{options.tag} order {order} x drift: Q = {np.mean(jump):.4e} +- {np.std(jump):.1e}, "
        f"E_self / Q max {np.abs(ratio).max():.5f} (model {np.abs(model).max():.5f}), "
        f"|measured - model| max {np.abs(ratio - model).max():.1e}, "
        f"cell 1 / cell 2: max {np.abs(ratio[:200]).max():.5f} / {np.abs(ratio[200:]).max():.5f}, "
        f"mean {ratio[:200].mean():+.1e} / {ratio[200:].mean():+.1e}, "
        f"run {seconds:.1f} s"
    )
    for target in TARGETS:
        n = np.argmin(np.abs((delta - target + 0.5) % 1 - 0.5))
        print(
            f"   delta {delta[n]:.3f}: E_self / Q {ratio[n]:+.6f} model {model[n]:+.6f}, "
            f"- model {ratio[n] - model[n]:+.1e}"
        )

    used = set()
    for target in TARGETS:
        cell, _ = cell_for(target, deltas, avoid=used)
        used.add(cell)
        delta, rows, seconds = drift_y(order, cell, options)
        ratio, jump = rows.T
        half = len(ratio) // 2
        print(
            f"{options.tag} order {order} y drift, delta {delta:.3f}: B_self / jump "
            f"{ratio.mean():+.5f} +- {ratio.std():.1e} "
            f"(model {predicted(delta, order, passes):+.5f}), "
            f"halves {ratio[:half].mean():+.5f} / {ratio[half:].mean():+.5f}, "
            f"jump {np.mean(jump):.3e}, run {seconds:.1f} s"
        )


if __name__ == "__main__":
    main()
