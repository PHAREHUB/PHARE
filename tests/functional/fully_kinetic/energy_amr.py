#!/usr/bin/env python3
"""
Energy of a 2D thermal plasma with one static refined box (P4 T11), against the same run on one
level, PHARE units. The setup of energy.py --short (16^2 cells, 24 particles per cell and species,
400 steps, order 1) at the refined runs' options (one filter pass, Marder s = 0.5, on both runs),
box over coarse cells 3..13 in x and y (121 of 256 cells).

Composite budget (plan section 12, pW's design): the fine level replaces the coarse one in the box.
  W = K + sum over levels of sum_nodes a (E^- . E^+ / c^2 + |B|^2) / 2 dV_l
- E^- . E^+ at each level's own half steps around its integer time t_n: the dump at t is E at
  t - dt_l / 2, so E^- is the dump at t_n and E^+ the dump at t_n + dt_l: t_{n+1} on the root,
  t_n + dt/2 on the fine level (fine dumps on every substep, fine_dump_lvl_max = 1); B at t_n
- a node's share a: the part of its control volume (one cell of its level, centred on it) inside
  the box for the fine level, outside it for the root: a node on the box face counts half on
  each side, a corner node a quarter in the box; fine ghost nodes get 0
- K = sum m w |v|^2 / 2 dV over the fine domain particles and the root particles whose cell is
  outside the box (the root particles under the box and the level ghost particles are left out)
No exact conservation is expected (band swap, erased leaving particles, reflection at the
interface): the crossings make W random-walk, ~ sqrt(number of crossings).

Runs per seed: one level (fields and particles), and the box run twice, bitwise the same physics:
fields with fine dumps on every substep, particles at the samples only (fine dumps write every
diagnostic at every substep: ~0.5 GB of particles per run otherwise). The two box runs must give
the same root E at the samples: checked.

Printed per seed and over seeds:
- W / W(0) drift per unit time (fit after the first sample), box and one level
- B_z energy growth per unit time per unit area, in the strip within STRIP coarse cells of the box
  boundary and elsewhere, over W(0) per unit area; their ratio
- the random walk: rms of the detrended increments of W and K between samples, over
  sqrt(expected crossings per interval) (Maxwellian flux through the box boundary, both ways, in
  root particles), in units of the mean root particle energy
Gates, placeholders to pin after the first run (--check asserts them): composite drift per unit
time <= one level mean + 3 sigma (3 seeds); strip B_z growth per unit area <= 2 x elsewhere.

Time step: energy.py's, c dt / dx = 0.5 per direction, c dt sqrt(2) / dx = 0.7071 = sin(pi / 4),
exactly the AMR coupling limit at ratio 2 (docs/superpowers/specs/p4/t10_t11_courant.md: finer
waves faster than the coarse step's Nyquist frequency alias in the synchronization): a strip
growth there may be that limit, not the particles crossing; --dt-scale 0.9 (0.64) separates them.
Not run yet.
"""

import argparse
from dataclasses import dataclass

import numpy as np
import waves
from energy import Te, Ti, masses
from waves import c_over_va, dl, mass_e

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

CELLS = 16
PPC = 24
STEPS = 400
SAMPLE_EVERY = 40
BOX = (3, 14)  # coarse cells, hi excluded, both directions
RATIO = 2
FILTER_PASSES, CLEANING = 1, 0.5
STRIP = 2  # coarse cells on each side of the box boundary
DRIFT_SIGMAS = 3
STRIP_BOUND = 2.0
diag_root = "phare_outputs/fully_kinetic_energy_amr"


@dataclass(frozen=True)
class Setup:
    order: int = 1
    dt_scale: float = 1.0  # time step over energy.py's

    @property
    def time_step(self):
        return waves.time_step * self.dt_scale

    def times(self, steps):
        return [round(s * self.time_step, 10) for s in steps]


def samples():
    return list(range(0, STEPS + 1, SAMPLE_EVERY))


def run(setup, seed, kind):
    """kind: "one" (one level), "fields" (box, E and B, fine dumps on every substep),
    "particles" (box, particles and E at the samples)"""
    order, times = setup.order, setup.times
    diag_dir = f"{diag_root}/order{order}_dt{setup.time_step:g}/seed{seed}_{kind}"
    lo, hi = BOX
    options = {"dir": diag_dir, "mode": "overwrite"}
    refinement = {}
    if kind != "one":
        refinement = {"refinement_boxes": {"L0": {"B0": [(lo, lo), (hi - 1, hi - 1)]}}}
    if kind == "fields":
        options["fine_dump_lvl_max"] = 1
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=setup.time_step,
        time_step_nbr=STEPS + 1,
        cells=[CELLS, CELLS],
        dl=[dl, dl],
        interp_order=order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=FILTER_PASSES,
        divergence_cleaning=CLEANING,
        diag_options={"format": "phareh5", "options": options},
        strict=True,
        **refinement,
    )
    ph.MaxwellianFluidModel(
        bx=waves.uniform(0.0),
        by=waves.uniform(0.0),
        bz=waves.uniform(0.0),
        protons=waves.maxwellian(1, 1.0, 1.0, Ti, ppc=PPC, seed=1 + 10 * seed),
        electrons=waves.maxwellian(-1, mass_e, 1.0, Te, ppc=PPC, seed=2 + 10 * seed),
    )
    s = samples()
    if kind in ("one", "fields"):
        field_times = times(sorted(set(s + [n + 1 for n in s])))
        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=field_times)
    if kind in ("one", "particles"):
        if kind == "particles":
            ph.ElectromagDiagnostics(quantity="E", write_timestamps=times(s))
        for pop in masses:
            ph.ParticleDiagnostics(
                quantity="domain", write_timestamps=times(s), population_name=pop
            )
    Simulator(sim).run().reset()
    return diag_dir


def inside_fraction(x, h):
    """per direction: the part of [x - h/2, x + h/2] (physical) inside the box"""
    lo, hi = BOX[0] * dl, BOX[1] * dl
    return np.clip((np.minimum(x + h / 2, hi) - np.maximum(x - h / 2, lo)) / h, 0, 1)


def boundary_distance(x, y):
    """distance (coarse cells) of a point to the box boundary"""
    lo, hi = BOX[0] * dl, BOX[1] * dl
    dx = np.maximum(lo - x, x - hi)  # > 0 outside along x
    dy = np.maximum(lo - y, y - hi)
    outside = np.hypot(np.maximum(dx, 0), np.maximum(dy, 0))
    inside = np.minimum(np.abs(np.minimum(dx, 0)), np.abs(np.minimum(dy, 0)))
    return np.where((dx > 0) | (dy > 0), outside, inside) / dl


def level_field(diag_dir, quantity, time, level):
    """{node key: (value, share, x, y)} of one component on a level: each physical node once
    (periodic root), share = its weight in the composite budget"""
    h = dl / RATIO**level
    L = CELLS * dl
    hier = hierarchy_from(
        h5_filename=f"{diag_dir}/EM_{quantity[0]}.h5", times=f"{time:.10f}"
    )
    out = {}
    for patch in hier.level(level, time).patches:
        pdata = patch.patch_datas[quantity]
        x, y = pdata.meshgrid()
        fraction = inside_fraction(x, h) * inside_fraction(y, h)
        share = fraction if level else 1 - fraction
        keep = share > 0
        if level == 0:
            keep &= (x > -1e-6) & (x < L - 1e-6) & (y > -1e-6) & (y < L - 1e-6)
        keys = np.round(2 * np.stack([x[keep], y[keep]], axis=1) / h).astype(int)
        values = pdata.dataset[:][keep].astype(np.float64)
        for key, v, a, px, py in zip(keys, values, share[keep], x[keep], y[keep]):
            out[tuple(key)] = (v, a, px, py)
    return out


def field_energy(setup, diag_dir, step, refined):
    """electric and magnetic composite energies, and the B_z energy in the strip and elsewhere
    with the areas they cover"""
    t, t_next = setup.times([step, step + 1])
    levels = [(0, t, t_next)]
    if refined:
        levels.append((1, t, round(t + setup.time_step / RATIO, 10)))
    electric = magnetic = 0.0
    bz = {"strip": 0.0, "elsewhere": 0.0}
    area = {"strip": 0.0, "elsewhere": 0.0}
    for level, t0, t1 in levels:
        volume = (dl / RATIO**level) ** 2
        for c in "xyz":
            before = level_field(diag_dir, f"E{c}", t0, level)
            after = level_field(diag_dir, f"E{c}", t1, level)
            electric += (
                sum(v * after[k][0] * a for k, (v, a, _, _) in before.items()) * volume
            )
            b = level_field(diag_dir, f"B{c}", t0, level)
            magnetic += sum(v * v * a for v, a, _, _ in b.values()) / 2 * volume
            if c == "z":
                for v, a, x, y in b.values():
                    region = (
                        "strip" if boundary_distance(x, y) <= STRIP else "elsewhere"
                    )
                    bz[region] += v * v * a / 2 * volume
                    area[region] += a * volume
    return electric / c_over_va**2 / 2, magnetic, bz, area


def kinetic_energy(setup, diag_dir, step, refined):
    (time,) = setup.times([step])
    lo, hi = BOX
    energy = 0.0
    for pop, mass in masses.items():
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/ions_pop_{pop}_domain.h5", times=f"{time:.10f}"
        )
        for level in (0, 1) if refined else (0,):
            volume = (dl / RATIO**level) ** 2
            for patch in hier.level(level, time).patches:
                particles = patch.patch_datas[f"{pop}_domain"].dataset
                if particles.size() == 0:
                    continue
                keep = np.ones(particles.size(), bool)
                if refined and level == 0:
                    cells = particles.iCells
                    keep = ~np.all((cells >= lo) & (cells < hi), axis=1)
                v2 = np.sum(particles.v[keep].astype(np.float64) ** 2, axis=1)
                energy += (
                    0.5 * mass * np.sum(particles.weights.ravel()[keep] * v2) * volume
                )
    return energy


def same_physics(setup, fields_dir, particles_dir):
    """the two box runs give the same root E at the samples (float32 dumps)"""
    for step in samples():
        (t,) = setup.times([step])
        for c in "xyz":
            a = level_field(fields_dir, f"E{c}", t, 0)
            b = level_field(particles_dir, f"E{c}", t, 0)
            if any(a[k][0] != b[k][0] for k in a):
                raise RuntimeError(f"the two box runs differ at step {step} (E{c})")


def crossings_per_interval(setup):
    """expected root particles crossing the box boundary per sample interval, both ways:
    sum over species of 2 n_p v_th / sqrt(2 pi) x perimeter x interval"""
    density = PPC / dl**2  # root particles per unit area and species
    perimeter = 4 * (BOX[1] - BOX[0]) * dl
    interval = SAMPLE_EVERY * setup.time_step
    return sum(
        2 * density * np.sqrt(T / m) / np.sqrt(2 * np.pi) * perimeter * interval
        for T, m in ((Ti, masses["protons"]), (Te, masses["electrons"]))
    )


def measure(setup, seed):
    s = samples()
    t = np.array(setup.times(s))

    def rate(y):
        return np.polyfit(t[1:], np.asarray(y)[1:], 1)[0]

    out = {}
    for refined in (False, True):
        if refined:
            fields_dir = run(setup, seed, "fields")
            particles_dir = run(setup, seed, "particles")
            same_physics(setup, fields_dir, particles_dir)
        else:
            fields_dir = particles_dir = run(setup, seed, "one")
        W, K, strip, elsewhere = [], [], [], []
        for step in s:
            k = kinetic_energy(setup, particles_dir, step, refined)
            e, m, bz, area = field_energy(setup, fields_dir, step, refined)
            W.append(k + e + m)
            K.append(k)
            strip.append(bz["strip"] / area["strip"])
            elsewhere.append(bz["elsewhere"] / area["elsewhere"])
        W, K = np.array(W), np.array(K)
        per_area = W[0] / (CELLS * dl) ** 2
        dW, dK = np.diff(W[1:]), np.diff(K[1:])
        out[refined] = {
            "drift": rate(W) / W[0],
            "strip": rate(strip) / per_area,
            "elsewhere": rate(elsewhere) / per_area,
            "rms dW": np.std(dW - dW.mean()),
            "rms dK": np.std(dK - dK.mean()),
            "particle energy": K[0] / (2 * PPC * CELLS**2),
        }
        r = out[refined]
        print(
            f"seed {seed} {'box' if refined else 'one level'}: W drift {r['drift']:+.3e} per "
            f"unit time; B_z growth per unit area / (W(0) per unit area) strip "
            f"{r['strip']:+.3e}, elsewhere {r['elsewhere']:+.3e}; rms increments W "
            f"{r['rms dW']:.2e}, K {r['rms dK']:.2e}",
            flush=True,
        )
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("seeds", nargs="*", type=int, default=[0, 1, 2])
    parser.add_argument("--order", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument(
        "--check", action="store_true", help="assert the placeholder gates"
    )
    parser.add_argument(
        "--dt-scale",
        type=float,
        default=1.0,
        help="time step over energy.py's (0.9: c dt sqrt(2) / dx = 0.64, below the AMR limit)",
    )
    args = parser.parse_args()
    setup = Setup(args.order, args.dt_scale)
    results = [measure(setup, seed) for seed in args.seeds]
    one = np.array([r[False]["drift"] for r in results])
    box = np.array([r[True]["drift"] for r in results])
    ratio = np.array([r[True]["strip"] / r[True]["elsewhere"] for r in results])
    n = crossings_per_interval(setup)
    walk_W = np.array([r[True]["rms dW"] for r in results]) / np.sqrt(n)
    walk_K = np.array([r[True]["rms dK"] for r in results]) / np.sqrt(n)
    unit = np.mean([r[True]["particle energy"] for r in results])
    bound = one.mean() + DRIFT_SIGMAS * one.std()
    print(
        f"drift per unit time: box {box.mean():+.3e} +- {box.std():.1e}, one level "
        f"{one.mean():+.3e} +- {one.std():.1e} (bound one level + {DRIFT_SIGMAS} sigma "
        f"{bound:+.3e}); strip / elsewhere B_z growth per unit area {ratio.mean():.2f} "
        f"(seeds {', '.join(f'{x:.2f}' for x in ratio)}, bound {STRIP_BOUND}); random walk "
        f"per crossing ({n:.0f} crossings per interval): W {walk_W.mean() / unit:.2e}, K "
        f"{walk_K.mean() / unit:.2e} of a root particle's mean energy; one level rms "
        f"increments of W {np.mean([r[False]['rms dW'] for r in results]):.2e}"
    )
    failures = []
    if np.any(box > bound):
        failures.append(f"drift above one level + {DRIFT_SIGMAS} sigma")
    if np.any(ratio > STRIP_BOUND):
        failures.append(f"strip B_z growth above {STRIP_BOUND} x elsewhere")
    if failures:
        message = "T11 (placeholders, not pinned): " + "; ".join(failures)
        if args.check:
            raise RuntimeError(message)
        print("report only, " + message)


if __name__ == "__main__":
    main()
