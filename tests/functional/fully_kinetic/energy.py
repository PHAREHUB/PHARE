#!/usr/bin/env python3
"""
Energy of a 2D thermal plasma (unmagnetized, kinetic ions and electrons), PHARE units:

    W^n = sum_p m w |v^n|^2 / 2 dV + sum_nodes (E^{n-1/2} . E^{n+1/2} / c^2 + |B^n|^2) / 2 dV

with w = n / ppc the particle weight. In vacuum the leapfrog Yee scheme conserves exactly this
field part (E at half steps, B at integer steps); with particles, momentum conserving PIC heats
(finite grid instability, noise), and W drifts. The dumped E is half a step behind its time
label: the dump at t_n is E^{n-1/2}, so E is dumped at t_n and t_{n+1}.

Setup: dx = lambda_De, w_pe dt = 0.1, 64 particles per cell and species, 1000 steps.
Measurement only: prints the relative drift of W per seed, then mean and mean + 3 sigma, over
the whole run and from the first sample on (E starts at 0, the field fluctuations take a few
plasma periods to fill: in 1D most of the drift is there, +3.7e-4 then flat within 1.5e-4).
Usage: energy.py [--1d] [seed ...] (--1d: 64 cells, a check of the script itself)
"""

import sys

import numpy as np
import waves
from waves import c_over_va, dl, mass_e, plasma_frequency, time_step

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

cells = (64,) if "--1d" in sys.argv else (32, 32)
ppc = 64
nsteps = 1000
sample_every = 100
# v_th,e = dx w_pe: dx = lambda_De
Te = mass_e * (dl * plasma_frequency(1, -1, mass_e)) ** 2
Ti = Te
masses = {"protons": 1.0, "electrons": mass_e}
cell_volume = dl ** len(cells)


def times(steps):
    return [round(step * time_step, 10) for step in steps]


def run(seed):
    diag_dir = f"phare_outputs/fully_kinetic_energy/seed_{seed}"
    samples = list(range(0, nsteps + 1, sample_every))
    field_times = times(sorted(set(samples + [s + 1 for s in samples])))
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=nsteps + 1,
        cells=list(cells),
        dl=[dl] * len(cells),
        interp_order=1,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        # v_th,e / c = 0.2: the speed of light warning is expected
        strict=False,
    )
    ph.MaxwellianFluidModel(
        bx=waves.uniform(0.0),
        by=waves.uniform(0.0),
        bz=waves.uniform(0.0),
        protons=waves.maxwellian(1, 1.0, 1.0, Ti, ppc=ppc, seed=1 + 10 * seed),
        electrons=waves.maxwellian(-1, mass_e, 1.0, Te, ppc=ppc, seed=2 + 10 * seed),
    )
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=field_times)
    for pop in masses:
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=times(samples), population_name=pop
        )
    Simulator(sim).run().reset()
    return diag_dir, samples


def nodes(diag_dir, quantity, component, time):
    """physical nodes of a component, each once (periodic): {node key: value}"""
    L = np.asarray(cells) * dl
    hier = hierarchy_from(
        h5_filename=f"{diag_dir}/EM_{quantity}.h5", times=f"{time:.10f}"
    )
    values = {}
    for patch in hier.level(0, time).patches:
        pdata = patch.patch_datas[f"{quantity}{component}"]
        mesh = pdata.meshgrid()
        inside = np.all(
            [(x > -1e-6) & (x < length - 1e-6) for x, length in zip(mesh, L)], axis=0
        )
        keys = np.round(2 * np.stack([x[inside] for x in mesh], axis=1) / dl)
        for key, v in zip(keys.astype(int), pdata.dataset[:][inside]):
            values[tuple(key)] = float(v)
    return values


def field_energy(diag_dir, step):
    t0, t1 = times([step, step + 1])
    electric = 0.0
    magnetic = 0.0
    for component in "xyz":
        before = nodes(diag_dir, "E", component, t0)
        after = nodes(diag_dir, "E", component, t1)
        electric += sum(before[key] * after[key] for key in before)
        magnetic += sum(v * v for v in nodes(diag_dir, "B", component, t0).values())
    return electric / c_over_va**2 / 2 * cell_volume, magnetic / 2 * cell_volume


def kinetic_energy(diag_dir, step):
    (time,) = times([step])
    energy = 0.0
    for pop, mass in masses.items():
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/ions_pop_{pop}_domain.h5", times=f"{time:.10f}"
        )
        for patch in hier.level(0, time).patches:
            particles = patch.patch_datas[f"{pop}_domain"].dataset
            v2 = np.sum(particles.v.astype(np.float64) ** 2, axis=1)
            energy += 0.5 * mass * np.sum(particles.weights.ravel() * v2)
    return energy * cell_volume


def measure(seed):
    diag_dir, samples = run(seed)
    energy = []
    for step in samples:
        kinetic = kinetic_energy(diag_dir, step)
        electric, magnetic = field_energy(diag_dir, step)
        energy.append(kinetic + electric + magnetic)
        print(
            f"seed {seed} step {step}: kinetic {kinetic:.6e} electric {electric:.3e}"
            f" magnetic {magnetic:.3e} total {energy[-1]:.6e}"
        )
    drift = (energy[-1] - energy[0]) / energy[0]
    heating = (energy[-1] - energy[1]) / energy[0]
    print(
        f"seed {seed}: W({nsteps} steps) / W(0) - 1 = {drift:+.3e},"
        f" from step {samples[1]}: {heating:+.3e}"
    )
    return drift, heating


def main():
    seeds = [int(arg) for arg in sys.argv[1:] if arg != "--1d"] or [0]
    drifts = np.array([measure(seed) for seed in seeds])
    for name, values in (
        ("total", drifts[:, 0]),
        (f"from step {sample_every}", drifts[:, 1]),
    ):
        print(
            f"relative energy drift {name}, {len(seeds)} seeds: mean {values.mean():+.3e},"
            f" std {values.std():.3e}, mean + 3 sigma {values.mean() + 3 * values.std():+.3e}"
        )


if __name__ == "__main__":
    main()
