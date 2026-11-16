#!/usr/bin/env python3
"""
Energy of a 2D thermal plasma (unmagnetized, kinetic ions and electrons), PHARE units:

    W^n = sum_p m w |v^n|^2 / 2 dV + sum_nodes (E^{n-1/2} . E^{n+1/2} / c^2 + |B^n|^2) / 2 dV

with w = n / ppc the particle weight. In vacuum the leapfrog Yee scheme conserves exactly this
field part (E at half steps, B at integer steps). The dumped E is half a step behind its time
label: the dump at t_n is E^{n-1/2}, so E is dumped at t_n and t_{n+1}.

Setup (default): v_th,e = 0.09 c (dx = 2.2 lambda_De), w_pe dt = 0.1, 32^2 cells, 64 particles
per cell and species, 1000 steps (t = 2), strict setup checks. Measurement only, not a
registered test: prints per seed the relative drift of W, total and after the first sample (E
starts at 0, the field fluctuations take a few plasma periods to fill), and the B_z energy growth
rate per unit time over W(0); then mean, std and mean + 3 sigma. Options: --help.

Measured, seed 0, Debug build; 2D runs with the C++ at a15c436cc, 1D at 524f91e6d or a15c436cc
(the 1D module build time was not recorded). In 2D the grid scale B_z noise grows linearly in
time, driven by the in-plane charge conserving currents; B_x, B_y stay flat, and so does 1D. The
deposit matches the published Esirkepov scheme to 1e-15 (tests, db4b884a8): a property of
unfiltered charge conserving deposition, the current filter (--filter-passes) is the remedy.
- 1D, v_th,e = 0.2 c: +3.06e-4 total, +9.1e-6 after step 100
- 2D, v_th,e = 0.2 c: +1.9e-3 total, +1.6e-3 after step 100, B_z 1.7e-4 -> 1.1e-3
- 2D, default setup: +3.3e-3 total, +2.7e-3 after step 100; B_z energy grows 6.7e-5 per unit
  time (4.2e-4 of W(0)), from 2.0e-5 at t = 0.2 to 1.4e-4 at t = 2; kinetic energy +1.5e-3
- 2D, default setup, --half-dt: +3.5e-3, +2.8e-3; B_z grows 6.3e-5 per unit time
"""

import argparse
from dataclasses import dataclass

import numpy as np
import waves
from waves import c_over_va, dl, mass_e

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

# (|v| + 3 v_th,e) / c = 0.27, below the 0.3 guard of the setup checks
Te = mass_e * 0.9**2
Ti = Te
masses = {"protons": 1.0, "electrons": mass_e}


@dataclass(frozen=True)
class Setup:
    ndim: int = 2
    interp_order: int = 1
    substeps: int = 1  # dt / substeps over the same physical time
    filter_passes: int = 0
    cells_per_direction: int = 32
    ppc: int = 64
    physical_steps: int = 1000  # steps of the reference dt
    sample_every: int = 100  # steps of the reference dt

    @property
    def cells(self):
        return [self.cells_per_direction] * self.ndim

    @property
    def time_step(self):
        return waves.time_step / self.substeps

    @property
    def nsteps(self):
        return self.physical_steps * self.substeps

    @property
    def samples(self):
        return list(range(0, self.nsteps + 1, self.sample_every * self.substeps))

    @property
    def cell_volume(self):
        return dl**self.ndim

    @property
    def variant(self):
        """one output directory per variant: runs never overwrite each other"""
        return (
            f"{self.ndim}d_{self.cells_per_direction}cells_{self.ppc}ppc"
            f"_{self.physical_steps}steps_order{self.interp_order}"
            f"_dt_over_{self.substeps}_filter{self.filter_passes}"
        )

    def times(self, steps):
        return [round(step * self.time_step, 10) for step in steps]


def run(setup, seed):
    diag_dir = f"phare_outputs/fully_kinetic_energy/{setup.variant}/seed_{seed}"
    samples = setup.samples
    field_times = setup.times(sorted(set(samples + [s + 1 for s in samples])))
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=setup.time_step,
        time_step_nbr=setup.nsteps + 1,
        cells=setup.cells,
        dl=[dl] * setup.ndim,
        interp_order=setup.interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=setup.filter_passes,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
    )
    ppc = setup.ppc
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
            quantity="domain",
            write_timestamps=setup.times(samples),
            population_name=pop,
        )
    Simulator(sim).run().reset()
    return diag_dir


def nodes(setup, diag_dir, quantity, component, time):
    """physical nodes of a component, each once (periodic): {node key: value}"""
    L = np.asarray(setup.cells) * dl
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


def field_energy(setup, diag_dir, step):
    """electric (E^{n-1/2} . E^{n+1/2}) and magnetic energies, and the B_z part of the latter"""
    t0, t1 = setup.times([step, step + 1])
    electric = 0.0
    magnetic = {}
    for component in "xyz":
        before = nodes(setup, diag_dir, "E", component, t0)
        after = nodes(setup, diag_dir, "E", component, t1)
        electric += sum(before[key] * after[key] for key in before)
        b = nodes(setup, diag_dir, "B", component, t0).values()
        magnetic[component] = sum(v * v for v in b) / 2 * setup.cell_volume
    electric *= setup.cell_volume / c_over_va**2 / 2
    return electric, sum(magnetic.values()), magnetic["z"]


def kinetic_energy(setup, diag_dir, step):
    (time,) = setup.times([step])
    energy = 0.0
    for pop, mass in masses.items():
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/ions_pop_{pop}_domain.h5", times=f"{time:.10f}"
        )
        for patch in hier.level(0, time).patches:
            particles = patch.patch_datas[f"{pop}_domain"].dataset
            v2 = np.sum(particles.v.astype(np.float64) ** 2, axis=1)
            energy += 0.5 * mass * np.sum(particles.weights.ravel() * v2)
    return energy * setup.cell_volume


def measure(setup, seed):
    """relative drift of W (total, after the first sample) and the B_z energy growth rate per
    unit time over W(0), fitted after the first sample"""
    diag_dir = run(setup, seed)
    samples = setup.samples
    times = np.array(setup.times(samples))
    energy, bz = [], []
    for step in samples:
        kinetic = kinetic_energy(setup, diag_dir, step)
        electric, magnetic, magnetic_z = field_energy(setup, diag_dir, step)
        energy.append(kinetic + electric + magnetic)
        bz.append(magnetic_z)
        print(
            f"seed {seed} step {step}: kinetic {kinetic:.6e} electric {electric:.3e}"
            f" magnetic {magnetic:.3e} (B_z {magnetic_z:.3e}) total {energy[-1]:.6e}"
        )
    drift = (energy[-1] - energy[0]) / energy[0]
    heating = (energy[-1] - energy[1]) / energy[0]
    bz_rate = np.polyfit(times[1:], np.array(bz[1:]) / energy[0], 1)[0]
    print(
        f"seed {seed}: W / W(0) - 1 = {drift:+.3e} total, {heating:+.3e} after step"
        f" {samples[1]}; B_z energy growth / W(0) = {bz_rate:.3e} per unit time"
    )
    return drift, heating, bz_rate


def parse_args():
    parser = argparse.ArgumentParser(description="FullyKineticModel energy measurement")
    parser.add_argument(
        "seeds", nargs="*", type=int, default=[0], help="random seeds, one run each"
    )
    parser.add_argument(
        "--1d",
        dest="one_d",
        action="store_true",
        help="64 cells in 1D, a check of the script itself",
    )
    parser.add_argument(
        "--half-dt",
        action="store_true",
        help="dt / 2 over the same physical time (c dt / dx = 0.25)",
    )
    parser.add_argument(
        "--order", type=int, choices=[1, 2, 3], default=1, help="interpolation order"
    )
    parser.add_argument(
        "--filter-passes",
        type=int,
        default=0,
        help="binomial current filter passes (Simulation current_filter_passes)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    setup = Setup(
        ndim=1 if args.one_d else 2,
        cells_per_direction=64 if args.one_d else 32,
        interp_order=args.order,
        substeps=2 if args.half_dt else 1,
        filter_passes=args.filter_passes,
    )
    results = np.array([measure(setup, seed) for seed in args.seeds])
    for name, values in (
        ("energy drift, total", results[:, 0]),
        (f"energy drift after step {setup.sample_every}", results[:, 1]),
        ("B_z energy growth / W(0) per unit time", results[:, 2]),
    ):
        print(
            f"{name}, {len(args.seeds)} seeds: mean {values.mean():+.3e},"
            f" std {values.std():.3e}, mean + 3 sigma {values.mean() + 3 * values.std():+.3e}"
        )


if __name__ == "__main__":
    main()
