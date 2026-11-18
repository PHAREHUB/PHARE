#!/usr/bin/env python3
"""
Fully kinetic, 2D, one level, periodic: every patch holding a node of E or B must hold the same
value for it.

A node on a patch border is physical in every patch around it (the Yee primal nodes of a border,
the corner nodes in 2D in four patches). The solver updates each copy on its own patch and the
ghost fills never overwrite them, so the copies only stay equal if their inputs are equal bit for
bit. In the 2D Harris run (run140) the copies of Ez and Bx drifted apart from t ~ 1.8 and grew
exponentially, while the copies of Ey stayed identical.

Uniform Maxwellian plasma, 32 x 32 cells cut into four 16 x 16 patches meeting at a corner
(smallest and largest patch size 16), the Harris run's dx, c, courant, mass ratio, temperatures,
current filter, no cleaning. Every --check-every steps it prints, per component, the largest
difference between two copies of a shared node over the four patches. It then asserts (i) the
copies agree on every shared node, (ii) unless --tiled-only, the E energy matches a single patch
run of the same seed (the particles do not depend on the tiling, per-cell seeds).

Red on 4 MPI ranks, one patch per rank (Debug, 20 ppc, C++ at 46c8b5d73):
    mpirun -n 4 python3 -u test_fully_kinetic_patches.py --tiled-only   (2000 steps, 90 s)
the copies of Ez, Bx and By drift apart and grow exponentially, rate ~8 per 1/Omega_ci (one run:
Ez 3.4e-8 at step 1200, 6.2e-7 at 1400, 1.2e-5 at 1600, 1.9e-3 at 2000, 0.63 at 2400, then
electrons over dx/dt and the Esirkepov displacement abort at step ~2818, as run140), while the
copies of Ex, Ey and Bz stay identical. The onset changes from run to run on 4 ranks (another
run: Ez 1.5e-8 at step 1600), the same input not giving the same result on 3+ ranks (#1526). On 1 rank no copy differs up to step 3200, no abort;
on 2 ranks single last-bit differences of Ez and By, not growing, up to step 2800.
Cause: the copies of a shared node of J were sums of the same terms in different orders, and a
corner node has four, so they differed in the last bits; the border sum is now followed by a
max over the copies.

The fully kinetic modules have no DataWrangler (live patch data), so the copies are read from
the h5 dumps, which hold each patch with its ghosts in single precision: a difference shows
once it reaches float32 resolution, a few 1e-8 of the field.
"""

import argparse
import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest, parse_cli_args

ph.NO_GUI()

parser = argparse.ArgumentParser(description="fully kinetic shared nodes of E and B")
parser.add_argument("--steps", type=int, default=2000)
parser.add_argument("--check-every", type=int, default=200)
parser.add_argument("--ppc", type=int, default=20)
parser.add_argument(
    "--tiled-only", action="store_true", help="no single patch run (no energy check)"
)
cli = parser.parse_args(parse_cli_args())

ndim = 2
cells = 32
tile = 16  # four patches, one corner inside the domain, three more through the periodicity
dl = 0.05
c_over_va = 16.0
courant = 0.75  # c dt sqrt(2) / dx, as tests/functional/fully_kinetic/harris_2d.py
time_step = np.floor(1e6 * courant * dl / (c_over_va * np.sqrt(2))) / 1e6
mass_e = 1.0 / 25.0
Te, Ti = 0.5 / 6.0, 2.5 / 6.0  # T_i / T_e = 5, as the Harris run
filter_passes = 1

COMPONENTS = [f"{q}{c}" for q in ["E", "B"] for c in ["x", "y", "z"]]
diag_dir = "phare_outputs/fully_kinetic_patches"

# relative E energy difference with the single patch run
MAX_ENERGY_DIFFERENCE = 1e-2


def uniform(value):
    def fn(x, y):
        return value + 0.0 * x

    return fn


def check_steps():
    return sorted(set(range(0, cli.steps + 1, cli.check_every)) | {cli.steps})


def config(patch_size):
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=cli.steps,
        cells=[cells] * ndim,
        dl=[dl] * ndim,
        interp_order=1,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=filter_passes,
        divergence_cleaning=0.0,
        smallest_patch_size=patch_size,
        largest_patch_size=patch_size,
        strict=False,  # dx = 2.8 Debye lengths, as the Harris run
        diag_options={
            "format": "phareh5",
            "options": {"dir": f"{diag_dir}/{patch_size}", "mode": "overwrite"},
        },
    )

    def population(charge, mass, T, seed):
        zero, vth = uniform(0.0), uniform(np.sqrt(T / mass))
        return {
            "charge": charge,
            "mass": mass,
            "density": uniform(1.0),
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": cli.ppc,
            "init": {"seed": seed},
        }

    ph.MaxwellianFluidModel(
        bx=uniform(1.0),
        by=uniform(0.0),
        bz=uniform(0.0),
        protons=population(1, 1.0, Ti, 2026),
        electrons=population(-1, mass_e, Te, 2027),
    )
    timestamps = [step * time_step for step in check_steps()]
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    return sim


def copies(patch_size, component, step):
    """the copies of the domain nodes of one component, over every patch: node keys
    (2 x / dl, 2 y / dl modulo the periodic domain, one integer) and values"""
    time = step * time_step
    hier = hierarchy_from(
        h5_filename=f"{diag_dir}/{patch_size}/EM_{component[0]}.h5", times=f"{time:.10f}"
    )
    keys, values = [], []
    patches = hier.level(0, time).patches
    for patch in patches:
        pdata = patch.patch_datas[component]
        lower = np.asarray(patch.box.lower) * dl
        upper = (np.asarray(patch.box.upper) + 1) * dl
        x, y = np.meshgrid(pdata.x, pdata.y, indexing="ij")
        inside = (
            (x > lower[0] - 1e-6 * dl)
            & (x < upper[0] + 1e-6 * dl)
            & (y > lower[1] - 1e-6 * dl)
            & (y < upper[1] + 1e-6 * dl)
        )
        ix = np.rint(2 * x[inside] / dl).astype(int) % (2 * cells)
        iy = np.rint(2 * y[inside] / dl).astype(int) % (2 * cells)
        keys.append(ix * 2 * cells + iy)
        values.append(np.asarray(pdata.dataset[:], dtype=np.float64)[inside])
    return len(patches), np.concatenate(keys), np.concatenate(values)


def mismatch(patch_size, component, step):
    """largest max - min over the copies of a node, that node (x, y), the field max"""
    _, keys, values = copies(patch_size, component, step)
    order = np.argsort(keys, kind="stable")
    keys, values = keys[order], values[order]
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    spread = np.maximum.reduceat(values, starts) - np.minimum.reduceat(values, starts)
    worst = int(np.argmax(spread))
    ix, iy = divmod(int(keys[starts[worst]]), 2 * cells)
    return spread[worst], (ix / 2, iy / 2), np.abs(values).max()


def electric_energy(patch_size, step):
    """sum of E^2 dx dy / 2 over the nodes of the domain, each node once"""
    energy = 0.0
    for component in ["Ex", "Ey", "Ez"]:
        _, keys, values = copies(patch_size, component, step)
        _, first = np.unique(keys, return_index=True)
        energy += np.sum(values[first] ** 2)
    return 0.5 * energy * dl**2


class FullyKineticPatchesTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None
        super().tearDown()

    def run_with(self, patch_size):
        ph.global_vars.sim = None
        Simulator(config(patch_size)).run().reset()

    def test_shared_nodes_and_single_patch(self):
        print(f"dt = {time_step}, c dt / dx = {c_over_va * time_step / dl:.4f}")
        self.run_with(tile)
        if not cli.tiled_only:
            self.run_with(cells)

        worst = {component: 0.0 for component in COMPONENTS}
        for step in check_steps():
            nbr_patches, keys, _ = copies(tile, "Ez", step)
            self.assertEqual(nbr_patches, (cells // tile) ** ndim)
            # Ez is primal-primal: 4 borders of 32 nodes, minus the 4 corners counted twice
            self.assertEqual(np.unique(keys, return_counts=True)[1].tolist().count(1),
                             cells**2 - 4 * cells + 4)  # fmt: skip
            line = []
            for component in COMPONENTS:
                spread, node, scale = mismatch(tile, component, step)
                worst[component] = max(worst[component], spread / max(scale, 1e-300))
                line.append(f"{component} {spread:.2e}/{scale:.2e} at {node}")
            print(f"step {step}: " + ", ".join(line), flush=True)

        for component, relative in worst.items():
            print(f"{component}: largest copy difference / field max {relative:.2e}")
        if cli.tiled_only:
            return self.assertTrue(all(r == 0.0 for r in worst.values()), worst)

        tiled, single = electric_energy(tile, cli.steps), electric_energy(cells, cli.steps)
        energy_difference = abs(tiled - single) / single
        print(f"E energy at step {cli.steps}: tiled {tiled:.6e}, single patch {single:.6e},"
              f" relative difference {energy_difference:.2e}")  # fmt: skip

        for component, relative in worst.items():
            with self.subTest(component=component):
                self.assertEqual(relative, 0.0, f"{component} copies differ")
        self.assertLess(energy_difference, MAX_ENERGY_DIFFERENCE)


if __name__ == "__main__":
    unittest.main()
