#!/usr/bin/env python3
"""
FullyKineticModel end to end, one level, periodic: kinetic ions and electrons, E from
Ampere-Maxwell, Esirkepov current. Command line (--help): dimension (1 if omitted), the
interpolation orders (those of CASES if omitted), --filter-passes N for the binomial current
filter: the Python charge density of the Gauss check is then filtered the same way, since the
filter commutes with the Yee divergence and the conserved residual is div E - c^2 F^N rho.
--cleaning s for the Marder divergence cleaning: the residual F = div E - c^2 F^N rho is then
no longer conserved but diffused, each step multiplies its mode k by
g_k = 1 - s lambda_k / sum_i (2 / dx_i^2), lambda_k = sum_i 4 / dx_i^2 sin^2(k_i dx_i / 2),
exactly (the particles and the current keep F as it is between two cleanings). The check then
asks for F at step n to be the initial F with every mode multiplied by g_k^n, for its norm to
decrease at every dump, below a measured fraction at the end.

With a seed, the initial particles do not depend on the tiling, so the same run on one patch
and on several patches must give the same fields: this exercises the patch border sums of the
current and charge (in 2D and 3D, at nodes shared by 4 and 8 patches), the field ghosts and the
particle exchange. Measured, serial and on 3 MPI ranks, in 1D, 2D and 3D: no difference in
the single precision dumps. Summation order differs (one patch adds the contributions of a
node in particle order, several patches partial sums then border sums), so the check keeps
headroom rather than asking for equality.
"""

import argparse
import itertools
import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest, parse_cli_args

ph.NO_GUI()

parser = argparse.ArgumentParser(description="FullyKineticModel end to end")
parser.add_argument("ndim", nargs="?", type=int, choices=[1, 2, 3], default=1)
parser.add_argument(
    "orders",
    nargs="*",
    type=int,
    choices=[1, 2, 3],
    help="interpolation orders (those of CASES if omitted)",
)
parser.add_argument(
    "--filter-passes",
    type=int,
    default=0,
    help="binomial current filter passes (the Python charge density gets the same)",
)
parser.add_argument(
    "--cleaning",
    type=float,
    default=0.0,
    help="Marder divergence cleaning s in [0, 1] (Simulation divergence_cleaning)",
)
cli = parser.parse_args(parse_cli_args())
ndim = cli.ndim
filter_passes = cli.filter_passes
cleaning = cli.cleaning

# cells per direction, the patch size of the tiled run (4, 16 and 8 patches: in 2D and 3D,
# nodes shared by 4 and 8 patches), particles per cell and species, interpolation orders,
# time steps. 3D is cut for the PR tier: 4 ppc and 40 steps took 26.8 s at order 1 (Debug), now
# 11.7 s; orders 2 and 3 add 29 s (registered apart)
# Tiles of 8 at orders 2 and 3 are serial only: on several ranks SAMRAI aborts cutting them (field
# ghosts 6, P4 plan pitfall 26); tests/functional/fully_kinetic/compare_runs.py tiles them by 12.
CASES = {
    1: {"cells": 64, "tiled": 16, "ppc": 50, "orders": [1, 2, 3], "steps": 40},
    2: {"cells": 32, "tiled": 8, "ppc": 8, "orders": [1, 2, 3], "steps": 40},
    3: {"cells": 16, "tiled": 8, "ppc": 2, "orders": [1, 2, 3], "steps": 20},
}
cells = CASES[ndim]["cells"]
interp_orders = cli.orders or CASES[ndim]["orders"]
dl = 0.04
L = cells * dl
c_over_va = 10.0
mass_e = 0.04  # m_i / m_e = 25
T = 0.03  # dx = 2.3 Debye lengths, 3 v_th,e / c = 0.27: no setup warning under strict
time_step = 0.002  # omega_pe dt = 0.1, c dt / dx = 0.5
time_step_nbr = CASES[ndim]["steps"]
diag_dir = f"phare_outputs/fully_kinetic_{ndim}d_filter{filter_passes}" + (
    f"_cleaning{cleaning}" if cleaning else ""
)
# with cleaning, E and particles are dumped every SAMPLE_EVERY steps to follow the residual
SAMPLE_EVERY = time_step_nbr // 5

# relative difference of the tiled and single patch fields, and Gauss residual drift.
# Measured (Debug): 0 everywhere; drift 2.3e-8 to 3.5e-8 (1D), 3.2e-8 to 5.5e-8 (2D), 3.2e-8 to
# 4.5e-8 (3D), C++ at a15c436cc. Current filter, C++ at a7e231319, order 1: 0 everywhere, drift
# 8.6e-8 (2D, N=1), 1.2e-7 (2D, N=2), 8.2e-8 (3D, N=1), serial and on 3 ranks
MAX_TILING_DIFFERENCE = 1e-6
MAX_GAUSS_DRIFT = 1e-6
# with cleaning: |F(n) - G^n F(0)| / max |F(0)|, and the final over the initial norm of F.
# Measured, s = 0.5, C++ at c757493ea, Debug: error <= 3.3e-8 (1D), 3.3e-8 (2D), 1.2e-7 (2D,
# N=1), 3.0e-8 (3D order 1, not registered); norm after 40 steps 0.026, 0.034, 0.043 (1D orders
# 1-3), 0.019 (2D), 0.071 (2D, N=1: the filter already took most of the grid scale residual),
# 0.024 after 20 steps (3D); slowest mode g^n to 6 digits
MAX_CLEANING_ERROR = 1e-6
MAX_CLEANED_FRACTION = 0.1


def uniform(value):
    def fn(*xyz):
        return value + 0.0 * xyz[0]

    return fn


def wave(amplitude, direction):
    """amplitude sin(2 pi x_d / L), zero along a direction that is not simulated"""

    def fn(*xyz):
        if direction >= ndim:
            return 0.0 * xyz[0]
        return amplitude * np.sin(2 * np.pi * xyz[direction] / L)

    return fn


def config(largest_patch_size, diag_subdir, interp_order):
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=[cells] * ndim,
        dl=[dl] * ndim,
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=filter_passes,
        divergence_cleaning=cleaning,
        smallest_patch_size=8,
        largest_patch_size=largest_patch_size,
        diag_options={
            "format": "phareh5",
            "options": {"dir": f"{diag_dir}/{diag_subdir}", "mode": "overwrite"},
        },
        strict=True,
    )

    def vth(mass):
        return uniform(np.sqrt(T / mass))

    common = {
        "density": uniform(1.0),
        "nbr_part_per_cell": CASES[ndim]["ppc"],
    }
    # by and bz only vary along x and y: div B = 0. The electron velocity varies along every
    # simulated direction, a longitudinal current in each.
    ph.MaxwellianFluidModel(
        bx=uniform(1.0),
        by=wave(0.1, 0),
        bz=wave(0.1, 1),
        protons={
            "charge": 1,
            "mass": 1.0,
            "vbulkx": uniform(0.0),
            "vbulky": uniform(0.0),
            "vbulkz": uniform(0.0),
            "vthx": vth(1.0),
            "vthy": vth(1.0),
            "vthz": vth(1.0),
            "init": {"seed": 2026},
            **common,
        },
        electrons={
            "charge": -1,
            "mass": mass_e,
            "vbulkx": wave(0.05, 0),
            "vbulky": wave(0.05, 1),
            "vbulkz": wave(0.05, 2),
            "vthx": vth(mass_e),
            "vthy": vth(mass_e),
            "vthz": vth(mass_e),
            "init": {"seed": 2027},
            **common,
        },
    )

    timestamps = sample_times()
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for pop in ["protons", "electrons"]:
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=timestamps, population_name=pop
        )
    return sim


def sample_steps():
    if cleaning:
        return list(range(0, time_step_nbr + 1, SAMPLE_EVERY))
    return [0, time_step_nbr]


def sample_times():
    return [step * time_step for step in sample_steps()]


charges = {"protons": 1.0, "electrons": -1.0}


def shape(distance, interp_order):
    """b-spline of the particle shape on primal nodes, as core::Weighter"""
    d = np.abs(distance)
    if interp_order == 1:
        return np.maximum(0.0, 1.0 - d)
    if interp_order == 2:
        return np.where(
            d < 0.5, 0.75 - d**2, np.where(d < 1.5, 0.5 * (1.5 - d) ** 2, 0.0)
        )
    return np.where(
        d < 1.0,
        2.0 / 3.0 - d**2 + 0.5 * d**3,
        np.where(d < 2.0, (2.0 - d) ** 3 / 6.0, 0.0),
    )


def charge_density_at_half_step(diag_subdir, time, interp_order):
    """rho at x - (dt/2) v of the dumped particles, on the primal nodes 0..cells-1 of each
    direction (periodic): the charge density at the time of the dumped E, which is half a step
    behind. The shape is the product of the 1D shapes."""
    rho = np.zeros([cells] * ndim)
    for pop, charge in charges.items():
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/{diag_subdir}/ions_pop_{pop}_domain.h5",
            times=f"{time:.10f}",
        )
        for patch in hier.level(0, time).patches:
            particles = patch.patch_datas[f"{pop}_domain"].dataset
            pos = (
                particles.iCells
                + particles.deltas.astype(np.float64)
                - 0.5 * time_step * particles.v[:, :ndim].astype(np.float64) / dl
            )
            left = np.floor(pos).astype(int)
            weights = charge * particles.weights.ravel()
            for offsets in itertools.product(range(-1, 3), repeat=ndim):
                node = left + np.asarray(offsets)
                w = weights * np.prod(shape(pos - node, interp_order), axis=1)
                np.add.at(rho, tuple((node % cells).T), w)
    return binomial(rho, filter_passes)


def binomial(field, passes):
    """the current filter on primal periodic nodes: (1/4, 1/2, 1/4) along every direction"""
    for _ in range(passes):
        for axis in range(ndim):
            field = (
                0.25 * np.roll(field, 1, axis)
                + 0.5 * field
                + 0.25 * np.roll(field, -1, axis)
            )
    return field


def field_by_node(diag_subdir, quantity, component, time=time_step * time_step_nbr):
    """physical nodes of one component, {(2 x / dl, 2 y / dl, ...) mod 2 cells: value}"""
    hier = hierarchy_from(
        h5_filename=f"{diag_dir}/{diag_subdir}/EM_{quantity}.h5", times=f"{time:.10f}"
    )
    values = {}
    for patch in hier.level(0, time).patches:
        pdata = patch.patch_datas[f"{quantity}{component}"]
        axes = [getattr(pdata, axis) for axis in "xyz"[:ndim]]
        data = pdata.dataset[:]
        for index in itertools.product(*[range(len(axis)) for axis in axes]):
            position = [axis[i] for axis, i in zip(axes, index)]
            if all(-1e-6 < p < L + 1e-6 for p in position):
                key = tuple(round(2 * p / dl) % (2 * cells) for p in position)
                values[key] = data[index]
    return values


def gauss_residual(diag_subdir, time, interp_order):
    """div E - c^2 rho on the primal nodes, periodic"""
    div = np.zeros([cells] * ndim)
    nodes = np.indices([cells] * ndim).reshape(ndim, -1).T
    for d, component in enumerate("xyz"[:ndim]):
        Ed = field_by_node(diag_subdir, "E", component, time)
        # E_d is dual along d: its node i + 1/2 has key 2i + 1
        e = np.zeros([cells] * ndim)
        for node in nodes:
            key = tuple(2 * n + (1 if dd == d else 0) for dd, n in enumerate(node))
            e[tuple(node)] = Ed[key]
        div += (e - np.roll(e, 1, axis=d)) / dl
    return div - c_over_va**2 * charge_density_at_half_step(
        diag_subdir, time, interp_order
    )


class FullyKineticTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_with(self, largest_patch_size, diag_subdir, interp_order):
        ph.global_vars.sim = None
        Simulator(config(largest_patch_size, diag_subdir, interp_order)).run().reset()

    def test_tiling_and_gauss(self):
        for interp_order in interp_orders:
            with self.subTest(interp_order=interp_order):
                self._tiling_and_gauss(interp_order)

    def _tiling_and_gauss(self, interp_order):
        one_patch, tiled = f"one_patch_{interp_order}", f"tiled_{interp_order}"
        self.run_with(cells, one_patch, interp_order)
        self.run_with(CASES[ndim]["tiled"], tiled, interp_order)

        if cpp.mpi_rank() != 0:
            return

        self._check_tiling(one_patch, tiled, interp_order)
        if cleaning:
            self._check_cleaning(tiled, interp_order)
        else:
            self._check_gauss(tiled, interp_order)

    def _check_tiling(self, one_patch, tiled, interp_order):
        nbr_patches = len(
            hierarchy_from(
                h5_filename=f"{diag_dir}/{tiled}/EM_B.h5", times=f"{0.0:.10f}"
            )
            .level(0, 0.0)
            .patches
        )
        # at least the requested tiling: with MPI the box splitting can cut more patches
        # (1D, 3 ranks: 6), which only adds borders. The single patch run is not one patch
        # either then, the check still compares two tilings.
        self.assertGreaterEqual(nbr_patches, (cells // CASES[ndim]["tiled"]) ** ndim)
        print(f"{ndim}D order {interp_order}: {nbr_patches} patches in the tiled run")

        for quantity in ["E", "B"]:
            for component in ["x", "y", "z"]:
                one = field_by_node(one_patch, quantity, component)
                many = field_by_node(tiled, quantity, component)
                self.assertEqual(one.keys(), many.keys())

                a = np.array([one[k] for k in sorted(one)])
                b = np.array([many[k] for k in sorted(one)])
                self.assertTrue(np.all(np.isfinite(a)), f"{quantity}{component}")
                scale = max(np.abs(a).max(), 1e-12)
                diff = np.abs(a - b).max() / scale
                print(
                    f"{ndim}D order {interp_order} {quantity}{component}:"
                    f" max |one - tiled| / max |one| = {diff:.3e}"
                )
                self.assertLess(diff, MAX_TILING_DIFFERENCE, f"{quantity}{component}")

    def _check_gauss(self, subdir, interp_order):
        initial = gauss_residual(subdir, 0.0, interp_order)
        final = gauss_residual(subdir, time_step * time_step_nbr, interp_order)
        rho = charge_density_at_half_step(subdir, 0.0, interp_order)
        scale = c_over_va**2 * np.abs(rho).max()
        drift = np.abs(final - initial).max() / scale
        print(
            f"{ndim}D order {interp_order} Gauss residual: initial"
            f" {np.abs(initial).max() / scale:.3e}, drift {drift:.3e}"
        )
        # E diagnostics are single precision
        self.assertLess(drift, MAX_GAUSS_DRIFT)

    def _check_cleaning(self, subdir, interp_order):
        steps = sample_steps()
        residuals = [
            gauss_residual(subdir, step * time_step, interp_order) for step in steps
        ]
        initial = residuals[0]
        scale = np.abs(initial).max()
        norms = np.array([np.linalg.norm(F) for F in residuals])
        amplification = marder_amplification()
        # the slowest mode, along x
        slowest = (1,) + (0,) * (ndim - 1)
        initial_modes = np.fft.fftn(initial)
        for step, F, norm in zip(steps, residuals, norms):
            predicted = np.fft.ifftn(amplification**step * initial_modes).real
            error = np.abs(F - predicted).max() / scale
            ratio = np.fft.fftn(F)[slowest] / initial_modes[slowest]
            print(
                f"{ndim}D order {interp_order} step {step}: |F| / |F(0)| ="
                f" {norm / norms[0]:.4e}, |F - G^n F(0)| / max |F(0)| = {error:.3e},"
                f" slowest mode {ratio.real:.6f} (g^n {amplification[slowest] ** step:.6f})"
            )
            self.assertLess(error, MAX_CLEANING_ERROR)
            self.assertAlmostEqual(
                ratio.real, amplification[slowest] ** step, delta=1e-4
            )
        self.assertTrue(np.all(np.diff(norms) <= 1e-6 * norms[0]), norms)
        self.assertLess(norms[-1] / norms[0], MAX_CLEANED_FRACTION)


def marder_amplification():
    """g_k per periodic mode of the grid (numpy fft order), s = d dt sum_i 2 / dx_i^2"""
    modes = np.meshgrid(*[np.fft.fftfreq(cells) * cells] * ndim, indexing="ij")
    eigenvalue = sum(4 / dl**2 * np.sin(np.pi * m / cells) ** 2 for m in modes)
    return 1 - cleaning * eigenvalue / (ndim * 2 / dl**2)


if __name__ == "__main__":
    unittest.main()
