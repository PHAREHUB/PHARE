#!/usr/bin/env python3
"""
FullyKineticModel end to end, one level, periodic: kinetic ions and electrons, E from
Ampere-Maxwell, Esirkepov current. The dimension is the first command line argument (1 if
omitted).

With a seed, the initial particles do not depend on the tiling, so the same run on one patch
and on several patches must give the same fields: this exercises the patch border sums of the
current and charge (in 2D and 3D, at nodes shared by 4 and 8 patches), the field ghosts and the
particle exchange. Not bitwise: one patch sums the contributions of a node in particle order,
several patches as partial sums then border sums. The dumps are single precision.
"""

import itertools
import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest, parse_cli_args

ph.NO_GUI()

ndim = int((parse_cli_args() or ["1"])[0])

# cells per direction, the patch size of the tiled run (4, 16 and 8 patches: in 2D and 3D,
# nodes shared by 4 and 8 patches), particles per cell and species, interpolation orders,
# time steps. 3D is cut to stay in the PR tier: 4 ppc and 40 steps took 26.8 s (Debug)
CASES = {
    1: {"cells": 64, "tiled": 16, "ppc": 50, "orders": [1, 2, 3], "steps": 40},
    2: {"cells": 32, "tiled": 8, "ppc": 8, "orders": [1, 2, 3], "steps": 40},
    3: {"cells": 16, "tiled": 8, "ppc": 2, "orders": [1], "steps": 20},
}
cells = CASES[ndim]["cells"]
interp_orders = CASES[ndim]["orders"]
dl = 0.04
L = cells * dl
c_over_va = 10.0
mass_e = 0.04  # m_i / m_e = 25
T = 0.03  # dx = 2.3 Debye lengths, 3 v_th,e / c = 0.27: no setup warning under strict
time_step = 0.002  # omega_pe dt = 0.1, c dt / dx = 0.5
time_step_nbr = CASES[ndim]["steps"]
diag_dir = f"phare_outputs/fully_kinetic_{ndim}d"

# relative difference of the tiled and single patch fields, and Gauss residual drift
# (1D, measured: 0 and 2.3e-8 to 3.5e-8)
MAX_TILING_DIFFERENCE = 1e-6
MAX_GAUSS_DRIFT = 1e-6


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

    timestamps = [0.0, time_step * time_step_nbr]
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for pop in ["protons", "electrons"]:
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=timestamps, population_name=pop
        )
    return sim


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
    return rho


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


if __name__ == "__main__":
    unittest.main()
