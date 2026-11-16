#!/usr/bin/env python3
"""
FullyKineticModel end to end, 1D, one level, periodic: kinetic ions and electrons, E from
Ampere-Maxwell, Esirkepov current.

With a seed, the initial particles do not depend on the tiling, so the same run on one patch
and on four patches must give the same fields: this exercises the patch border sums of the
current and charge, the field ghosts and the particle exchange.
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = 64
dl = 0.04
L = cells * dl
c_over_va = 10.0
mass_e = 0.04  # m_i / m_e = 25
T = 0.03  # dx = 2.3 Debye lengths, 3 v_th,e / c = 0.27: no setup warning under strict
time_step = 0.002  # omega_pe dt = 0.1, c dt / dx = 0.5
time_step_nbr = 40
diag_dir = "phare_outputs/fully_kinetic"


interp_orders = [1, 2, 3]


def config(largest_patch_size, diag_subdir, interp_order):
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=dl,
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

    def density(x):
        return 1.0 + 0.0 * x

    def zero(x):
        return 0.0 * x

    def bx(x):
        return 1.0 + 0.0 * x

    def by(x):
        return 0.1 * np.sin(2 * np.pi * x / L)

    def vx_electrons(x):
        return 0.05 * np.sin(2 * np.pi * x / L)

    def vth(mass):
        def fn(x):
            return np.sqrt(T / mass) + 0.0 * x

        return fn

    common = {
        "density": density,
        "vbulky": zero,
        "vbulkz": zero,
        "nbr_part_per_cell": 50,
    }
    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=zero,
        protons={
            "charge": 1,
            "mass": 1.0,
            "vbulkx": zero,
            "vthx": vth(1.0),
            "vthy": vth(1.0),
            "vthz": vth(1.0),
            "init": {"seed": 2026},
            **common,
        },
        electrons={
            "charge": -1,
            "mass": mass_e,
            "vbulkx": vx_electrons,
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
    """rho at x - (dt/2) v of the dumped particles, on the primal nodes 0..cells-1 (periodic):
    the charge density at the time of the dumped E, which is half a step behind"""
    rho = np.zeros(cells)
    for pop, charge in charges.items():
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/{diag_subdir}/ions_pop_{pop}_domain.h5",
            times=f"{time:.10f}",
        )
        for patch in hier.level(0, time).patches:
            particles = patch.patch_datas[f"{pop}_domain"].dataset
            pos = (
                particles.iCells[:, 0]
                + particles.deltas[:, 0].astype(np.float64)
                - 0.5 * time_step * particles.v[:, 0].astype(np.float64) / dl
            )
            left = np.floor(pos).astype(int)
            for node in (left - 1, left, left + 1, left + 2):
                np.add.at(
                    rho,
                    node % cells,
                    charge
                    * particles.weights.ravel()
                    * shape(pos - node, interp_order),
                )
    return rho


def gauss_residual(diag_subdir, time, interp_order):
    """div E - c^2 rho on the primal nodes, periodic"""
    Ex = field_by_node(diag_subdir, "E", "x", time)
    # dual node j at (j + 1/2) dl has key 2j+1
    ex = np.array([Ex[2 * j + 1] for j in range(cells)], dtype=np.float64)
    div = (ex - np.roll(ex, 1)) / dl
    return div - c_over_va**2 * charge_density_at_half_step(
        diag_subdir, time, interp_order
    )


def field_by_node(diag_subdir, quantity, component, time=time_step * time_step_nbr):
    """physical nodes of one component, {2 x / dl: value}"""
    hier = hierarchy_from(
        h5_filename=f"{diag_dir}/{diag_subdir}/EM_{quantity}.h5", times=f"{time:.10f}"
    )
    values = {}
    for patch in hier.level(0, time).patches:
        pdata = patch.patch_datas[f"{quantity}{component}"]
        for x, v in zip(pdata.x, pdata.dataset[:]):
            if -1e-6 < x < L + 1e-6:
                values[round(2 * x / dl)] = v
    return values


class FullyKineticTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_with(self, largest_patch_size, diag_subdir, interp_order):
        ph.global_vars.sim = None
        Simulator(config(largest_patch_size, diag_subdir, interp_order)).run().reset()

    def test_one_patch_and_four_patches_agree(self):
        for interp_order in interp_orders:
            with self.subTest(interp_order=interp_order):
                self._one_and_four_patches(interp_order)

    def _one_and_four_patches(self, interp_order):
        one_patch, four_patches = (
            f"one_patch_{interp_order}",
            f"four_patches_{interp_order}",
        )
        self.run_with(cells, one_patch, interp_order)
        self.run_with(cells // 4, four_patches, interp_order)

        if cpp.mpi_rank() != 0:
            return

        for quantity in ["E", "B"]:
            for component in ["x", "y", "z"]:
                one = field_by_node(one_patch, quantity, component)
                four = field_by_node(four_patches, quantity, component)
                self.assertEqual(one.keys(), four.keys())

                a = np.array([one[k] for k in sorted(one)])
                b = np.array([four[k] for k in sorted(one)])
                self.assertTrue(np.all(np.isfinite(a)), f"{quantity}{component}")
                scale = max(np.abs(a).max(), 1e-12)
                diff = np.abs(a - b).max() / scale
                print(
                    f"order {interp_order} {quantity}{component}:"
                    f" max |one - four| / max |one| = {diff:.3e}"
                )
                self.assertLess(diff, 1e-5, f"{quantity}{component}")

    def test_gauss_residual_does_not_drift(self):
        for interp_order in interp_orders:
            with self.subTest(interp_order=interp_order):
                self._gauss(interp_order)

    def _gauss(self, interp_order):
        subdir = f"gauss_{interp_order}"
        self.run_with(cells // 4, subdir, interp_order)

        if cpp.mpi_rank() != 0:
            return

        initial = gauss_residual(subdir, 0.0, interp_order)
        final = gauss_residual(subdir, time_step * time_step_nbr, interp_order)
        rho = charge_density_at_half_step(subdir, 0.0, interp_order)
        scale = c_over_va**2 * np.abs(rho).max()
        drift = np.abs(final - initial).max() / scale
        print(
            f"order {interp_order} Gauss residual: initial"
            f" {np.abs(initial).max() / scale:.3e}, drift {drift:.3e}"
        )
        # E diagnostics are single precision
        self.assertLess(drift, 1e-6)


if __name__ == "__main__":
    unittest.main()
