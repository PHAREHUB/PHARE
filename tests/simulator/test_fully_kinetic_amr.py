#!/usr/bin/env python3
"""
FullyKineticModel with one static refined box (P4 M0, tests T1 and T2), 1D, periodic root: the
plasma of test_fully_kinetic.py (uniform thermal ions and electrons, a small B and electron
velocity wave), the root of 64 cells and one refined box over its cells 24 to 39, 20 coarse steps.

Checks:
- the run goes through: the FK level initializer throws on a regrid, so a run without exception
  also says SAMRAI never re-initialized the fine level from an old one (static boxes are kept,
  GriddingAlgorithm.cpp:2916-2937);
- T1, fine time step: the fine level takes dt/2 (ratio 2, the light CFL is the same on both
  levels), so with a fine dump on every substep its dump times are k dt/2 (within 1e-12, SAMRAI
  evens out the substeps);
- E finite on both levels;
- Gauss on the root away from the box: the residual div E - c^2 rho drifts by less than 1e-6 of
  c^2 max |rho| on root nodes more than ROOT_MARGIN coarse cells outside the box, as on one level.

Reported, not asserted yet: the fine level Gauss residual at t=0 and at the end, on fine nodes
more than FINE_MARGIN fine cells inside the box.

Known hole until D5 (P4 plan, level-ghost current): level-ghost particles deposit no current, so
in the band next to the fine level boundary the start-up current J_back of E^{-dt_f/2}, and every
step's J, miss the plasma beyond the boundary; the fine residual is only reported inside the box.
The fine level E^{-dt_f/2} is built as the root one (E^0 = 0, then half a fine step back with the
split particles' current), so at t=0 the fine residual is that of the split particles.
"""

import unittest

import h5py
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core.box import Box1D
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest
from tests.simulator import test_fully_kinetic as fk

ph.NO_GUI()

assert fk.ndim == 1, "1D only for now"

fk.diag_dir = "phare_outputs/fully_kinetic_amr_1d"
interp_order = 1
time_step_nbr = 20
box = Box1D(24, 39)  # root cells, refined into fine cells 48 to 79
ratio = 2
ROOT_MARGIN = 3  # coarse cells: the particle shape and the E stencil around the box
FINE_MARGIN = 3  # fine cells away from the fine level boundary
MAX_GAUSS_DRIFT = 1e-6
TIME_TOLERANCE = 1e-12


def config():
    sim = ph.Simulation(
        time_step=fk.time_step,
        time_step_nbr=time_step_nbr,
        cells=[fk.cells],
        dl=[fk.dl],
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=fk.c_over_va,
        refinement_boxes={"L0": [box]},
        smallest_patch_size=8,
        largest_patch_size=16,
        diag_options={
            "format": "phareh5",
            "options": {
                "dir": fk.diag_dir,
                "mode": "overwrite",
                "fine_dump_lvl_max": 1,
            },
        },
        strict=True,
    )

    def vth(mass):
        return fk.uniform(np.sqrt(fk.T / mass))

    common = {"density": fk.uniform(1.0), "nbr_part_per_cell": fk.CASES[1]["ppc"]}
    ph.MaxwellianFluidModel(
        bx=fk.uniform(1.0),
        by=fk.wave(0.1, 0),
        bz=fk.wave(0.1, 1),
        protons={
            "charge": 1,
            "mass": 1.0,
            **{f"vbulk{k}": fk.uniform(0.0) for k in "xyz"},
            **{f"vth{k}": vth(1.0) for k in "xyz"},
            "init": {"seed": 2026},
            **common,
        },
        electrons={
            "charge": -1,
            "mass": fk.mass_e,
            "vbulkx": fk.wave(0.05, 0),
            "vbulky": fk.wave(0.05, 1),
            "vbulkz": fk.wave(0.05, 2),
            **{f"vth{k}": vth(fk.mass_e) for k in "xyz"},
            "init": {"seed": 2027},
            **common,
        },
    )

    every_step = [step * fk.time_step for step in range(time_step_nbr + 1)]
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=every_step)
    for pop in ["protons", "electrons"]:
        ph.ParticleDiagnostics(
            quantity="domain",
            write_timestamps=[0.0, time_step_nbr * fk.time_step],
            population_name=pop,
        )
    return sim


def level_times(level):
    with h5py.File(f"{fk.diag_dir}/EM_E.h5", "r") as f:
        return sorted(float(t) for t in f["t"] if f"pl{level}" in f["t"][t])


def fine_gauss_residual(time):
    """div E - c^2 rho on the fine primal nodes more than FINE_MARGIN cells inside the box,
    rho from the fine domain particles moved half a fine step back: {node: residual}"""
    dl = fk.dl / ratio
    dt = fk.time_step / ratio
    lower, upper = ratio * box.lower[0], ratio * (box.upper[0] + 1)
    nodes = range(lower + FINE_MARGIN, upper - FINE_MARGIN + 1)

    rho = {node: 0.0 for node in range(lower - 2, upper + 3)}
    for pop, charge in fk.charges.items():
        hier = hierarchy_from(
            h5_filename=f"{fk.diag_dir}/ions_pop_{pop}_domain.h5", times=f"{time:.10f}"
        )
        for patch in hier.level(1, time).patches:
            particles = patch.patch_datas[f"{pop}_domain"].dataset
            pos = (
                particles.iCells[:, 0]
                + particles.deltas[:, 0].astype(np.float64)
                - 0.5 * dt * particles.v[:, 0].astype(np.float64) / dl
            )
            left = np.floor(pos).astype(int)
            weights = charge * particles.weights.ravel()
            for offset in range(-1, 3):
                node = left + offset
                w = weights * fk.shape(pos - node, interp_order)
                for n, wn in zip(node, w):
                    if n in rho:
                        rho[n] += wn

    # E_x is dual: its value at x = (i + 1/2) dl, keyed by i
    Ex = {}
    hier = hierarchy_from(h5_filename=f"{fk.diag_dir}/EM_E.h5", times=f"{time:.10f}")
    for patch in hier.level(1, time).patches:
        pdata = patch.patch_datas["Ex"]
        for x, value in zip(pdata.x, pdata.dataset[:]):
            Ex[int(round(x / dl - 0.5))] = value

    return {
        node: (Ex[node] - Ex[node - 1]) / dl - fk.c_over_va**2 * rho[node]
        for node in nodes
    }


class FullyKineticAMRTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def test_static_box(self):
        ph.global_vars.sim = None
        Simulator(config()).run().reset()

        if cpp.mpi_rank() != 0:
            return

        self._check_fine_time_step()
        self._check_finite()
        self._check_root_gauss()
        self._report_fine_gauss()

    def _check_fine_time_step(self):
        expected = np.arange(2 * time_step_nbr + 1) * fk.time_step / ratio
        times = np.asarray(level_times(1))
        print(f"fine level dump times: {len(times)}, expected {len(expected)}")
        self.assertEqual(len(times), len(expected))
        self.assertLess(np.abs(times - expected).max(), TIME_TOLERANCE)

    def _check_finite(self):
        time = time_step_nbr * fk.time_step
        hier = hierarchy_from(h5_filename=f"{fk.diag_dir}/EM_E.h5", times=f"{time:.10f}")
        for level in [0, 1]:
            for patch in hier.level(level, time).patches:
                for component in "xyz":
                    data = patch.patch_datas[f"E{component}"].dataset[:]
                    self.assertTrue(np.all(np.isfinite(data)), f"L{level} E{component}")

    def _check_root_gauss(self):
        final_time = time_step_nbr * fk.time_step
        initial = fk.gauss_residual("", 0.0, interp_order)
        final = fk.gauss_residual("", final_time, interp_order)
        rho = fk.charge_density_at_half_step("", 0.0, interp_order)
        scale = fk.c_over_va**2 * np.abs(rho).max()

        nodes = np.arange(fk.cells)
        outside = (nodes < box.lower[0] - ROOT_MARGIN) | (
            nodes > box.upper[0] + 1 + ROOT_MARGIN
        )
        drift = np.abs(final - initial)[outside].max() / scale
        inside = np.abs(final - initial)[~outside].max() / scale
        print(
            f"root Gauss residual drift: outside the box {drift:.3e},"
            f" inside (reported) {inside:.3e}"
        )
        # E diagnostics are single precision
        self.assertLess(drift, MAX_GAUSS_DRIFT)

    def _report_fine_gauss(self):
        scale = fk.c_over_va**2 * np.abs(
            fk.charge_density_at_half_step("", 0.0, interp_order)
        ).max()
        for time in [0.0, time_step_nbr * fk.time_step]:
            residual = np.array(list(fine_gauss_residual(time).values()))
            print(
                f"fine Gauss residual at t = {time:g}, inside the box:"
                f" max {np.abs(residual).max() / scale:.3e} (relative, not asserted)"
            )


if __name__ == "__main__":
    unittest.main()
