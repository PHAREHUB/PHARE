#!/usr/bin/env python3
"""
FullyKineticModel with one static refined box (P4 M0, tests T1 and T2), periodic root, the
plasma of test_fully_kinetic.py (uniform thermal ions and electrons, a small B and electron
velocity wave) and its command line: dimension, interpolation order, --filter-passes N,
--cleaning s (given explicitly, 0 included: required with refinement). Cases:
    1 1 --cleaning 0                         1D, 64 cells, box over root cells 24 to 39
    2 1 --filter-passes 1 --cleaning 0.5     2D, 32^2 cells, box over root cells 12 to 19
--one-level runs the same case without the box: the reference the bounds below come from.

Checks:
- no regrid: SAMRAI keeps a level whose user boxes did not change (GriddingAlgorithm.cpp:
  2916-2937) unless a data-dependent load balancer is set (refused for FK): the regrid count of
  every level stays 0;
- T1, fine time step: the fine level takes dt/2 (ratio 2, the light CFL is the same on both
  levels), so with a fine dump on every substep its dump times are k dt/2 (within 1e-12, SAMRAI
  evens out the substeps);
- E finite on both levels;
- without cleaning, Gauss on the root away from the box: the residual div E - c^2 F^N rho drifts
  by less than MAX_ROOT_DRIFT (relative to c^2 max |rho|) on root nodes more than ROOT_MARGIN
  coarse cells outside the box. Nothing the fine level does breaks it there: the reflux changes
  B by a curl, which has no divergence, and E only on covered nodes. With cleaning the residual
  is diffused, not kept: reported only.

Reported, not asserted yet: the residual inside and around the box, and in 1D without filter
the fine level residual on fine nodes more than FINE_MARGIN fine cells inside the box.

Known hole until D5 (P4 plan, level-ghost current): level-ghost particles deposit no current, so
in the band next to the fine level boundary the start-up current J_back of E^{-dt_f/2}, and every
step's J, miss the plasma beyond the boundary. The fine level E^{-dt_f/2} is built as the root
one (E^0 = 0, then half a fine step back with the split particles' current), so at t=0 the fine
residual away from that band is that of the split particles.

Numbers: TIME_STEP_NBR = 20 coarse steps gives 10 regrid attempts (every 2 coarse steps with
refined time stepping, TimeRefinementIntegrator.cpp:124-146), 20 swaps of the level-ghost
particles, and moves a thermal electron by 0.9 coarse cell (v_th dt / dx = 0.043 per step), so
the particles within a cell or two of the level boundary cross it. MAX_ROOT_DRIFT is twice the
drift measured on one level for the same case (--one-level), the floor of the single precision
E dumps; ROOT_MARGIN covers the particle shape and the coarsened E stencils around the box.
"""

import sys
import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core.box import Box1D, Box2D
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

one_level = "--one-level" in sys.argv
if one_level:
    sys.argv.remove("--one-level")

# parses the command line
from tests.simulator import test_fully_kinetic as fk  # noqa: E402

ph.NO_GUI()

ndim = fk.ndim
interp_order = fk.interp_orders[0]
BOXES = {1: Box1D(24, 39), 2: Box2D(12, 19)}
box = BOXES[ndim]
ratio = 2
TIME_STEP_NBR = 20
ROOT_MARGIN = 3  # coarse cells around the box
FINE_MARGIN = 3  # fine cells away from the fine level boundary
# twice the drift measured with --one-level (Debug, at cf2cd16b9), max over the root of
# |F(20 dt) - F(0)| / c^2 max |rho|: 1D order 1 N=0: 2.76e-8; 2D order 1 N=1: 8.36e-8
MAX_ROOT_DRIFT = {1: 6e-8, 2: 2e-7}[ndim]
TIME_TOLERANCE = 1e-12

fk.diag_dir = (
    f"phare_outputs/fully_kinetic_amr_{ndim}d_filter{fk.filter_passes}"
    f"_cleaning{fk.cleaning}" + ("_one_level" if one_level else "")
)


def config():
    sim = ph.Simulation(
        time_step=fk.time_step,
        time_step_nbr=TIME_STEP_NBR,
        cells=[fk.cells] * ndim,
        dl=[fk.dl] * ndim,
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=fk.c_over_va,
        current_filter_passes=fk.filter_passes,
        divergence_cleaning=fk.cleaning,
        refinement_boxes=None if one_level else {"L0": [box]},
        smallest_patch_size=8,
        largest_patch_size=fk.CASES[ndim]["tiled"],
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

    common = {"density": fk.uniform(1.0), "nbr_part_per_cell": fk.CASES[ndim]["ppc"]}
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

    every_step = [step * fk.time_step for step in range(TIME_STEP_NBR + 1)]
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=every_step)
    for pop in ["protons", "electrons"]:
        ph.ParticleDiagnostics(
            quantity="domain",
            write_timestamps=[0.0, TIME_STEP_NBR * fk.time_step],
            population_name=pop,
        )
    return sim


def level_times(level):
    import h5py

    with h5py.File(f"{fk.diag_dir}/EM_E.h5", "r") as f:
        return sorted(float(t) for t in f["t"] if f"pl{level}" in f["t"][t])


def outside_box():
    """root primal nodes more than ROOT_MARGIN cells away from the box, along any direction"""
    nodes = np.indices([fk.cells] * ndim)
    outside = np.zeros([fk.cells] * ndim, dtype=bool)
    for d in range(ndim):
        outside |= (nodes[d] < box.lower[d] - ROOT_MARGIN) | (
            nodes[d] > box.upper[d] + 1 + ROOT_MARGIN
        )
    return outside


def fine_gauss_residual_1d(time):
    """1D, no filter: div E - c^2 rho on the fine primal nodes more than FINE_MARGIN cells
    inside the box, rho from the fine domain particles moved half a fine step back"""
    dl = fk.dl / ratio
    dt = fk.time_step / ratio
    lower, upper = ratio * box.lower[0], ratio * (box.upper[0] + 1)
    nodes = range(lower + FINE_MARGIN, upper - FINE_MARGIN + 1)

    rho = dict.fromkeys(range(lower - 2, upper + 3), 0.0)
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
            Ex[round(x / dl - 0.5)] = value

    return np.array(
        [(Ex[n] - Ex[n - 1]) / dl - fk.c_over_va**2 * rho[n] for n in nodes]
    )


class FullyKineticAMRTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def test_static_box(self):
        ph.global_vars.sim = None
        sim = config()
        run = Simulator(sim).run()
        regrids = [run.regrid_count(level) for level in range(sim.max_nbr_levels)]
        run.reset()
        self.assertEqual(regrids, [0] * sim.max_nbr_levels, "regrids per level")

        if cpp.mpi_rank() != 0:
            return

        if not one_level:
            self._check_fine_time_step()
        self._check_finite()
        self._check_root_gauss()
        if not one_level and ndim == 1 and fk.filter_passes == 0:
            self._report_fine_gauss()

    def _check_fine_time_step(self):
        expected = np.arange(2 * TIME_STEP_NBR + 1) * fk.time_step / ratio
        times = np.asarray(level_times(1))
        print(f"fine level dump times: {len(times)}, expected {len(expected)}")
        self.assertEqual(len(times), len(expected))
        self.assertLess(np.abs(times - expected).max(), TIME_TOLERANCE)

    def _check_finite(self):
        time = TIME_STEP_NBR * fk.time_step
        hier = hierarchy_from(
            h5_filename=f"{fk.diag_dir}/EM_E.h5", times=f"{time:.10f}"
        )
        for level in [0] if one_level else [0, 1]:
            for patch in hier.level(level, time).patches:
                for component in "xyz":
                    data = patch.patch_datas[f"E{component}"].dataset[:]
                    self.assertTrue(np.all(np.isfinite(data)), f"L{level} E{component}")

    def _check_root_gauss(self):
        final_time = TIME_STEP_NBR * fk.time_step
        initial = fk.gauss_residual("", 0.0, interp_order)
        final = fk.gauss_residual("", final_time, interp_order)
        rho = fk.charge_density_at_half_step("", 0.0, interp_order)
        scale = fk.c_over_va**2 * np.abs(rho).max()

        change = np.abs(final - initial) / scale
        outside = outside_box()
        print(
            f"{ndim}D N={fk.filter_passes} s={fk.cleaning}"
            f"{' one level' if one_level else ''}: root Gauss residual,"
            f" initial {np.abs(initial).max() / scale:.3e}, change outside the box"
            f" {change[outside].max():.3e}, around and inside {change[~outside].max():.3e}"
        )
        if fk.cleaning == 0:
            # E diagnostics are single precision
            self.assertLess(change[outside].max(), MAX_ROOT_DRIFT)

    def _report_fine_gauss(self):
        scale = (
            fk.c_over_va**2
            * np.abs(fk.charge_density_at_half_step("", 0.0, interp_order)).max()
        )
        for time in [0.0, TIME_STEP_NBR * fk.time_step]:
            residual = fine_gauss_residual_1d(time)
            print(
                f"fine Gauss residual at t = {time:g}, inside the box:"
                f" max {np.abs(residual).max() / scale:.3e} (relative, not asserted)"
            )


if __name__ == "__main__":
    unittest.main()
