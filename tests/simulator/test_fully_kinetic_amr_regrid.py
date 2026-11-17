#!/usr/bin/env python3
"""
FullyKineticModel, tagging on a moving B, two levels (P5a, T15). A B_z bump of a thin plasma
(n = 0.01: the plasma frequency is far below c k at the bump scale) starts with E = 0
(FKLevelInitializer), so it splits into two light pulses moving at +-c, half a coarse cell per
coarse step at c dt / dx = 0.5. SAMRAI regrids level 1 every two coarse steps and its boxes follow
the pulses. The command line is test_fully_kinetic.py's (dimension, interpolation order,
--filter-passes N, --cleaning s), the plasma its temperature and particles per cell. Cases:
    1 1 --cleaning 0.5                       1D, 64 cells, 24 coarse steps
    2 1 --filter-passes 1 --cleaning 0.5     2D, 32^2 cells, bump uniform in y, 12 coarse steps

Checks:
- the run goes to the end (no SAMRAI assertion in Debug, no exception);
- level 1 regridded at least MIN_REGRIDS times, and at least MIN_REGRIDS + 1 different sets of
  level 1 boxes in the dumps (every coarse step);
- E and B finite on every level of every dump;
- Gauss on the root away from every level 1 box of the run: on root nodes more than ROOT_MARGIN
  coarse cells outside all of them, the drift max |F(T) - F(0)| / c^2 max |rho| (F = div E -
  c^2 S^N rho) is at most twice the drift of the same run on one level on the same nodes, plus
  DRIFT_FLOOR (the float32 E dumps);
- 1D, no filter: the fine level residual max |F| / c^2 max |rho| on fine nodes more than
  FINE_MARGIN fine cells inside a level 1 box, over the dumps: reported, asserted below
  FINE_RESIDUAL_BOUND once measured (None: to be pinned from the first run).
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

# parses the command line
from tests.simulator import test_fully_kinetic as fk  # noqa: E402

ph.NO_GUI()

ndim = fk.ndim
interp_order = fk.interp_orders[0]
ratio = 2
STEPS = {1: 24, 2: 12}[ndim]
DENSITY = 0.01
BUMP, WIDTH = 0.4, 4 * fk.dl  # B_z amplitude and width: each pulse keeps half the amplitude
TAGGING_THRESHOLD = 0.02  # the pulses' 5-cell average slope over 1 + |B| is about 0.036
MIN_REGRIDS = 3
ROOT_MARGIN = 3  # coarse cells around every level 1 box of the run
FINE_MARGIN = 3  # fine cells inside a level 1 box
DRIFT_FLOOR = 1e-7
FINE_RESIDUAL_BOUND = None  # pinned from the first measurement

fk.diag_dir = (
    f"phare_outputs/fully_kinetic_amr_regrid_{ndim}d_o{interp_order}"
    f"_filter{fk.filter_passes}_cleaning{fk.cleaning}"
)
times = [step * fk.time_step for step in range(STEPS + 1)]


def config(refined):
    sim = ph.Simulation(
        time_step=fk.time_step,
        time_step_nbr=STEPS,
        cells=[fk.cells] * ndim,
        dl=[fk.dl] * ndim,
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=fk.c_over_va,
        current_filter_passes=fk.filter_passes,
        divergence_cleaning=fk.cleaning,
        **(
            {
                "refinement": "tagging",
                "max_nbr_levels": 2,
                "tagging_threshold": TAGGING_THRESHOLD,
            }
            if refined
            else {}
        ),
        diag_options={
            "format": "phareh5",
            "options": {
                "dir": f"{fk.diag_dir}/{'refined' if refined else 'one_level'}",
                "mode": "overwrite",
                "fine_dump_lvl_max": 1,
            },
        },
    )

    def vth(mass):
        return fk.uniform(np.sqrt(fk.T / mass))

    def bump(*xyz):
        return BUMP * np.exp(-(((xyz[0] - 0.5 * fk.L) / WIDTH) ** 2))

    def species(charge, mass, seed):
        return {
            "charge": charge,
            "mass": mass,
            "density": fk.uniform(DENSITY),
            **{f"vbulk{k}": fk.uniform(0.0) for k in "xyz"},
            **{f"vth{k}": vth(mass) for k in "xyz"},
            "nbr_part_per_cell": fk.CASES[ndim]["ppc"],
            "init": {"seed": seed},
        }

    ph.MaxwellianFluidModel(
        bx=fk.uniform(0.0),
        by=fk.uniform(0.0),
        bz=bump,
        protons=species(1, 1.0, 2026),
        electrons=species(-1, fk.mass_e, 2027),
    )
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=times)
    for pop in ["protons", "electrons"]:
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=times, population_name=pop
        )
    return sim


def level_boxes(subdir, time):
    """the level 1 patch boxes at that dump, in fine cells, as ((lower), (upper)) tuples"""
    hier = hierarchy_from(
        h5_filename=f"{fk.diag_dir}/{subdir}/EM_B.h5", times=f"{time:.10f}"
    )
    if 1 not in hier.levels(time):
        return []
    return sorted(
        (tuple(p.box.lower), tuple(p.box.upper)) for p in hier.level(1, time).patches
    )


def outside_every_box(all_boxes):
    """root primal nodes more than ROOT_MARGIN coarse cells away from every level 1 box"""
    nodes = np.indices([fk.cells] * ndim)
    inside = np.zeros([fk.cells] * ndim, dtype=bool)
    for lower, upper in all_boxes:
        near = np.ones([fk.cells] * ndim, dtype=bool)
        for d in range(ndim):
            lo = lower[d] // ratio - ROOT_MARGIN
            hi = (upper[d] + 1) // ratio + ROOT_MARGIN
            # periodic: distance along d modulo the domain
            near &= ((nodes[d] - lo) % fk.cells) <= (hi - lo)
        inside |= near
    return ~inside


def fine_residual_1d(subdir, time):
    """1D, no filter: max |div E - c^2 rho| on fine primal nodes more than FINE_MARGIN cells
    inside a level 1 box, rho from the fine domain particles moved half a fine step back, and
    c^2 max |rho| there"""
    dl, dt = fk.dl / ratio, fk.time_step / ratio
    worst, scale = 0.0, 0.0
    efield = hierarchy_from(
        h5_filename=f"{fk.diag_dir}/{subdir}/EM_E.h5", times=f"{time:.10f}"
    )
    particles = {
        pop: hierarchy_from(
            h5_filename=f"{fk.diag_dir}/{subdir}/ions_pop_{pop}_domain.h5",
            times=f"{time:.10f}",
        )
        for pop in fk.charges
    }
    for lower, upper in level_boxes(subdir, time):
        lo, hi = lower[0], upper[0] + 1  # primal nodes of the box
        nodes = range(lo + FINE_MARGIN, hi - FINE_MARGIN + 1)
        rho = dict.fromkeys(range(lo - 3, hi + 4), 0.0)
        for pop, charge in fk.charges.items():
            for patch in particles[pop].level(1, time).patches:
                part = patch.patch_datas[f"{pop}_domain"].dataset
                pos = (
                    part.iCells[:, 0]
                    + part.deltas[:, 0].astype(np.float64)
                    - 0.5 * dt * part.v[:, 0].astype(np.float64) / dl
                )
                left = np.floor(pos).astype(int)
                weights = charge * part.weights.ravel()
                for offset in range(-1, 3):
                    node = left + offset
                    w = weights * fk.shape(pos - node, interp_order)
                    for n, wn in zip(node, w):
                        if n in rho:
                            rho[n] += wn
        Ex = {}
        for patch in efield.level(1, time).patches:
            pdata = patch.patch_datas["Ex"]
            for x, value in zip(pdata.x, pdata.dataset[:]):
                Ex[round(x / dl - 0.5)] = value
        for n in nodes:
            if n in Ex and n - 1 in Ex:
                residual = (Ex[n] - Ex[n - 1]) / dl - fk.c_over_va**2 * rho[n]
                worst = max(worst, abs(residual))
        scale = max(scale, fk.c_over_va**2 * max(abs(v) for v in rho.values()))
    return worst, scale


class FullyKineticAMRRegridTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_case(self, refined):
        ph.global_vars.sim = None
        run = Simulator(config(refined)).run()
        regrids = run.regrid_count(1) if refined else 0
        run.reset()
        return regrids

    def test_moving_pulses(self):
        regrids = self.run_case(refined=True)
        self.run_case(refined=False)
        if cpp.mpi_rank() != 0:
            return

        print(f"{ndim}D order {interp_order}: level 1 regridded {regrids} times")
        self.assertGreaterEqual(regrids, MIN_REGRIDS)
        boxes = [level_boxes("refined", t) for t in times]
        self.assertGreaterEqual(len(set(map(tuple, boxes))), MIN_REGRIDS + 1)

        self._check_finite()
        self._check_root_drift([box for at_time in boxes for box in at_time])
        if ndim == 1 and fk.filter_passes == 0:
            self._check_fine_residual()

    def _check_finite(self):
        for quantity in ["E", "B"]:
            for time in times:
                hier = hierarchy_from(
                    h5_filename=f"{fk.diag_dir}/refined/EM_{quantity}.h5",
                    times=f"{time:.10f}",
                )
                for level in hier.levels(time):
                    for patch in hier.level(level, time).patches:
                        for component in "xyz":
                            data = patch.patch_datas[f"{quantity}{component}"].dataset[:]
                            self.assertTrue(
                                np.all(np.isfinite(data)),
                                f"{quantity}{component} L{level} t={time}",
                            )

    def _check_root_drift(self, all_boxes):
        outside = outside_every_box(all_boxes)
        if not outside.any():
            print("root drift: no root node outside every level 1 box, not checked")
            return
        drift = {}
        for subdir in ["refined", "one_level"]:
            initial = fk.gauss_residual(subdir, 0.0, interp_order)
            final = fk.gauss_residual(subdir, times[-1], interp_order)
            rho = fk.charge_density_at_half_step(subdir, 0.0, interp_order)
            scale = fk.c_over_va**2 * np.abs(rho).max()
            drift[subdir] = np.abs(final - initial)[outside].max() / scale
        print(
            f"root drift on {outside.sum()} nodes outside every level 1 box: refined"
            f" {drift['refined']:.3e}, one level {drift['one_level']:.3e}"
        )
        self.assertLessEqual(drift["refined"], 2 * drift["one_level"] + DRIFT_FLOOR)

    def _check_fine_residual(self):
        worst, scale = 0.0, 0.0
        for time in times:
            w, s = fine_residual_1d("refined", time)
            worst, scale = max(worst, w), max(scale, s)
        relative = worst / scale
        print(f"fine level residual max |F| / c^2 max |rho| = {relative:.3e}")
        if FINE_RESIDUAL_BOUND is not None:
            self.assertLess(relative, FINE_RESIDUAL_BOUND)


if __name__ == "__main__":
    unittest.main()
