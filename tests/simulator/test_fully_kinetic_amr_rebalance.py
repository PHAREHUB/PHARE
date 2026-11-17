#!/usr/bin/env python3
"""
FullyKineticModel, level 0 rebalanced during the run (P5b, T16): T15's moving light pulses (two
levels, level 1 regridded every two coarse steps) with a LoadBalancer that rebalances level 0
every EVERY coarse steps, against the same run without one. 1D, T15's command line:
    1 1 --cleaning 0.5

A rebalance re-creates level 0 (SAMRAI GriddingAlgorithm::makeCoarsestLevel: a new level object,
initialized from the old one) and FKLevelInitializer copies B, E and the particles, nothing is
recomputed. The particles are copied cell by cell (ParticleArray::export_particles, through the
cell map), so their order changes even where the boxes come back identical, and the deposits add
them in another order: the two runs agree to round-off, not bitwise, and the dumps are float32.

Checks, boundary_cleaning_mask off and on (on needs cleaning > 0):
- level 0 re-created (STEPS - 1) // EVERY times, never in the reference run;
- one rank: level 0 boxes identical to the reference's at every dump;
- E and B on every level and at every dump within TOLERANCE of the reference, relative to the
  reference's max (on several ranks the reference is not bitwise reproducible anyway, #1526);
- T15's root Gauss check: the drift of the rebalanced run on the root nodes away from every
  level 1 box at most twice the reference's, plus DRIFT_FLOOR.
With the mask on, a new level 0 object used to forget the cells level 1 synchronized
(SolverFK::skippedNodes_): one cleaning without the covered edge skip after each rebalance.
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

# parses the command line
from tests.simulator import test_fully_kinetic as fk
from tests.simulator import test_fully_kinetic_amr_regrid as t15

ph.NO_GUI()

# coarse steps between two rebalances: odd, so some follow a level 1 regrid and some not
EVERY = 3
# float32 dumps (6e-8) and round-off; pinned from the first measurement
TOLERANCE = 1e-6
MASKS = [False, True] if fk.cleaning > 0 else [False]

fk.diag_dir = (
    f"phare_outputs/fully_kinetic_amr_rebalance_{t15.ndim}d_o{t15.interp_order}"
    f"_filter{fk.filter_passes}_cleaning{fk.cleaning}"
)


def boxes(subdir, time, level):
    """the patch boxes of a level at that dump, as ((lower), (upper)) tuples"""
    hier = hierarchy_from(
        h5_filename=f"{fk.diag_dir}/{subdir}/EM_B.h5", times=f"{time:.10f}"
    )
    if level not in hier.levels(time):
        return []
    return sorted(
        (tuple(p.box.lower), tuple(p.box.upper))
        for p in hier.level(level, time).patches
    )


def level_nodes(subdir, quantity, component, time, level):
    """1D: {2 x / dl_level: value} on the nodes of the level's patch boxes, ghosts excluded"""
    hier = hierarchy_from(
        h5_filename=f"{fk.diag_dir}/{subdir}/EM_{quantity}.h5", times=f"{time:.10f}"
    )
    if level not in hier.levels(time):
        return {}
    dl = fk.dl / t15.ratio**level
    values = {}
    for patch in hier.level(level, time).patches:
        pdata = patch.patch_datas[f"{quantity}{component}"]
        lower, upper = patch.box.lower[0] * dl, (patch.box.upper[0] + 1) * dl
        for x, value in zip(pdata.x, pdata.dataset[:]):
            if lower - 1e-6 <= x <= upper + 1e-6:
                values[round(2 * x / dl)] = value
    return values


class FullyKineticAMRRebalanceTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_case(self, subdir, mask, rebalance):
        ph.global_vars.sim = None
        sim = t15.config(True, subdir=subdir, boundary_cleaning_mask=mask)
        if rebalance:
            ph.LoadBalancer(every=EVERY, on_init=False)
        run = Simulator(sim).initialize()
        for _ in range(t15.STEPS):  # not run(): it resets, the counts go with it
            run.advance()
        rebalances = run.regrid_count(0)
        run.reset()
        return rebalances

    def test_rebalanced_level_0(self):
        if t15.ndim != 1:
            self.skipTest("T16 is 1D")
        for mask in MASKS:
            with self.subTest(mask=mask):
                self._rebalanced(mask)

    def _rebalanced(self, mask):
        reference, rebalanced = (
            f"reference_mask{int(mask)}",
            f"rebalanced_mask{int(mask)}",
        )
        self.assertEqual(self.run_case(reference, mask, rebalance=False), 0)
        rebalances = self.run_case(rebalanced, mask, rebalance=True)
        if cpp.mpi_rank() != 0:
            return

        print(f"mask {mask}: level 0 re-created {rebalances} times")
        self.assertEqual(rebalances, (t15.STEPS - 1) // EVERY)
        if cpp.mpi_size() == 1:
            for time in t15.times:
                self.assertEqual(
                    boxes(reference, time, 0), boxes(rebalanced, time, 0), f"t={time}"
                )
        self._check_fields(reference, rebalanced, mask)
        self._check_root_drift(reference, rebalanced)

    def _check_fields(self, reference, rebalanced, mask):
        worst, where = 0.0, None
        for time in t15.times:
            for level in (0, 1):
                for quantity in ("E", "B"):
                    for component in "xyz":
                        ref = level_nodes(reference, quantity, component, time, level)
                        new = level_nodes(rebalanced, quantity, component, time, level)
                        label = f"{quantity}{component} L{level} t={time:.3f}"
                        self.assertEqual(ref.keys(), new.keys(), label)
                        if not ref:
                            continue
                        a = np.array([ref[k] for k in sorted(ref)])
                        b = np.array([new[k] for k in sorted(ref)])
                        self.assertTrue(np.all(np.isfinite(b)), label)
                        diff = np.abs(b - a).max() / max(np.abs(a).max(), 1e-12)
                        if diff >= worst:
                            worst, where = diff, label
        print(
            f"mask {mask}: max |rebalanced - reference| / max |reference|"
            f" = {worst:.3e} ({where})"
        )
        self.assertLess(worst, TOLERANCE, where)

    def _check_root_drift(self, reference, rebalanced):
        all_boxes = [
            box
            for subdir in (reference, rebalanced)
            for time in t15.times
            for box in boxes(subdir, time, 1)
        ]
        outside = t15.outside_every_box(all_boxes)
        if not outside.any():
            print("root drift: no root node outside every level 1 box, not checked")
            return
        drift = {}
        for subdir in (reference, rebalanced):
            initial = fk.gauss_residual(subdir, 0.0, t15.interp_order)
            final = fk.gauss_residual(subdir, t15.times[-1], t15.interp_order)
            rho = fk.charge_density_at_half_step(subdir, 0.0, t15.interp_order)
            scale = fk.c_over_va**2 * np.abs(rho).max()
            drift[subdir] = np.abs(final - initial)[outside].max() / scale
        print(
            f"root drift on {outside.sum()} nodes: rebalanced {drift[rebalanced]:.3e},"
            f" reference {drift[reference]:.3e}"
        )
        self.assertLessEqual(drift[rebalanced], 2 * drift[reference] + t15.DRIFT_FLOOR)


if __name__ == "__main__":
    unittest.main()
