#!/usr/bin/env python3
#
# polytropic electron closure with a non uniform entropy carried across a level
# border: the density and the flow are uniform, Pe is a bump that the flow convects
# from the coarse level into a refined box. s = Pe / n^gamma is then just Pe, moved
# at the flow speed, and its fine level ghosts at the border are interpolated in space
# and in time from the coarse level's entropy.
#
# This only checks that the transport through the border is finite, bounded by the
# initial extrema and brings the bump into the fine level: the first order upwind scheme
# diffuses the bump too much (more on the coarse level) to compare it with a profile
# translated by V t, and the effect of the time interpolation on its own is small
# (V dt / dx ~ 0.04 per coarse step).

import unittest

import numpy as np
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

cells = 100
dl = 0.5
L = cells * dl
V = 2.0
P0 = 0.5
amplitude = 0.5
width = 2.5
x0 = 8.0
gamma = 5.0 / 3.0
time_step = 0.01
final_time = 5.0
out = "phare_outputs/polytropic_border_advection/"


def Pe(x):
    return P0 * (1.0 + amplitude * np.exp(-(((x - x0) / width) ** 2)))


def config(diag_dir, **layout):
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        gamma=gamma,
        largest_patch_size=20,
        smallest_patch_size=20,
        diag_options={"format": "phareh5", "options": {"dir": diag_dir, "mode": "overwrite"}},
        **layout,
    )

    def one(x):
        return 1.0

    def zero(x):
        return 0.0

    def vth(x):
        return 0.3

    ph.MaxwellianFluidModel(
        bx=one,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": one,
            "vbulkx": lambda x: V,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": 200,
            "init": {"seed": 2026},
        },
    )
    ph.ElectronModel(closure="polytropic", Pe=Pe)
    ph.FluidDiagnostics(quantity="electron_pressure", write_timestamps=[0.0, final_time])
    ph.FluidDiagnostics(quantity="charge_density", write_timestamps=[0.0, final_time])
    return sim


def domain_nodes(hier, time, ilvl):
    xs, vs = [], []
    for patch in hier.level(ilvl, time).patches:
        pdata = next(iter(patch.patch_datas.values()))
        g = pdata.ghosts_nbr[0]
        domain = slice(g, len(pdata.x) - g)
        xs.append(np.asarray(pdata.x)[domain])
        vs.append(np.asarray(pdata.dataset[:])[domain])
    x, idx = np.unique(np.round(np.concatenate(xs), 9), return_index=True)
    return x, np.concatenate(vs)[idx]


def entropy(run, time, ilvl):
    x, pe = domain_nodes(run.GetPe(time, all_primal=False), time, ilvl)
    xn, n = domain_nodes(run.GetNi(time), time, ilvl)
    np.testing.assert_allclose(x, xn)
    return x, pe / n**gamma


class PolytropicBorderAdvectionTest(unittest.TestCase):
    def tearDown(self):
        ph.global_vars.sim = None

    def test_bump_crosses_the_level_border(self):
        # the bump starts on the coarse level, left of the box (15 < x < 35), and ends up
        # inside it
        diag_dir = out + "refined_box"
        Simulator(config(diag_dir, refinement_boxes={"L0": {"B0": [(30,), (70,)]}})).run()
        run = Run(diag_dir)

        x, s0 = entropy(run, 0.0, 1)
        self.assertLess(np.max(np.abs(s0 - Pe(x))), 1e-6 * P0)

        x, s = entropy(run, final_time, 1)
        self.assertTrue(np.all(np.isfinite(s)))
        scale = amplitude * P0
        moved = (np.max(s) - P0) / scale
        print(f"bump height {moved:.3f} of its initial one")

        # the bump got into the refined level, and the border made neither an excess nor
        # a deficit: an interpolation in time between two states cannot go beyond them
        self.assertGreater(moved, 0.4)
        self.assertLess(np.max(s), P0 * (1 + amplitude) * 1.02)
        self.assertGreater(np.min(s), P0 * 0.98)


if __name__ == "__main__":
    startMPI()
    unittest.main()
