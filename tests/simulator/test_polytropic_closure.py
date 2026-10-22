#!/usr/bin/env python3
#
# polytropic electron closure in the full solver: a pressure balanced structure
# convected by the ion drift V0,
#
#     n = 1, Vi = (V0, 0, 0), B = (0, 0, Bz(x)), Ti uniform
#     Pe(x) = P0 + dP sin(2 pi x / L),  Bz^2 = B0^2 + 2 (P0 + dP - Pe(x))
#
# The exact solution is the initial state translated at V0, so the entropy
# s = Pe / n^gamma, recovered from the dumped electron pressure and ion density,
# must follow s0(x - V0 t). It is compared to the profiles that a pressure left in
# place (0) or carried three times too far (3 V0, one push per PPC stage) would
# give. Runs on one level and with a refined box, whose borders the structure
# crosses. Kept short: a plasma drifting across the grid slowly grows fluctuations
# in PHARE whatever the electron closure, which would dominate a long run.

import unittest

import numpy as np
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

cells = 200
dl = 0.5
L = cells * dl
V0 = 0.5
P0, dP, B0 = 0.5, 0.2, 1.0
gamma = 5.0 / 3.0
time_step = 0.01
final_time = 4.0  # translation 2 = 4 cells
out = "phare_outputs/polytropic_closure/"


def Pe(x):
    return P0 + dP * np.sin(2 * np.pi * x / L)


def Bz(x):
    return np.sqrt(B0**2 + 2 * (P0 + dP - Pe(x)))


def s0(x):  # n = 1 initially
    return Pe(np.mod(x, L))


def config(diag_dir, **layout):
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        largest_patch_size=20,
        smallest_patch_size=20,
        diag_options={"format": "phareh5", "options": {"dir": diag_dir, "mode": "overwrite"}},
        **layout,
    )

    def one(x):
        return 1.0

    def zero(x):
        return 0.0

    def drift(x):
        return V0

    def vth(x):
        return 0.3

    ph.MaxwellianFluidModel(
        bx=zero,
        by=zero,
        bz=Bz,
        protons={
            "charge": 1,
            "density": one,
            "vbulkx": drift,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": 200,
            "init": {"seed": 2026},
        },
    )
    ph.ElectronModel(closure="polytropic", Pe=Pe, Gamma=gamma)

    timestamps = [0.0, final_time]
    ph.FluidDiagnostics(quantity="electron_pressure", write_timestamps=timestamps)
    ph.FluidDiagnostics(quantity="charge_density", write_timestamps=timestamps)
    return sim


def domain_nodes(hier, time, ilvl):
    """x and values on the domain nodes of level ilvl, overlaps removed"""
    xs, vs = [], []
    for patch in hier.level(ilvl, time).patches:
        pdata = list(patch.patch_datas.values())[0]
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


class PolytropicClosureTest(unittest.TestCase):
    def tearDown(self):
        ph.global_vars.sim = None

    def run_and_check(self, name, ilvl, **layout):
        diag_dir = out + name
        Simulator(config(diag_dir, **layout)).run()
        run = Run(diag_dir)

        x, s = entropy(run, final_time, ilvl)

        def rms(shift):
            return np.sqrt(np.mean((s - s0(x - shift)) ** 2))

        errors = {k: rms(k * V0 * final_time) for k in (0, 1, 3)}
        print(name, {f"{k} V0": f"{e:.2e}" for k, e in errors.items()})

        # the translated profile is matched, and much better than the wrong speeds
        self.assertLess(errors[1], 0.2 * min(errors[0], errors[3]))
        self.assertLess(errors[1], 0.05 * dP)

    def test_one_level(self):
        self.run_and_check("one_level", 0)

    def test_refined_box(self):
        # the structure crosses the box borders; checked on the refined level
        self.run_and_check(
            "refined_box", 1, refinement_boxes={"L0": {"B0": [(40,), (120,)]}}
        )


if __name__ == "__main__":
    startMPI()
    unittest.main()
