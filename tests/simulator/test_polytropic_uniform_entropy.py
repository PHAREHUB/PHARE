#!/usr/bin/env python3
#
# polytropic electron closure in the full solver, compressive flow: the entropy
# s = Pe / n^gamma, uniform at t=0 (Pe = P0 n(x)^gamma), stays uniform since s is
# only transported and carries no source, whatever the density does. It is
# recovered here from the dumped electron pressure and ion density.
#
# The ions converge on the middle of the box (Vi = V1 sin(2 pi x / L), B uniform),
# so that the density, and Pe with it, change a lot: a test with a uniform density
# would pass with Pe left untouched. A particle density is noisy: if s were
# seeded from the deposited density rather than from the analytic profiles, it
# would not be uniform at all but carry gamma dn/n of that noise, ~10 % here.
# Runs on one level and with a refined box.

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

cells = 100
dl = 0.5
L = cells * dl
n1 = 0.3
V1 = 0.5
P0 = 0.5
gamma = 5.0 / 3.0
time_step = 0.01
final_time = 1.0
out = "phare_outputs/polytropic_uniform_entropy/"


def density(x):
    return 1.0 + n1 * np.sin(2 * np.pi * x / L)


def Pe(x):
    return P0 * density(x) ** gamma


def vx(x):
    return V1 * np.sin(2 * np.pi * x / L)


def bump(x):
    # a localised structure in B: what the tagging refines, and that the compressive flow moves
    return 1.0 * np.exp(-(((x - 0.35 * L) / 3.0) ** 2))


def config(diag_dir, bz=None, timestamps=(0.0, final_time), **layout):
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        gamma=gamma,
        largest_patch_size=20,
        smallest_patch_size=20,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
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
        bz=bz or zero,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": vx,
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

    ph.FluidDiagnostics(quantity="electron_pressure", write_timestamps=list(timestamps))
    ph.FluidDiagnostics(quantity="charge_density", write_timestamps=list(timestamps))
    return sim


def domain_nodes(hier, time, ilvl):
    """values on the domain nodes of level ilvl, overlaps removed"""
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
    return n, pe / n**gamma


def level_boxes(hier, time, ilvl):
    return sorted(
        (tuple(p.box.lower), tuple(p.box.upper)) for p in hier.level(ilvl, time).patches
    )


class PolytropicUniformEntropyTest(unittest.TestCase):
    def tearDown(self):
        ph.global_vars.sim = None

    def run_and_check(self, name, ilvl, **layout):
        diag_dir = out + name
        Simulator(config(diag_dir, **layout)).run()
        run = Run(diag_dir)

        n0, s_init = entropy(run, 0.0, ilvl)
        n, s = entropy(run, final_time, ilvl)

        # the flow does compress: the density, and so Pe, changed
        compression = np.max(np.abs(n - n0) / n0)
        deviation = {
            t: np.max(np.abs(si / P0 - 1)) for t, si in ((0.0, s_init), (final_time, s))
        }
        print(
            name,
            f"max dn/n {compression:.2e}",
            {t: f"{d:.2e}" for t, d in deviation.items()},
        )

        self.assertGreater(compression, 0.05)
        # s is the user's, uniform, to round off: and not that of the noisy density
        for t, d in deviation.items():
            self.assertLess(d, 1e-6, f"entropy not uniform at t={t}")

    def test_one_level(self):
        self.run_and_check("one_level", 0)

    def test_refined_box(self):
        self.run_and_check(
            "refined_box", 1, refinement_boxes={"L0": {"B0": [(30,), (70,)]}}
        )

    def test_regridding(self):
        # the refined level follows a structure convected by the flow: it is regridded, and
        # the entropy of the new fine nodes must be the (uniform) one of the coarse level, not
        # that of a pressure refined then divided by a noisy particle density
        times = (0.0, 0.5 * final_time, final_time)
        diag_dir = out + "regridding"
        Simulator(
            config(
                diag_dir,
                bz=lambda x: 1.0 + bump(x),
                timestamps=times,
                refinement="tagging",
                max_nbr_levels=2,
                tagging_threshold=0.01,
                tag_buffer=2,
                nesting_buffer=1,
            )
        ).run()
        run = Run(diag_dir)

        boxes = [level_boxes(run.GetNi(t), t, 1) for t in times]
        print("level 1 boxes", dict(zip(times, boxes)))
        self.assertTrue(all(len(b) > 0 for b in boxes), "no refined level")
        self.assertNotEqual(
            boxes[0], boxes[-1], "the refined level did not move: no regrid"
        )

        for t in times:
            _, s = entropy(run, t, 1)
            d = np.max(np.abs(s / P0 - 1))
            print(f"regridding t={t}", f"{d:.2e}")
            self.assertLess(d, 1e-6, f"entropy not uniform on level 1 at t={t}")


if __name__ == "__main__":
    startMPI()
    unittest.main()
