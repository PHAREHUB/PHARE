#!/usr/bin/env python3
"""
FullyKineticModel, one static refined box, plasma: no transverse instability (P4, pV's reproducer
waves.py 1 --refinement-box 8,23 --seeds 1, found 11-17).

The 1D Langmuir case of tests/functional/fully_kinetic/waves.py (32 cells, a standing electron
velocity wave v0 sin(kx), T_e = 0.04, T_i = 0.01, no B), order 1, with a box over root cells 8 to
23 and the refined runs' options (one filter pass, Marder s = 0.5, boundary_cleaning_mask on),
and the same run on one level. Seed 0: protons 1, electrons 2 (seed s: 1 + 10 s, 2 + 10 s), 100
particles per cell and species, 1200 coarse steps (t = 2.4).
Check: on each level of the box run, the growth rate of the box region's rms B_z (sqrt of its
mean B_z^2 over the box cells; least squares fit of its log against t in [0.3, 2.4]) minus the
same rate on the one level run is below MAX_GROWTH_RATE = 0.5 per unit time. This is pV's
border characterisation metric ("rms box", build/db/border/summary.py, 11-18), on the box cells:
pV's harness (bq.py) took 2 ghost values per side where the dumps have 4, so its window sits 2
cells left of the box (on the fine level, 2 ghost cells in, the last 2 box cells out). The fields
are the same (this test's seed 0 runs are bitwise pV's); on pV's window this test's fit gives pV's
seed 0 numbers exactly (one level +0.169, box +0.111 / +0.103), on the box cells +0.179 and
+0.080 / +0.081: -0.10 / -0.10 (root / fine) relative to one level.

Why a fitted rate, not the former ratio of maxima (max |B_z| after t = 0.3 over its value at
0.3, 600 steps, bound 2): at seed 0, mask off, that ratio is 2.59 on the fine level, 2.03 times
the one level run's, but its window ends at t = 1.2 on a temporary bump (100 particles per cell,
mask off: +25% amplitude at t ~ 1.2, gone by 2.4, absent with the mask on), not on a growth.
pV's fitted rates (its window), box minus one level, t in [0.3, 2.4], 3 seeds (mean +- seed s.e.):

    particles per cell, mask    fine level       root level
    100, on                     +0.01 +- 0.04    0.00 +- 0.03   (seed 0: -0.07, -0.06)
    100, off                    +0.05 +- 0.13    +0.03 +- 0.13  (worst seed: +0.30, fine)
    400, on                     +0.08 +- 0.04    +0.07 +- 0.05

No border instability at order 1 with one filter pass: the box holds 1.6 to 7 times less fine
scale power than a uniform fine run, and nothing grows past t ~ 1.6.
What the bound catches: the stale fine B ghosts after the reflux (before 4103718fb; red test
ee59c5225) gave a fitted rate of max |B_z| of 4.45 (fine) / 4.20 (root), one level 0.11, at 400
particles per cell (4103718fb's message): 8 times the bound. That is max |B_z|, not this rms;
the rms rate of that regime was not measured (UNVERIFIED). This test's former ratio, 9.2 by
t <= 1.2, alone means ln(9.2) / 0.9 >= 2.5 per unit time, 5 times the bound. No runtime option
brings that regime back (the B ghost refill is not switchable): no red run.
Runtime 133 s, Debug (box and one level runs; nightly).

E and B are dumped at every coarse step on both levels (the mode may flip sign every coarse step),
the domain particles at the start and the end; --diagnose dumps the particles at every coarse step.
"""

import sys
import unittest

import h5py
import numpy as np

import pyphare.pharein as ph
from pyphare.simulator import simulator
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

diagnose = "--diagnose" in sys.argv
if diagnose:
    sys.argv.remove("--diagnose")

ph.NO_GUI()

cells, dl, time_step = 32, 0.04, 0.002
c_over_va, mass_e = 10.0, 0.04
L = cells * dl
BOX = (8, 23)
TIME_STEP_NBR = 1200
PPC = 100
FIT_FROM = 0.3
MAX_GROWTH_RATE = 0.5  # per unit time, box run level minus the one level run
diag_root = "phare_outputs/fully_kinetic_amr_stability"


def uniform(value):
    return lambda x: value + 0.0 * x


def maxwellian(charge, mass, T, vx, seed):
    vth = uniform(np.sqrt(T / mass))
    return {
        "charge": charge,
        "mass": mass,
        "density": uniform(1.0),
        "vbulkx": vx,
        "vbulky": uniform(0.0),
        "vbulkz": uniform(0.0),
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
        "nbr_part_per_cell": PPC,
        "init": {"seed": seed},
    }


def config(diag_dir, refined):
    refinement = (
        {"refinement_boxes": {"L0": {"B0": [(BOX[0],), (BOX[1],)]}}} if refined else {}
    )
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=TIME_STEP_NBR,
        cells=[cells],
        dl=[dl],
        interp_order=1,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=1,
        divergence_cleaning=0.5,
        boundary_cleaning_mask=True,
        **refinement,
        strict=False,  # v_th,e / c = 0.1 with a 0.15 drift: warned regime, as waves.py
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite", "fine_dump_lvl_max": 10},
        },
    )
    k = 2 * np.pi / L
    ph.MaxwellianFluidModel(
        bx=uniform(0.0),
        by=uniform(0.0),
        bz=uniform(0.0),
        protons=maxwellian(1, 1.0, 0.01, uniform(0.0), 1),
        electrons=maxwellian(-1, mass_e, 0.04, lambda x: 0.15 * np.sin(k * x), 2),
    )
    every_step = time_step * np.arange(TIME_STEP_NBR + 1)
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=every_step)
    for pop in ["protons", "electrons"]:
        ph.ParticleDiagnostics(
            quantity="domain",
            write_timestamps=every_step if diagnose else every_step[[0, -1]],
            population_name=pop,
        )
    return sim


def box_cells(level):
    """first and last cell of the box region in the level's index space"""
    return BOX[0] * 2**level, (BOX[1] + 1) * 2**level - 1


def growth_rates(diag_dir):
    """{level: growth rate of the box region's rms B_z}, the least squares slope of its log
    against t, for t >= FIT_FROM (pV's "rms box", build/db/border/summary.py)"""
    sums = {}  # {level: {t: [sum of B_z^2, cells]}}
    with h5py.File(f"{diag_dir}/EM_B.h5", "r") as f:
        for t in f["t"]:
            for level, patches in f["t"][t].items():
                lv = int(level[2:])
                lo, hi = box_cells(lv)
                for p in patches.values():
                    lower, n = int(p.attrs["lower"][0]), int(p.attrs["nbrCells"][0])
                    bz = p["EM_B_z"][()]  # dual in 1D: one value per cell, ghosts both sides
                    ghosts = (len(bz) - n) // 2
                    first, last = max(lo, lower), min(hi, lower + n - 1)
                    if first > last:
                        continue
                    values = bz[ghosts + first - lower : ghosts + last - lower + 1]
                    s = sums.setdefault(lv, {}).setdefault(float(t), [0.0, 0])
                    s[0] += np.sum(values**2)
                    s[1] += len(values)
    rates = {}
    for lv, by_time in sorted(sums.items()):
        times = np.array(sorted(by_time))
        cells = np.array([by_time[t][1] for t in times])
        assert np.all(cells == box_cells(lv)[1] - box_cells(lv)[0] + 1), lv
        rms = np.sqrt(np.array([by_time[t][0] for t in times]) / cells)
        fit = (times >= FIT_FROM) & (rms > 0)
        rates[lv] = np.polyfit(times[fit], np.log(rms[fit]), 1)[0]
    return rates


class FullyKineticAMRStabilityTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None
        simulator.exit_on_exception = True

    def run_case(self, refined):
        diag_dir = f"{diag_root}/{'box' if refined else 'one_level'}"
        ph.global_vars.sim = None
        simulator.exit_on_exception = False  # an electron passing c raises
        Simulator(config(diag_dir, refined)).run().reset()
        return growth_rates(diag_dir)

    def test_transverse_fields_stay_at_noise_level(self):
        one = self.run_case(refined=False)[0]
        box = self.run_case(refined=True)
        print(f"one level: rms B_z growth rate {one:+.3f} per unit time")
        for level, rate in box.items():
            print(
                f"box level {level}: {rate:+.3f}, minus one level {rate - one:+.3f}"
                f" (bound {MAX_GROWTH_RATE})"
            )
        self.assertEqual(sorted(box), [0, 1])
        for level, rate in box.items():
            self.assertLess(rate - one, MAX_GROWTH_RATE, f"level {level}")


if __name__ == "__main__":
    unittest.main()
