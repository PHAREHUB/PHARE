#!/usr/bin/env python3
"""
FullyKineticModel, one static refined box, plasma: no transverse instability (P4, pV's reproducer
waves.py 1 --refinement-box 8,23 --seeds 1, found 11-17).

The 1D Langmuir case of tests/functional/fully_kinetic/waves.py (32 cells, a standing electron
velocity wave v0 sin(kx), T_e = 0.04, T_i = 0.01, no B), order 1, with a box over root cells 8 to
23 and the refined runs' options (one filter pass, Marder s = 0.5), and the same run on one level.
On one level the transverse fields stay at the noise level; with the box and stale fine B ghosts
at the first substep after the reflux, B_z and E_y grew at the coarse grid scale (gamma about 4.4
per unit time, 0.09 w_pe) until an electron passed c.
Checks: on each level of the box run, the growth of max |B_z| after t = 0.3 (its max after over
its value at t = 0.3) is below MAX_RELATIVE_GROWTH times that of the one level run.

Calibration (Debug, C++ at bd07c2150, 100 particles per cell and species, 600 coarse steps,
seed s: protons 1 + 10 s, electrons 2 + 10 s; this test is seed 0), growth then ratio:

    seed   one level   box level 0       box level 1
    0      1.28        1.86 (1.45)       2.83 (2.21)
    1      1.28        1.77 (1.38)       1.90 (1.48)
    2      1.76        1.27 (0.72)       1.35 (0.77)

Stale ghosts (at ee59c5225): 9.2 on the box run, ratio above 7. The seed 0 level 1 ratio is above
the bound: the residual instability of the fixed scheme (vacuum seeded growth rate 1.3, pV's
linear model 1.43) or the noise of a ratio of maxima, not settled. Unregistered until then.
Runtime about 65 s (box 49 s, one level 16 s).

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
TIME_STEP_NBR = 600
PPC = 100
CHECK_FROM = 0.3
MAX_RELATIVE_GROWTH = 2.0
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


def growths(diag_dir):
    """{level: max |B_z| over the level's patches after CHECK_FROM, over its value at CHECK_FROM}"""
    series = {}
    with h5py.File(f"{diag_dir}/EM_B.h5", "r") as f:
        for t in sorted(f["t"], key=float):
            for level, patches in f["t"][t].items():
                values = [np.abs(p["EM_B_z"][()]).max() for p in patches.values()]
                if values:
                    series.setdefault(int(level[2:]), []).append((float(t), max(values)))
    out = {}
    for level, s in sorted(series.items()):
        times, bz = np.array(s).T
        out[level] = bz[times > CHECK_FROM].max() / bz[np.argmin(np.abs(times - CHECK_FROM))]
    return out


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
        return growths(diag_dir)

    def test_transverse_fields_stay_at_noise_level(self):
        one = self.run_case(refined=False)[0]
        box = self.run_case(refined=True)
        print(f"one level: max |B_z| up to {one:.2f}x its value at t = {CHECK_FROM}")
        for level, growth in box.items():
            print(f"box level {level}: {growth:.2f}x, ratio {growth / one:.2f}")
        self.assertEqual(sorted(box), [0, 1])
        for level, growth in box.items():
            self.assertLess(growth / one, MAX_RELATIVE_GROWTH, f"level {level}")


if __name__ == "__main__":
    unittest.main()
