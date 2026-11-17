#!/usr/bin/env python3
"""
FullyKineticModel, one static refined box, vacuum: a seeded B_z does not grow (P4 T9 residual,
docs/superpowers/specs/p4/t9_residual.md).

1D vacuum (n = 1e-6, 4 particles per cell and species, as T8), 32 coarse cells, box over root
cells 8 to 23, order 1, the refined runs' options (one filter pass, Marder s = 0.5), c dt / dx =
0.5. B_z is seeded at 1e-6 on the coarse modes 8 to 16 with random phases (fixed seed). On one
level it stays bounded; with the box, the fine B ghosts of the second substep interpolated in time
towards a coarse B the sync then replaced, and the interface grew at gamma 1.33 per unit time.
Deterministic (no plasma: the particles carry no current that matters).
Checks: on each level, max |B_z| over the physical nodes at t = T_END over its value at
t = T_START is at most MAX_GROWTH.

pV's linear model (5 phase sets) and PHARE (this phase set), growth from 0.2 to 8.0:

    one level            model 0.87 - 1.31    PHARE 0.87 (L0)
    HEAD (before B-x)    model 1.9e3 - 6.7e3  PHARE 5.67e3 (L0), 1.36e4 (L1)
    B-x partial reflux   model 0.28 - 0.95 (L0), 0.14 - 0.27 (L1)

The bound 2 catches gamma >= ln 2 / 7.8 = 0.09 with a 1.5x margin over the beating of the 9 modes.
"""

import unittest

import h5py
import numpy as np

import pyphare.pharein as ph
from pyphare.simulator import simulator
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells, dl, time_step = 32, 0.04, 0.002
c_over_va, mass_e = 10.0, 0.04
L = cells * dl
BOX = (8, 23)
MODES = np.arange(8, 17)
PHASE_SEED = 7
T_START, T_END = 0.2, 8.0
MAX_GROWTH = 2.0
diag_dir = "phare_outputs/fully_kinetic_amr_vacuum_stability"


def uniform(value):
    return lambda x: value + 0.0 * x


def seeded_bz():
    phases = np.random.default_rng(PHASE_SEED).uniform(0, 2 * np.pi, MODES.size)
    return lambda x: 1e-6 * sum(
        np.cos(2 * np.pi * m * x / L + phase) for m, phase in zip(MODES, phases)
    )


def vacuum(charge, mass, T, seed):
    vth = uniform(np.sqrt(T / mass))
    return {
        "charge": charge,
        "mass": mass,
        "density": uniform(1e-6),
        "vbulkx": uniform(0.0),
        "vbulky": uniform(0.0),
        "vbulkz": uniform(0.0),
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
        "nbr_part_per_cell": 4,
        "init": {"seed": seed},
    }


def config():
    sim = ph.Simulation(
        time_step=time_step,
        final_time=T_END,
        cells=[cells],
        dl=[dl],
        interp_order=1,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=1,
        divergence_cleaning=0.5,
        refinement_boxes={"L0": {"B0": [(BOX[0],), (BOX[1],)]}},
        strict=False,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite", "fine_dump_lvl_max": 10},
        },
    )
    ph.MaxwellianFluidModel(
        bx=uniform(0.0),
        by=uniform(0.0),
        bz=seeded_bz(),
        protons=vacuum(1, 1.0, 0.01, 1),
        electrons=vacuum(-1, mass_e, 0.04, 2),
    )
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=[T_START, T_END])
    return sim


def max_bz(f, t):
    """{level: max |B_z| over the physical nodes of the level's patches} at the dump nearest t"""
    key = min(f["t"], key=lambda k: abs(float(k) - t))
    assert abs(float(key) - t) < 0.5 * time_step, f"no dump at t = {t}"
    out = {}
    for level, patches in f["t"][key].items():
        for patch in patches.values():
            ds = patch["EM_B_z"]
            g = int(ds.attrs["ghosts"])
            value = np.abs(ds[()][g:-g]).max()
            ilvl = int(level[2:])
            out[ilvl] = max(out.get(ilvl, 0.0), value)
    return out


class FullyKineticAMRVacuumStabilityTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None
        simulator.exit_on_exception = True

    def test_seeded_vacuum_does_not_grow(self):
        simulator.exit_on_exception = False  # a run stopped by the growth fails below
        Simulator(config()).run().reset()
        with h5py.File(f"{diag_dir}/EM_B.h5", "r") as f:
            start, end = max_bz(f, T_START), max_bz(f, T_END)
        self.assertEqual(sorted(start), [0, 1])
        growth = {level: end[level] / start[level] for level in sorted(start)}
        for level, g in growth.items():
            print(f"level {level}: max |B_z| x{g:.3g} from t = {T_START} to {T_END}")
        for level, g in growth.items():
            self.assertLessEqual(g, MAX_GROWTH, f"level {level}")


if __name__ == "__main__":
    unittest.main()
