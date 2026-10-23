#!/usr/bin/env python3
"""
Two ion species of different q/m drifting relative to each other in a uniform
plasma, with a uniform magnetic field and no current. With J = 0 the total force
on the ions is (sum_s q_s n_s V_s - Ne Ve) x B, which vanishes only if the
electron velocity is the charge-weighted ion velocity. The total ion momentum,
i.e. the mass-weighted bulk velocity, must then stay constant on every level.
Using the mass-weighted velocity in Ohm's law instead makes it rotate (#1466).
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.core.box import Box1D
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

time_step = 0.005
time_step_nbr = 100
ppc = 100

# protons drifting across B through sodium at rest. With the mass-weighted
# velocity in Ohm's law, dVi_z/dt ~ 0.2 here, so Vi_z ~ 0.1 at the end of the run
protons = {"mass": 1.0, "density": 1.0, "vbulky": 1.0}
sodium = {"mass": 23.0, "density": 0.1, "vbulky": 0.0}
vth = 0.1

tolerance = 0.01


def config():
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=50,
        dl=0.2,
        refinement_boxes={"L0": {"B0": Box1D(15, 34)}},
        hyper_resistivity=0.002,
        resistivity=0.001,
        strict=True,
    )

    def constant(value):
        return lambda x: value + 0.0 * x

    def population(species, seed):
        return {
            "charge": 1,
            "mass": species["mass"],
            "density": constant(species["density"]),
            "vbulkx": constant(0.0),
            "vbulky": constant(species["vbulky"]),
            "vbulkz": constant(0.0),
            "vthx": constant(vth),
            "vthy": constant(vth),
            "vthz": constant(vth),
            "nbr_part_per_cell": ppc,
            "init": {"seed": seed},
        }

    ph.MaxwellianFluidModel(
        bx=constant(1.0),
        by=constant(0.0),
        bz=constant(0.0),
        protons=population(protons, seed=1),
        sodium=population(sodium, seed=2),
    )
    ph.ElectronModel(closure="isothermal", Te=0.1)
    return sim


def level_bulk_velocity(simulator, ilvl):
    """mean of the ion bulk velocity over the physical nodes of a level"""
    patch_level = simulator.data_wrangler().getPatchLevel(ilvl)
    mean = []
    for getter in (patch_level.getVix, patch_level.getViy, patch_level.getViz):
        values = [
            np.asarray(patch.data)[
                patch.nGhosts : patch.nGhosts + patch.upper[0] - patch.lower[0] + 2
            ]
            for patch in getter()
        ]
        mean.append(np.concatenate(values).mean())
    return np.array(mean)


class MultiSpeciesMomentumTest(SimulatorTest):
    def tearDown(self):
        if hasattr(self, "simulator"):
            self.simulator.reset()
        super().tearDown()

    def test_ion_momentum_is_conserved(self):
        self.simulator = Simulator(config()).initialize()
        nbr_levels = self.simulator.data_wrangler().getNumberOfLevels()
        self.assertEqual(nbr_levels, 2)

        initial = [level_bulk_velocity(self.simulator, ilvl) for ilvl in range(2)]
        for _ in range(time_step_nbr):
            self.simulator.advance()
        final = [level_bulk_velocity(self.simulator, ilvl) for ilvl in range(2)]

        expected_vy = (
            protons["mass"] * protons["density"] * protons["vbulky"]
            + sodium["mass"] * sodium["density"] * sodium["vbulky"]
        ) / (protons["mass"] * protons["density"] + sodium["mass"] * sodium["density"])

        for ilvl in range(2):
            print(f"L{ilvl} initial {initial[ilvl]} final {final[ilvl]}")
            np.testing.assert_allclose(initial[ilvl][1], expected_vy, atol=tolerance)
            np.testing.assert_allclose(final[ilvl], initial[ilvl], atol=tolerance)


if __name__ == "__main__":
    unittest.main()
