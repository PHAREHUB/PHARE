import unittest
import warnings

import numpy as np

from pyphare.pharein import global_vars, simulation
from pyphare.pharein.fully_kinetic import check_fully_kinetic
from pyphare.pharein.maxwellian_fluid_model import MaxwellianFluidModel

cells = 64
dl = 0.04
L = cells * dl


def uniform(value):
    def fn(x):
        return value + 0.0 * x

    return fn


def make_sim(**overrides):
    kwargs = {
        "time_step": 0.002,
        "time_step_nbr": 10,
        "cells": cells,
        "dl": dl,
        "interp_order": 1,
        "model_options": ["FullyKineticModel"],
        "c_over_va": 10.0,
    }
    kwargs.update(overrides)
    return simulation.Simulation(**kwargs)


def pop(charge, mass, density, T, vx=0.0):
    return {
        "charge": charge,
        "mass": mass,
        "density": uniform(density),
        "vbulkx": uniform(vx),
        "vbulky": uniform(0.0),
        "vbulkz": uniform(0.0),
        "vthx": uniform(np.sqrt(T / mass)),
        "vthy": uniform(np.sqrt(T / mass)),
        "vthz": uniform(np.sqrt(T / mass)),
    }


def make_model(
    electron_density=1.0,
    electron_charge=-1,
    electron_mass=0.04,
    bx=1.0,
    Te=0.04,
    vx=0.0,
):
    return MaxwellianFluidModel(
        bx=uniform(bx),
        by=uniform(0.0),
        bz=uniform(0.0),
        protons=pop(1, 1.0, 1.0, 0.01),
        electrons=pop(electron_charge, electron_mass, electron_density, Te, vx),
    )


class TestFullyKineticChecks(unittest.TestCase):
    def setUp(self):
        global_vars.sim = None

    def tearDown(self):
        global_vars.sim = None

    def test_a_neutral_well_resolved_setup_passes(self):
        sim = make_sim()
        make_model()
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_fully_kinetic(sim)

    def test_c_over_va_is_required(self):
        with self.assertRaisesRegex(ValueError, "requires c_over_va"):
            make_sim(c_over_va=None)

    def test_c_over_va_is_only_for_fully_kinetic(self):
        with self.assertRaisesRegex(ValueError, "only used by FullyKineticModel"):
            make_sim(model_options=["HybridModel"])

    def test_light_cfl_is_enforced(self):
        # c dt / dx = 10 * 0.005 / 0.04 > 1
        with self.assertRaisesRegex(ValueError, "light CFL"):
            make_sim(time_step=0.005)

    def test_adaptive_time_step_is_refused(self):
        with self.assertRaisesRegex(ValueError, "constant time step"):
            make_sim(
                time_step={"mode": "adaptive", "cfl_wave": 0.5, "cfl_diffusive": 0.5},
                time_step_nbr=None,
                final_time=0.02,
            )

    def test_refinement_is_refused(self):
        with self.assertRaisesRegex(ValueError, "one level only"):
            make_sim(refinement="tagging", max_nbr_levels=2)
        global_vars.sim = None
        with self.assertRaisesRegex(ValueError, "one level only"):
            make_sim(refinement_boxes={"L0": {"B0": [(10,), (20,)]}})

    def test_restarts_are_refused(self):
        with self.assertRaisesRegex(ValueError, "restarts"):
            make_sim(
                restart_options={
                    "dir": "restarts",
                    "mode": "overwrite",
                    "timestamps": [0.01],
                }
            )

    def test_net_charge_is_refused(self):
        sim = make_sim()
        make_model(electron_density=0.9)
        with self.assertRaisesRegex(ValueError, "net charge"):
            check_fully_kinetic(sim)

    def test_light_and_plasma_limits_combine(self):
        # c dt / dx = 0.9 and w_pe dt = 1: each fine alone, unstable together
        sim = make_sim(time_step=0.0036)
        make_model(electron_mass=(10 * 0.0036) ** 2)
        with self.assertRaisesRegex(ValueError, "unstable"):
            check_fully_kinetic(sim)

    def test_the_plasma_frequency_is_the_total_one(self):
        # two beams of half density: w_pb dt = 1.5 each, sqrt(2) times more together,
        # (2.12 / 2)^2 + (c dt / dx)^2 = 1.12 + 0.25 > 1
        sim = make_sim()
        mass = (10 * np.sqrt(0.5) * 0.002 / 1.5) ** 2
        MaxwellianFluidModel(
            bx=uniform(1.0),
            by=uniform(0.0),
            bz=uniform(0.0),
            protons=pop(1, 1.0, 1.0, 0.01),
            beam1=pop(-1, mass, 0.5, 1e-6),
            beam2=pop(-1, mass, 0.5, 1e-6),
        )
        with self.assertRaisesRegex(ValueError, "unstable"):
            check_fully_kinetic(sim)

    def test_bulk_velocity_beyond_a_cell_per_step_is_refused(self):
        # 25 * 0.002 = 0.05 > dx
        sim = make_sim()
        make_model(vx=25.0)
        with self.assertRaisesRegex(ValueError, "bulk velocity"):
            check_fully_kinetic(sim)

    def test_tail_warning_counts_the_particles(self):
        # 6400 electrons: sqrt(2 ln 6400) + 1 = 5.2 thermal velocities,
        # v_th = sqrt(0.6 / 0.04) = 3.9: 5.2 * 3.9 * 0.002 = 0.040 >= dx
        sim = make_sim()
        make_model(Te=0.6)
        with self.assertWarns(UserWarning) as caught:
            check_fully_kinetic(sim)
        self.assertIn("thermal", " ".join(str(w.message) for w in caught.warnings))

    def test_fast_particles_compared_to_light_warn(self):
        # v_th = sqrt(0.16 / 0.04) = 2, 3 v_th / c = 0.6
        sim = make_sim()
        make_model(Te=0.16)
        with self.assertWarns(UserWarning) as caught:
            check_fully_kinetic(sim)
        self.assertIn(
            "speed of light", " ".join(str(w.message) for w in caught.warnings)
        )

    def test_strict_turns_warnings_into_errors(self):
        sim = make_sim(strict=True)
        make_model(Te=0.16)
        with self.assertRaisesRegex(RuntimeError, "speed of light"):
            check_fully_kinetic(sim)

    def test_unstable_plasma_period_is_refused(self):
        # omega_pe = c sqrt(n / me) = 10 / sqrt(1e-6) = 1e4, omega_pe dt = 20
        sim = make_sim()
        make_model(electron_mass=1e-6)
        with self.assertRaisesRegex(ValueError, "unstable"):
            check_fully_kinetic(sim)

    def test_population_charge_and_mass_are_checked(self):
        sim = make_sim()
        make_model(electron_charge=0)
        with self.assertRaisesRegex(ValueError, "non zero charge"):
            check_fully_kinetic(sim)

    def test_poorly_resolved_setups_warn(self):
        # me = 0.01: omega_pe dt = 100 * 0.0035 = 0.35, Omega_ce dt = 20 / 0.01 * 0.0035 = 7,
        # lambda_De = sqrt(1e-5 / 0.01) / 100 = 3e-4 << dx
        sim = make_sim(time_step=0.0035)
        make_model(electron_mass=0.01, bx=20.0, Te=1e-5)
        with self.assertWarns(UserWarning) as caught:
            check_fully_kinetic(sim)
        messages = " ".join(str(w.message) for w in caught.warnings)
        for word in ["plasma frequency", "cyclotron", "Debye"]:
            self.assertIn(word, messages)

    def test_fast_thermal_electrons_warn(self):
        # 4 v_th,e dt = 4 * sqrt(4 / 0.04) * 0.002 = 0.08 > dx
        sim = make_sim()
        make_model(Te=4.0)
        with self.assertWarns(UserWarning) as caught:
            check_fully_kinetic(sim)
        self.assertIn("thermal", " ".join(str(w.message) for w in caught.warnings))


if __name__ == "__main__":
    unittest.main()
