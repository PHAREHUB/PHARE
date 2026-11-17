import unittest
import warnings
from unittest.mock import patch

import numpy as np

from pyphare.pharein import global_vars, simulation
from pyphare.pharein.fully_kinetic import check_fully_kinetic
from pyphare.pharein.initialize import fully_kinetic as initialize_fully_kinetic
from pyphare.pharein.initialize.general import populate_load_balancer
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


def uniform_nd(value):
    def fn(*xyz):
        return value + 0.0 * xyz[0]

    return fn


def anisotropic_pop(charge, mass, Tx, Ty):
    return {
        "charge": charge,
        "mass": mass,
        "density": uniform_nd(1.0),
        "vbulkx": uniform_nd(0.0),
        "vbulky": uniform_nd(0.0),
        "vbulkz": uniform_nd(0.0),
        "vthx": uniform_nd(np.sqrt(Tx / mass)),
        "vthy": uniform_nd(np.sqrt(Ty / mass)),
        "vthz": uniform_nd(np.sqrt(Tx / mass)),
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

    def test_tagging_is_accepted_with_the_regrid_warning(self):
        # P5-1: the tagger is there, the regrid of the fine levels it creates is not (P5a)
        with self.assertWarnsRegex(UserWarning, "regrid"):
            sim = make_sim(
                refinement="tagging", max_nbr_levels=2, divergence_cleaning=0.5
            )
        self.assertEqual(sim.refinement, "tagging")

    def test_tagging_writes_fk_method_and_leaves_load_balancing_off(self):
        # the C++ tagger factory reads fk_method; tagging turns load balancing on for the
        # hybrid, which the fully kinetic model does not support yet (P5b)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sim = make_sim(
                refinement="tagging", max_nbr_levels=2, divergence_cleaning=0.5
            )
        make_model()
        written = {}

        class Recorder:
            def __getattr__(self, adder):
                return lambda path, value, *args: written.__setitem__(path, value)

        populate_load_balancer(Recorder(), sim)
        with patch.object(initialize_fully_kinetic, "DictPopulator", Recorder):
            initialize_fully_kinetic.populateDict(sim)
        self.assertEqual(
            written["simulation/AMR/refinement/tagging/fk_method"], "default"
        )
        self.assertFalse(written["simulation/AMR/loadbalancing/active"])

    def test_refinement_boxes_are_accepted(self):
        sim = make_sim(
            refinement_boxes={"L0": {"B0": [(10,), (20,)]}}, divergence_cleaning=0.5
        )
        self.assertEqual(sim.max_nbr_levels, 2)

    def test_divergence_cleaning_must_be_set_with_refinement(self):
        # Gauss's law breaks at the level boundary (band swap, erased particles, splitting):
        # the user chooses the cleaning, 0 included, rather than getting none silently
        with self.assertRaisesRegex(
            ValueError, "divergence_cleaning must be set with refinement"
        ):
            make_sim(refinement_boxes={"L0": {"B0": [(10,), (20,)]}})
        for s in [0, 0.5]:
            global_vars.sim = None
            sim = make_sim(
                refinement_boxes={"L0": {"B0": [(10,), (20,)]}}, divergence_cleaning=s
            )
            self.assertEqual(sim.divergence_cleaning, float(s))

    def test_current_filter_is_off_by_default(self):
        self.assertEqual(make_sim().current_filter_passes, 0)

    def test_current_filter_passes_are_kept(self):
        self.assertEqual(make_sim(current_filter_passes=2).current_filter_passes, 2)

    def test_current_filter_passes_must_be_a_non_negative_integer(self):
        for passes in [-1, 1.5, True]:
            global_vars.sim = None
            with self.assertRaisesRegex(ValueError, "current_filter_passes"):
                make_sim(current_filter_passes=passes)

    def test_current_filter_is_only_for_fully_kinetic(self):
        with self.assertRaisesRegex(ValueError, "only used by FullyKineticModel"):
            make_sim(
                model_options=["HybridModel"], c_over_va=None, current_filter_passes=1
            )

    def test_divergence_cleaning_is_off_by_default(self):
        self.assertEqual(make_sim().divergence_cleaning, 0.0)

    def test_divergence_cleaning_is_kept(self):
        self.assertEqual(make_sim(divergence_cleaning=0.5).divergence_cleaning, 0.5)
        global_vars.sim = None
        self.assertEqual(make_sim(divergence_cleaning=1).divergence_cleaning, 1.0)

    def test_divergence_cleaning_must_be_in_the_stable_range(self):
        # s = d dt sum_i 2 / dx_i^2 > 1 grows the grid scale mode of div E - c^2 rho
        for s in [-0.1, 1.01, True, "0.5"]:
            global_vars.sim = None
            with self.assertRaisesRegex(ValueError, "divergence_cleaning"):
                make_sim(divergence_cleaning=s)

    def test_divergence_cleaning_is_only_for_fully_kinetic(self):
        with self.assertRaisesRegex(ValueError, "only used by FullyKineticModel"):
            make_sim(
                model_options=["HybridModel"], c_over_va=None, divergence_cleaning=0.5
            )

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
        with self.assertRaisesRegex(ValueError, "frozen background"):
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

    def test_local_charge_separation_is_refused(self):
        # globally neutral, but the ions are modulated and the electrons are not
        sim = make_sim()
        protons = pop(1, 1.0, 1.0, 0.01)
        protons["density"] = lambda x: 1.0 + 0.1 * np.cos(2 * np.pi * x / L)
        MaxwellianFluidModel(
            bx=uniform(1.0),
            by=uniform(0.0),
            bz=uniform(0.0),
            protons=protons,
            electrons=pop(-1, 0.04, 1.0, 0.04),
        )
        with self.assertRaisesRegex(ValueError, "frozen background"):
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

    def test_debye_length_uses_the_coldest_simulated_direction(self):
        # 2D, v_th,x = sqrt(0.02 / 0.04) = 0.71: dx / lambda_D,x = 2.8 < 3,
        # v_th,y = sqrt(1e-5 / 0.04) = 0.016: dy / lambda_D,y = 127
        sim = make_sim(cells=(32, 32), dl=(dl, dl))
        electrons = anisotropic_pop(-1, 0.04, Tx=0.02, Ty=1e-5)
        MaxwellianFluidModel(
            bx=uniform_nd(1.0),
            by=uniform_nd(0.0),
            bz=uniform_nd(0.0),
            protons=anisotropic_pop(1, 1.0, Tx=0.01, Ty=0.01),
            electrons=electrons,
        )
        with self.assertWarns(UserWarning) as caught:
            check_fully_kinetic(sim)
        self.assertIn("Debye", " ".join(str(w.message) for w in caught.warnings))

    def test_debye_length_ignores_directions_not_simulated(self):
        # 1D: a cold y direction is not resolved by the grid, no heating warning
        sim = make_sim()
        MaxwellianFluidModel(
            bx=uniform(1.0),
            by=uniform(0.0),
            bz=uniform(0.0),
            protons=anisotropic_pop(1, 1.0, Tx=0.01, Ty=0.01),
            electrons=anisotropic_pop(-1, 0.04, Tx=0.02, Ty=1e-5),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            check_fully_kinetic(sim)

    def test_fast_thermal_electrons_warn(self):
        # 4 v_th,e dt = 4 * sqrt(4 / 0.04) * 0.002 = 0.08 > dx
        sim = make_sim()
        make_model(Te=4.0)
        with self.assertWarns(UserWarning) as caught:
            check_fully_kinetic(sim)
        self.assertIn("thermal", " ".join(str(w.message) for w in caught.warnings))


if __name__ == "__main__":
    unittest.main()
