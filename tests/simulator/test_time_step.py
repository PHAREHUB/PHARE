#
# Config-level validation tests for the time_step option (constant scalar vs adaptive dict).
# These only construct ph.Simulation (pharein) and never run the simulator, so they are cheap
# and need no cpp module / MPI / HighFive.
#

import unittest
import numpy as np
import pyphare.pharein as ph

# minimal valid geometry; time parameters are supplied per-test
baseArgs = dict(
    boundary_types="periodic",
    cells=np.array([20]),
    dl=0.3,
)


# pure pharein config validation: never runs the simulator, so no MPI / cpp module needed
class TimeStepValidation(unittest.TestCase):
    def setUp(self):
        ph.global_vars.sim = None

    def tearDown(self):
        ph.global_vars.sim = None

    # ---- constant (scalar time_step, default) ----------------------------------------------

    def test_constant_is_the_default(self):
        sim = ph.Simulation(time_step=0.001, time_step_nbr=10, **baseArgs)
        self.assertEqual(sim.time_stepper.mode, "constant")
        self.assertEqual(sim.time_step, 0.001)
        self.assertEqual(sim.time_step_nbr, 10)

    def test_constant_final_time_and_step(self):
        sim = ph.Simulation(time_step=0.001, final_time=1.0, **baseArgs)
        self.assertEqual(sim.time_stepper.mode, "constant")
        self.assertEqual(sim.time_step, 0.001)

    def test_constant_final_time_and_nbr_has_no_time_step_kwarg(self):
        sim = ph.Simulation(time_step_nbr=10, final_time=1.0, **baseArgs)
        self.assertEqual(sim.time_stepper.mode, "constant")
        self.assertEqual(sim.time_step_nbr, 10)

    def test_constant_dict_with_value(self):
        sim = ph.Simulation(
            time_step={"mode": "constant", "value": 0.001},
            time_step_nbr=10,
            **baseArgs,
        )
        self.assertEqual(sim.time_stepper.mode, "constant")
        self.assertEqual(sim.time_step, 0.001)
        self.assertEqual(sim.time_step_nbr, 10)

    def test_constant_dict_rejects_unknown_key(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "constant", "dt": 0.001},  # key is "value"
                final_time=1.0,
                **baseArgs,
            )

    # ---- adaptive (dict time_step) ---------------------------------------------------------

    def test_adaptive_accepts_final_time_and_cfl(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl": 0.4},
            final_time=1.0,
            **baseArgs,
        )
        self.assertEqual(sim.time_stepper.mode, "adaptive")
        self.assertEqual(sim.time_stepper.cfl, 0.4)
        self.assertEqual(sim.final_time, 1.0)
        # with adaptive dt these are unknown ahead of the run
        self.assertIsNone(sim.time_step)
        self.assertIsNone(sim.time_step_nbr)

    def test_adaptive_fourier_defaults_to_cfl(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl": 0.4}, final_time=1.0, **baseArgs
        )
        self.assertEqual(sim.time_stepper.fourier, 0.4)

    def test_adaptive_fourier_explicit(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl": 0.4, "fourier": 0.2},
            final_time=1.0,
            **baseArgs,
        )
        self.assertEqual(sim.time_stepper.fourier, 0.2)

    def test_adaptive_requires_cfl(self):
        with self.assertRaises(ValueError):
            ph.Simulation(time_step={"mode": "adaptive"}, final_time=1.0, **baseArgs)

    def test_adaptive_requires_final_time(self):
        with self.assertRaises(ValueError):
            ph.Simulation(time_step={"mode": "adaptive", "cfl": 0.4}, **baseArgs)

    def test_adaptive_rejects_time_step_nbr(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl": 0.4},
                final_time=1.0,
                time_step_nbr=10,
                **baseArgs,
            )

    def test_adaptive_rejects_non_positive_cfl(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl": 0.0},
                final_time=1.0,
                **baseArgs,
            )

    def test_adaptive_rejects_unknown_dict_key(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl": 0.4, "courant": 0.2},
                final_time=1.0,
                **baseArgs,
            )

    def test_adaptive_fourier_rejects_non_positive(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl": 0.4, "fourier": 0.0},
                final_time=1.0,
                **baseArgs,
            )

    def test_adaptive_rejects_final_time_before_restart_time(self):
        with self.assertRaises(RuntimeError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl": 0.4},
                final_time=1.0,
                restart_options={"mode": "overwrite", "restart_time": 2.0},
                **baseArgs,
            )

    def test_adaptive_accepts_final_time_after_restart_time(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl": 0.4},
            final_time=5.0,
            restart_options={"mode": "overwrite", "restart_time": 2.0},
            **baseArgs,
        )
        self.assertEqual(sim.final_time, 5.0)

    # ---- unknown mode ----------------------------------------------------------------------

    def test_unknown_mode_raises(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "variable", "cfl": 0.4},  # only "adaptive" is valid
                final_time=1.0,
                **baseArgs,
            )


if __name__ == "__main__":
    unittest.main()
