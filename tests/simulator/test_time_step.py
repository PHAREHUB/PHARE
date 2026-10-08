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
        self.assertEqual(sim.time_stepper.time_step, 0.001)
        self.assertEqual(sim.time_stepper.time_step_nbr, 10)

    def test_constant_final_time_and_step(self):
        sim = ph.Simulation(time_step=0.001, final_time=1.0, **baseArgs)
        self.assertEqual(sim.time_stepper.mode, "constant")
        self.assertEqual(sim.time_stepper.time_step, 0.001)

    def test_constant_final_time_and_nbr_has_no_time_step_kwarg(self):
        sim = ph.Simulation(time_step_nbr=10, final_time=1.0, **baseArgs)
        self.assertEqual(sim.time_stepper.mode, "constant")
        self.assertEqual(sim.time_stepper.time_step_nbr, 10)

    def test_constant_dict_with_value(self):
        sim = ph.Simulation(
            time_step={"mode": "constant", "value": 0.001},
            time_step_nbr=10,
            **baseArgs,
        )
        self.assertEqual(sim.time_stepper.mode, "constant")
        self.assertEqual(sim.time_stepper.time_step, 0.001)
        self.assertEqual(sim.time_stepper.time_step_nbr, 10)

    def test_constant_dict_rejects_unknown_key(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "constant", "dt": 0.001},  # key is "value"
                final_time=1.0,
                **baseArgs,
            )

    # ---- adaptive (dict time_step) ---------------------------------------------------------

    def test_adaptive_accepts_final_time_and_cfl_wave(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl_wave": 0.4},
            final_time=1.0,
            **baseArgs,
        )
        self.assertEqual(sim.time_stepper.mode, "adaptive")
        self.assertEqual(sim.time_stepper.cfl_wave, 0.4)
        self.assertEqual(sim.time_stepper.final_time, 1.0)
        # with adaptive dt these are unknown ahead of the run
        self.assertFalse(hasattr(sim.time_stepper, "time_step"))
        self.assertFalse(hasattr(sim.time_stepper, "time_step_nbr"))

    def test_adaptive_cfl_diffusive_defaults_to_cfl_wave(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl_wave": 0.4}, final_time=1.0, **baseArgs
        )
        self.assertEqual(sim.time_stepper.cfl_diffusive, 0.4)

    def test_adaptive_cfl_diffusive_explicit(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl_wave": 0.4, "cfl_diffusive": 0.2},
            final_time=1.0,
            **baseArgs,
        )
        self.assertEqual(sim.time_stepper.cfl_diffusive, 0.2)

    def test_adaptive_requires_cfl_wave(self):
        with self.assertRaises(ValueError):
            ph.Simulation(time_step={"mode": "adaptive"}, final_time=1.0, **baseArgs)

    def test_adaptive_requires_final_time(self):
        with self.assertRaises(ValueError):
            ph.Simulation(time_step={"mode": "adaptive", "cfl_wave": 0.4}, **baseArgs)

    def test_adaptive_rejects_time_step_nbr(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl_wave": 0.4},
                final_time=1.0,
                time_step_nbr=10,
                **baseArgs,
            )

    def test_adaptive_rejects_non_positive_cfl_wave(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl_wave": 0.0},
                final_time=1.0,
                **baseArgs,
            )

    def test_adaptive_rejects_unknown_dict_key(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl_wave": 0.4, "courant": 0.2},
                final_time=1.0,
                **baseArgs,
            )

    def test_adaptive_cfl_diffusive_rejects_non_positive(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl_wave": 0.4, "cfl_diffusive": 0.0},
                final_time=1.0,
                **baseArgs,
            )

    def test_adaptive_rejects_final_time_before_restart_time(self):
        with self.assertRaises(RuntimeError):
            ph.Simulation(
                time_step={"mode": "adaptive", "cfl_wave": 0.4},
                final_time=1.0,
                restart_options={"mode": "overwrite", "restart_time": 2.0},
                **baseArgs,
            )

    def test_adaptive_accepts_final_time_after_restart_time(self):
        sim = ph.Simulation(
            time_step={"mode": "adaptive", "cfl_wave": 0.4},
            final_time=5.0,
            restart_options={"mode": "overwrite", "restart_time": 2.0},
            **baseArgs,
        )
        self.assertEqual(sim.time_stepper.final_time, 5.0)

    # ---- unknown mode ----------------------------------------------------------------------

    def test_unknown_mode_raises(self):
        with self.assertRaises(ValueError):
            ph.Simulation(
                time_step={"mode": "variable", "cfl_wave": 0.4},  # only "adaptive" is valid
                final_time=1.0,
                **baseArgs,
            )

    # ---- zero-duration run ------------------------------------------------------------------

    def test_zero_duration_accepts_initial_timestamp(self):
        sim = ph.Simulation(time_step=0.005, final_time=0.0, **baseArgs)
        self.assertEqual(sim.time_stepper.time_step_nbr, 0)
        self.assertEqual(sim.time_stepper.time_step, 1.0)
        diag = ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([0.0]))
        np.testing.assert_array_equal(diag.write_timestamps, [0.0])

    def test_zero_duration_from_time_step_nbr(self):
        sim = ph.Simulation(time_step_nbr=0, final_time=0.0, **baseArgs)
        self.assertEqual(sim.time_stepper.time_step_nbr, 0)
        self.assertEqual(sim.time_stepper.time_step, 1.0)

    def test_zero_duration_after_restart_accepts_restart_timestamp(self):
        ph.Simulation(
            time_step=0.005,
            final_time=2.5,
            restart_options={"mode": "overwrite", "restart_time": 2.5},
            **baseArgs,
        )
        diag = ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([2.5]))
        np.testing.assert_array_equal(diag.write_timestamps, [2.5])

    def test_zero_duration_rejects_timestamp_after_final_time(self):
        ph.Simulation(time_step=0.005, final_time=0.0, **baseArgs)
        with self.assertRaises(RuntimeError):
            ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([0.005]))

    def test_off_grid_timestamp_still_rejected(self):
        ph.Simulation(time_step=0.005, final_time=1.0, **baseArgs)
        with self.assertRaises(RuntimeError):
            ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([0.0012]))


if __name__ == "__main__":
    unittest.main()
