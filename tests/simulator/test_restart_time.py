#!/usr/bin/env python3

import subprocess
import sys
import unittest
from pathlib import Path

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator import simulator
from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator import SimulatorTest

GAMMA = 5.0 / 3.0
V = (0.3, -0.2, 0.1)
B = (0.5, 0.4, 0.3)
P = 0.2

OUT = "phare_outputs/restart_time"
CONSTANT = dict(time_step=0.1, time_step_nbr=6)
ADAPTIVE = dict(time_step={"mode": "adaptive", "cfl_wave": 0.5}, final_time=1.0)


def uniform_model():
    def const(value):
        return lambda x, y: value + 0.0 * x

    ph.MHDModel(
        density=const(1.0),
        vx=const(V[0]),
        vy=const(V[1]),
        vz=const(V[2]),
        bx=const(B[0]),
        by=const(B[1]),
        bz=const(B[2]),
        p=const(P),
    )


class RecordingTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None
        super().tearDown()

    def out(self):
        out = f"{OUT}/{self._testMethodName}"
        if out not in self.diag_dirs:
            self.register_diag_dir_for_cleanup(out)
        return out

    def build(
        self,
        restart_options,
        diag_dir="diags",
        diags=False,
        diag_timestamps=None,
        **time_kwargs,
    ):
        out = self.out()
        ph.global_vars.sim = None
        ph.Simulation(
            interp_order=None,
            model_options=["MHDModel"],
            reconstruction="WENOZ",
            limiter="None",
            riemann="Rusanov",
            mhd_timestepper="TVDRK3",
            hall=False,
            eta=0.0,
            nu=0.0,
            gamma=GAMMA,
            max_mhd_level=1,
            cells=(20, 10),
            dl=(0.3, 0.2),
            restart_options=dict(
                dir=f"{out}/restarts", mode="overwrite", **restart_options
            ),
            diag_options={
                "format": "phareh5",
                "options": {"dir": f"{out}/{diag_dir}", "mode": "overwrite"},
            },
            **time_kwargs,
        )
        uniform_model()
        if diags:
            from pyphare.pharein.diagnostics import all_timestamps

            ph.MHDDiagnostics(
                quantity="rho", write_timestamps=all_timestamps(ph.global_vars.sim)
            )
        if diag_timestamps is not None:
            ph.MHDDiagnostics(quantity="rho", write_timestamps=diag_timestamps)
        return ph.global_vars.sim

    def run_recording(self, restart_options, advance_past_end=False, **kwargs):
        self.build(restart_options, **kwargs)
        simulator.exit_on_exception = False
        try:
            self.simulator = Simulator(ph.global_vars.sim).initialize()
            times = [self.simulator.currentTime()]
            end_time = self.simulator.cpp_sim.endTime()
            while self.simulator.currentTime() < end_time:
                self.simulator.advance()
                times.append(self.simulator.currentTime())
            if advance_past_end:
                self.assertEqual(self.simulator.currentTime(), end_time)
                with self.assertRaisesRegex(RuntimeError, "already at its final time"):
                    self.simulator.advance()
        finally:
            simulator.exit_on_exception = True
            if self.simulator is not None:
                self.simulator.reset()
            self.simulator = None
        return times

    def diag_times(self, diag_dir):
        import h5py

        path = next(Path(f"{self.out()}/{diag_dir}").glob("*.h5"))
        with h5py.File(path, "r") as f:
            return sorted(f["t"].keys(), key=float)


class RestartTimeTest(RecordingTest):
    def assert_ulp_close(self, times, reference):
        times, reference = np.array(times), np.array(reference)
        self.assertEqual(len(times), len(reference))
        np.testing.assert_array_less(
            np.abs(times - reference),
            4 * np.spacing(reference) + np.finfo(float).tiny,
        )

    def test_constant_restart_resumes_at_restored_time(self):
        tA = self.run_recording(dict(timestamps=[0.3]), **CONSTANT)
        tB = self.run_recording(
            dict(restart_time=0.300004), time_step=0.1, time_step_nbr=3
        )
        self.assertEqual(tB[0], tA[3])
        self.assertEqual(len(tB), 4)
        self.assert_ulp_close(tB, tA[3:])

    def test_constant_restart_with_new_dt(self):
        tA = self.run_recording(dict(timestamps=[0.3]), **CONSTANT)
        tB = self.run_recording(dict(restart_time=0.3), time_step=0.05, time_step_nbr=4)
        self.assertEqual(tB[0], tA[3])
        self.assertEqual(len(tB), 5)
        self.assertLessEqual(abs(tB[-1] - (tA[3] + 0.05 * 4)), np.spacing(tB[-1]))

    def test_adaptive_restart_resumes_at_restored_time(self):
        tA = self.run_recording(dict(timestamps=[0.5]), **ADAPTIVE)
        tB = self.run_recording(dict(restart_time="auto"), **ADAPTIVE)
        self.assertIn(tB[0], tA)
        self.assertEqual(tB[-1], 1.0)
        self.assert_ulp_close(tB, tA[tA.index(tB[0]) :])

    def test_adaptive_to_constant_restart(self):
        tA = self.run_recording(dict(timestamps=[0.5]), **ADAPTIVE)
        tB = self.run_recording(
            dict(restart_time="auto"), time_step=0.01, time_step_nbr=5
        )
        self.assertIn(tB[0], tA)
        self.assertEqual(len(tB), 6)
        self.assertLessEqual(abs(tB[-1] - (tB[0] + 0.01 * 5)), np.spacing(tB[-1]))

    def test_adaptive_restart_lands_on_final_time(self):
        restart_time = 63.0 / 64.0
        self.run_recording(
            dict(timestamps=[restart_time]),
            time_step=restart_time / 32,
            time_step_nbr=32,
        )
        for cfl in (0.507, 0.173, 0.365):
            with self.subTest(cfl=cfl):
                tB = self.run_recording(
                    dict(restart_time=restart_time),
                    time_step={"mode": "adaptive", "cfl_wave": cfl},
                    final_time=2.0,
                )
                steps = np.diff(tB)
                self.assertEqual(tB[0], restart_time)
                self.assertEqual(tB[-1], 2.0)
                self.assertGreater(steps.min(), 1e-6 * steps.max())

    def test_restart_with_zero_steps(self):
        tA = self.run_recording(dict(timestamps=[0.3]), **CONSTANT)
        tB = self.run_recording(dict(restart_time=0.3), time_step=0.1, time_step_nbr=0)
        self.assertEqual(tB, [tA[3]])

    def test_typed_restart_time_dumps_at_the_restored_time(self):
        tA = self.run_recording(
            dict(timestamps=[0.009]), time_step=0.001, time_step_nbr=9
        )
        self.assertNotEqual(tA[9], 0.009)
        tB = self.run_recording(
            dict(restart_time=0.009),
            diag_dir="diag_B",
            diag_timestamps=[0.009],
            time_step=0.001,
            final_time=0.009,
        )
        self.assertEqual(tB, [tA[9]])
        from pyphare import cpp

        if cpp.mpi_rank() == 0:
            self.assertEqual(len(self.diag_times("diag_B")), 1)

    def test_typed_restart_timestamps_at_the_restored_time(self):
        tA = self.run_recording(
            dict(timestamps=[0.009]), time_step=0.001, time_step_nbr=9
        )
        tB = self.run_recording(
            dict(restart_time=0.009, timestamps=[0.009, 0.012]),
            time_step=0.001,
            time_step_nbr=3,
        )
        self.assertEqual(tB[0], tA[9])
        self.assertEqual(len(tB), 4)

    def test_constant_restart_dumps_on_the_uninterrupted_steps(self):
        dt = 1.0 / 30.0
        self.run_recording(
            dict(timestamps=[10 * dt]),
            diag_dir="diag_A",
            diags=True,
            time_step=dt,
            time_step_nbr=30,
        )
        self.run_recording(
            dict(restart_time="auto"),
            diag_dir="diag_B",
            diags=True,
            time_step=dt,
            time_step_nbr=20,
        )
        from pyphare import cpp

        if cpp.mpi_rank() == 0:
            self.assertEqual(self.diag_times("diag_B"), self.diag_times("diag_A")[10:])

    def test_restart_timestamps_tolerate_float64_rounding_of_tiny_steps(self):
        tA = self.run_recording(dict(timestamps=[0.3]), **CONSTANT)
        sim = self.build(
            dict(restart_time=0.3, timestamps=[0.30001]),
            time_step=1e-9,
            time_step_nbr=10000,
        )
        self.assertEqual(sim.start_time(), tA[3])
        self.assertEqual(
            list(sim.restart_options["timestamps"]), [sim.time_stepper.final_time]
        )

    def test_constant_restart_off_the_grid_of_the_new_dt(self):
        dt, nbr = 0.07, 3
        tA = self.run_recording(dict(timestamps=[0.3]), **CONSTANT)
        start = tA[3]
        with self.assertRaisesRegex(
            RuntimeError, "inconsistent with simulation.time_step"
        ):
            self.build(
                dict(restart_time=0.3),
                diag_timestamps=[start + dt / 2],
                time_step=dt,
                time_step_nbr=nbr,
            )
        with self.assertRaisesRegex(
            RuntimeError, "inconsistent with simulation.time_step"
        ):
            self.build(
                dict(restart_time=0.3, timestamps=[start + 1.5 * dt]),
                time_step=dt,
                time_step_nbr=nbr,
            )
        tB = self.run_recording(
            dict(restart_time=0.3, timestamps=[start + nbr * dt]),
            diag_dir="diag_B",
            diags=True,
            time_step=dt,
            time_step_nbr=nbr,
        )
        self.assertEqual(tB[0], start)
        self.assertEqual(len(tB), nbr + 1)
        from pyphare import cpp

        if cpp.mpi_rank() == 0:
            self.assertEqual(len(self.diag_times("diag_B")), nbr + 1)

    def test_adaptive_zero_step_restart_dumps_at_the_restored_time(self):
        tA = self.run_recording(
            dict(timestamps=[0.009]), time_step=0.001, time_step_nbr=9
        )
        tB = self.run_recording(
            dict(restart_time=0.009),
            diag_dir="diag_B",
            diag_timestamps=[0.009],
            time_step={"mode": "adaptive", "cfl_wave": 0.5},
            final_time=0.009,
        )
        self.assertEqual(tB, [tA[9]])
        if cpp.mpi_rank() == 0:
            self.assertEqual(len(self.diag_times("diag_B")), 1)

    def test_constant_final_time_and_step_count_lands_on_final_time(self):
        t = self.run_recording({}, final_time=0.9, time_step_nbr=7)
        self.assertEqual(len(t), 8)
        self.assertEqual(t[-1], 0.9)

    def test_adaptive_sliver_is_merged_into_the_last_step(self):
        tR = self.run_recording({}, **ADAPTIVE)
        final = tR[5] + 0.5e-6 * (tR[5] - tR[4])
        t = self.run_recording(
            {}, time_step={"mode": "adaptive", "cfl_wave": 0.5}, final_time=final
        )
        self.assertEqual(len(t), 6)
        self.assertEqual(t[-1], final)

    def strip_phare_time(self):
        import h5py

        for path in Path(f"{self.out()}/restarts").glob("*/restore.*/nodes.*/proc.*"):
            with h5py.File(path, "a") as f:
                del f["phare"].attrs["time"]

    def test_legacy_restart_resumes_at_the_exact_time(self):
        dt = 1.0 / 30.0
        tA = self.run_recording(
            dict(timestamps=[10 * dt]), time_step=dt, time_step_nbr=30
        )
        if cpp.mpi_rank() == 0:
            self.strip_phare_time()
        cpp.mpi_barrier()
        tB = self.run_recording(
            dict(restart_time="auto"), time_step=dt, time_step_nbr=20
        )
        self.assertEqual(tB[0], tA[10])
        self.assertEqual(len(tB), 21)
        self.assert_ulp_close(tB, tA[10:])

    def test_restored_time_without_mpi_layout_match(self):
        tA = self.run_recording(dict(timestamps=[0.3]), **CONSTANT)
        folder = cpp.cpp_etc_lib().restart_path_for_time(f"{self.out()}/restarts", 0.3)
        self.assertEqual(cpp.cpp_etc_lib().restored_time(folder), tA[3])
        if cpp.mpi_rank() == 0:
            code = (
                "from pyphare import cpp; "
                f"print(repr(cpp.cpp_etc_lib().restored_time({folder!r})))"
            )
            out = subprocess.run(
                [sys.executable, "-c", code], capture_output=True, text=True, check=True
            )
            self.assertEqual(float(out.stdout.strip().splitlines()[-1]), tA[3])
        with self.assertRaisesRegex(RuntimeError, "no restart file"):
            cpp.cpp_etc_lib().restored_time(f"{self.out()}/nowhere")


if __name__ == "__main__":
    startMPI()
    unittest.main()
