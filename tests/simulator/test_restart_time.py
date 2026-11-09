#!/usr/bin/env python3

import unittest
from pathlib import Path

import numpy as np

import pyphare.pharein as ph
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

    def run_recording(
        self,
        restart_options,
        diag_dir="diags",
        diags=False,
        advance_past_end=False,
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

        simulator.exit_on_exception = False
        try:
            self.simulator = Simulator(ph.global_vars.sim).initialize()
            times = [self.simulator.currentTime()]
            end_time = self.simulator.cpp_sim.endTime()
            while self.simulator.currentTime() < end_time:
                self.simulator.advance()
                times.append(self.simulator.currentTime())
            if advance_past_end:
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


if __name__ == "__main__":
    startMPI()
    unittest.main()
