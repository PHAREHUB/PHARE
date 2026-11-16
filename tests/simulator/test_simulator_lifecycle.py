#!/usr/bin/env python3
#
# Simulator.run() keeps the simulator alive: its data wrangler works afterwards, the
#  diagnostic files are complete without a reset, and setting up another simulator
#  releases it. Setting up while one is live but unfinished raises.

import unittest

from ddt import data, ddt

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest, populate_simulation

time_step_nbr = 2
final_time = 0.002
timestamps = [0, final_time]
diag_dir = "phare_outputs/simulator_lifecycle"


def config(diag_format="phareh5", **diag_kwargs):
    def diags(model):
        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(
                quantity=quantity, write_timestamps=timestamps, **diag_kwargs
            )

    return populate_simulation(
        1,
        1,
        time_step_nbr=time_step_nbr,
        final_time=final_time,
        diag_options={
            "format": diag_format,
            "options": {"dir": f"{diag_dir}/{diag_format}", "mode": "overwrite"},
        },
        diags_fn=diags,
    )


@ddt
class SimulatorLifecycleTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None
        self.diag_dirs = [diag_dir]
        self.success = False

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        super().tearDown()

    def test_data_wrangler_after_run(self):
        self.simulator = Simulator(config()).run()
        self.assertIsNotNone(self.simulator.cpp_sim)
        self.assertAlmostEqual(self.simulator.currentTime(), final_time)
        self.assertTrue(self.simulator.data_wrangler().lvl0EM())
        self.success = True

    @data("phareh5", "pharevtkhdf")
    def test_diags_complete_after_run_without_reset(self, diag_format):
        import h5py

        flush_never = 0  # files stay open until the end of run()
        self.simulator = Simulator(config(diag_format, flush_every=flush_never)).run()

        ext = "h5" if diag_format == "phareh5" else "vtkhdf"
        with h5py.File(f"{diag_dir}/{diag_format}/EM_B.{ext}", "r") as h5:
            if diag_format == "phareh5":
                n_steps = len(h5["t"])
            else:
                n_steps = int(h5["VTKHDF/Steps"].attrs["NSteps"])
        self.assertEqual(n_steps, len(timestamps))
        self.success = True

    def test_next_setup_releases_previous_simulator(self):
        first = Simulator(config()).run()
        self.simulator = Simulator(config()).run()
        self.assertIsNone(first.cpp_sim)
        self.assertIsNotNone(self.simulator.cpp_sim)
        first = None  # its __del__ must not reset SAMRAI under the live simulator
        self.assertTrue(self.simulator.data_wrangler().lvl0EM())
        self.success = True

    def test_setup_raises_while_another_is_unfinished(self):
        first = Simulator(config()).initialize()  # live, not at its final time
        with self.assertRaises(RuntimeError):
            Simulator(config()).setup()
        self.assertIsNotNone(first.cpp_sim)  # untouched
        first.reset()
        self.simulator = Simulator(config()).run()
        self.success = True

    def test_dropped_simulator_is_released(self):
        first = Simulator(config()).initialize()
        first = None  # noqa: F841, no reference left: collected and reset
        self.simulator = Simulator(config()).run()
        self.success = True


if __name__ == "__main__":
    unittest.main()
