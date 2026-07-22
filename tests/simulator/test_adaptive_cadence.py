#
# Runtime coverage for the dump-cadence catch-up scheduler and adaptive-dt dumping:
#  - the catch-up scheduler keeps a timestamp-based cadence from freezing when several
#    scheduled times fall within a single coarse step (dt > gap between timestamps);
#  - a diagnostic scheduled exactly at final_time still fires under adaptive dt.
#
# test_time_step.py only checks config-level validation (never runs the simulator); the tests
# here actually run a few coarse steps and inspect the dumped h5 files.
#

import os
import unittest
import numpy as np

from pathlib import Path

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy.fromh5 import h5_filename_from, h5_time_grp_key
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest


def _h5_time_group_count(h5_filepath):
    import h5py  # see doc/conventions.md section 2.1.1

    with h5py.File(h5_filepath, "r") as h5_file:
        return len(h5_file[h5_time_grp_key].keys())


def setup_model(sim, ppc=10):
    def density(x):
        return 1.0

    def bx(x):
        return 1.0

    def by(x):
        return 0.0

    def bz(x):
        return 0.0

    def v(x):
        return 0.0

    def vth(x):
        return 0.1

    vvv = dict(vbulkx=v, vbulky=v, vbulkz=v, vthx=vth, vthy=vth, vthz=vth)

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1337},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)


out = "phare_outputs/adaptive_cadence"


def base_args(diagdir, **extra):
    args = dict(
        interp_order=1,
        time_step={"mode": "adaptive", "cfl": 0.8},
        final_time=100.0,  # generous: we control the number of steps manually
        boundary_types="periodic",
        cells=20,
        dl=0.3,
        diag_options=dict(
            format="phareh5", options=dict(dir=diagdir, mode="overwrite")
        ),
    )
    args.update(extra)
    return args


class AdaptiveCadenceTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(AdaptiveCadenceTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super().tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_write_timestamps_catch_up_does_not_freeze(self):
        """Regression test for the dump-cadence-freeze bug: when several scheduled write
        timestamps fall within a single coarse step, the scheduler must catch up (advance past
        all of them) instead of consuming one and freezing. A timestamp array much finer than
        the coarse step forces many due times per step; the cadence must keep firing every step
        rather than stalling after the first."""
        n_advances = 4
        # a timestamp grid far finer than any plausible dt forces the catch-up loop to consume
        # several scheduled times per step. final_time stays large enough that n_advances steps
        # don't reach it. Fed directly as write_timestamps (no period option involved).
        sim = self.simulation(**base_args(out + "/timestamps_catchup", final_time=1.0))
        setup_model(sim)
        ph.ElectromagDiagnostics(
            quantity="B", write_timestamps=np.arange(0.0, 1.0, 1e-4)
        )

        self.simulator = Simulator(sim).initialize()
        for _ in range(n_advances):
            self.simulator.advance()
        self.simulator.reset()

        diag_dir = sim.diag_options["options"]["dir"]
        diagInfo = next(iter(sim.diagnostics.values()))
        h5_filepath = os.path.join(diag_dir, h5_filename_from(diagInfo))
        # one dump per step (init + each advance): if the cadence had frozen after the
        # first step, this would be 1 instead of n_advances + 1
        self.assertEqual(_h5_time_group_count(h5_filepath), n_advances + 1)

    def test_dump_scheduled_exactly_at_final_time_is_not_dropped(self):
        """Regression test: under adaptive dt, timeStep() clamps to 0 on the very last step
        (currentTime()==endTime()). If dump() re-derived its timestep from a fresh timeStep()
        query instead of the dt that actually produced the current state, a diagnostic scheduled
        exactly at final_time would silently never fire (0 < 0 is always false)."""
        final_time = 0.2  # small: run() drives the whole thing to completion, few steps expected

        sim = self.simulation(
            **base_args(out + "/final_time_dump", final_time=final_time)
        )
        setup_model(sim)
        ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([final_time]))

        self.simulator = Simulator(sim)
        self.simulator.run()
        self.simulator = None  # run() already reset()

        diag_dir = sim.diag_options["options"]["dir"]
        diagInfo = next(iter(sim.diagnostics.values()))
        h5_filepath = os.path.join(diag_dir, h5_filename_from(diagInfo))
        self.assertEqual(_h5_time_group_count(h5_filepath), 1)


if __name__ == "__main__":
    startMPI()
    unittest.main()
