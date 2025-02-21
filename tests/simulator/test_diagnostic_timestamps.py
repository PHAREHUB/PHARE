#!/usr/bin/env python3

from pyphare.cpp import cpp_lib

cpp = cpp_lib()
import os
import unittest

import h5py
import numpy as np
import pyphare.pharein as ph
from ddt import data, ddt
from pyphare.core.box import Box1D
from pyphare.pharein import ElectromagDiagnostics, ElectronModel
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.hierarchy.fromh5 import h5_filename_from, h5_time_grp_key
from pyphare.pharesee.hierarchy.hierarchy import format_timestamp
from pyphare.simulator.simulator import Simulator


def setup_model(ppc):
    def density(x):
        return 1.0

    def by(x):
        return 0.0

    def bz(x):
        return 0.0

    def bx(x):
        return 1.0

    def vx(x):
        return 0.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vthx(x):
        return 1.00

    def vthy(x):
        return 1.00

    def vthz(x):
        return 1.00

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
    }

    model = ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "nbr_part_per_cell": ppc},
    )
    ElectronModel(closure="isothermal", Te=0.12)
    return model


out = "phare_outputs/diagnostic_ts_test/"
simArgs = {
    "smallest_patch_size": 10,
    "largest_patch_size": 10,
    "time_step_nbr": 1e5,  # is sufficient based on https://github.com/PHARCHIVE/test_snippets/blob/main/numeric/double/increment_error.cpp
    "time_step": 0.001,
    "boundary_types": "periodic",
    "cells": 10,
    "dl": 0.2,
    "diag_options": {"format": "phareh5", "options": {"dir": out, "mode": "overwrite"}},
    "strict": True,
}


@ddt
class DiagnosticsTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(DiagnosticsTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def setUp(self):
        from pyphare.simulator.simulator import startMPI

        startMPI()

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None

    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]

    @data(
        ({"L0": {"B0": Box1D(10, 14), "B1": Box1D(15, 19)}}),
    )
    def test_hierarchy_timestamp_cadence(self, refinement_boxes):
        """
        this test checks diagnostics are dumped at the correct timestamps

        - only B is dumped
        - we dump every 2dt, or every 3dt
        - we check:
            - diag file exists
            - expected number of times in diag h5
            - expected values of times in diag h5
            - that we can skip init dump (t=t0)


        This test only runs in 1D. other Dims are unnecessary
        """
        dim = refinement_boxes["L0"]["B0"].ndim

        time_step = 0.001
        # time_step_nbr chosen to force diagnostics dumping double imprecision cadence calculations accuracy testing
        time_step_nbr = 101
        final_time = time_step * time_step_nbr

        for skip_init in [0, 1]:  # 1 = skip init dumps
            for diag_cadence in [2, 3]:
                simInput = simArgs.copy()
                diag_outputs = f"phare_outputs_hierarchy_timestamp_cadence_{dim}_{self.ddt_test_id()}_{diag_cadence}"
                simInput["diag_options"]["options"]["dir"] = diag_outputs
                simInput["time_step_nbr"] = time_step_nbr

                ph.global_vars.sim = None
                simulation = ph.Simulation(**simInput)
                setup_model(10)

                timestamps = np.arange(0, final_time, time_step * i)[skip_init:]
                for quantity in ["B"]:
                    ElectromagDiagnostics(
                        quantity=quantity,
                        write_timestamps=timestamps,
                        flush_every=ElectromagDiagnostics.h5_flush_never,
                    )

                Simulator(simulation).run()

                for diagname, diagInfo in simulation.diagnostics.items():
                    h5_filename = os.path.join(diag_outputs, h5_filename_from(diagInfo))
                    self.assertTrue(os.path.exists(h5_filename))

                    hier = hierarchy_from(h5_filename=h5_filename)

                    time_hier_keys = list(hier.time_hier.keys())
                    self.assertEqual(len(time_hier_keys), len(timestamps))

                    for diag_cadence, timestamp in enumerate(time_hier_keys):
                        self.assertEqual(
                            format_timestamp(timestamps[diag_cadence]), timestamp
                        )


if __name__ == "__main__":
    unittest.main()
