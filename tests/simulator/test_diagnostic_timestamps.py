
#!/usr/bin/env python3

from pyphare.cpp import cpp_lib
cpp = cpp_lib()
from tests.diagnostic import dump_all_diags
from tests.simulator import populate_simulation
from pyphare.pharein import ElectronModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics, ParticleDiagnostics
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy import hierarchy_from, h5_filename_from, h5_time_grp_key
import pyphare.pharein as ph
import unittest
import os
import h5py
import numpy as np
from ddt import ddt, data

from tests.simulator.config import project_root


def setup_model(ppc):

    def density(x):
        return 1.

    def by(x):
        return 0.

    def bz(x):
        return 0.

    def bx(x):
        return 1.

    def vx(x):
        return 0.

    def vy(x):
        return 0.

    def vz(x):
        return 0.

    def vthx(x):
        return 1.00

    def vthy(x):
        return 1.00

    def vthz(x):
        return 1.00


    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz
    }

    model = ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, **vvv, "nbr_part_per_cell":ppc}
    )
    ElectronModel(closure="isothermal", Te=0.12)
    return model


out = "phare_outputs/diagnostic_ts_test/"
simArgs = {
  "smallest_patch_size": 10, "largest_patch_size": 10,
  "time_step_nbr":1e5, # is sufficient based on https://github.com/PHARCHIVE/test_snippets/blob/main/numeric/double/increment_error.cpp
  "time_step": .001,
  "boundary_types":"periodic",
  "cells":10,
  "dl":0.2,
  "diag_options": {"format": "phareh5", "options": {"dir": out, "mode":"overwrite"}}
}


@ddt
class DiagnosticsTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(DiagnosticsTest, self).__init__(*args, **kwargs)
        startMPI()
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None

    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]

    def test_dump_diags_timestamps(self):
        print("test_dump_diags dim/interp:{}/{}".format(1, 1))

        simulation = ph.Simulation(**simArgs.copy())
        sim = simulation

        dump_every = 1
        timestamps = np.arange(0, sim.final_time + sim.time_step, dump_every*sim.time_step)
        setup_model(10)

        for quantity in ["B"]:
            ElectromagDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps,
                flush_every=ElectromagDiagnostics.h5_flush_never,
            )

        Simulator(simulation).run()

        def make_time(stamp):
            return "{:.10f}".format(stamp)

        for diagInfo in ph.global_vars.sim.diagnostics:
            h5_filename = os.path.join(out, h5_filename_from(diagInfo))
            self.assertTrue(os.path.exists(h5_filename))

            h5_file = h5py.File(h5_filename, "r")

            for timestamp in timestamps:
                self.assertIn(make_time(timestamp), h5_file[h5_time_grp_key])



if __name__ == "__main__":
    unittest.main()

