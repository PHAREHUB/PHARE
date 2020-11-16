
#!/usr/bin/env python3


from pybindlibs import cpp
from tests.diagnostic import dump_all_diags
from tests.simulator import populate_simulation
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy import hierarchy_from
import pyphare.pharein as ph
import unittest
import os
import h5py
import numpy as np
from ddt import ddt, data

from tests.simulator.config import project_root


def setup_model():
    def density(*xyz):
        return 1.

    def by(*xyz):
        x = xyz[0]
        L = ph.global_vars.sim.simulation_domain()
        return 0.1*np.cos(2*np.pi*x/L[0])

    def bz(*xyz):
        x = xyz[0]
        L = ph.global_vars.sim.simulation_domain()
        return 0.1*np.sin(2*np.pi*x/L[0])

    def bx(*xyz):
        return 1.

    def vx(*xyz):
        return 0.

    def vy(*xyz):
        x = xyz[0]
        L = ph.global_vars.sim.simulation_domain()
        return 0.1*np.cos(2*np.pi*x/L[0])

    def vz(*xyz):
        x = xyz[0]
        L = ph.global_vars.sim.simulation_domain()
        return 0.1*np.sin(2*np.pi*x/L[0])

    def vthx(*xyz):
        return 0.01

    def vthy(*xyz):
        return 0.01

    def vthz(*xyz):
        return 0.01

    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz
    }

    model = ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 1337}}
    )
    ElectronModel(closure="isothermal", Te=0.12)
    return model


out = "phare_outputs/diagnostic_test/"
simArgs = {
  "time_step_nbr":30000,
  "final_time":30.,
  "boundary_types":"periodic",
  "cells":40,
  "dl":0.3,
  "max_nbr_levels":2,
  "diag_options": {"format": "phareh5", "options": {"dir": out, "mode":"overwrite", "fine_dump_lvl_max": 10}}
}

def dup(dic):
    dic.update(simArgs.copy())
    return dic


@ddt
class DiagnosticsTest(unittest.TestCase):

    _test_cases = (
      dup({
        "smallest_patch_size": 10,
        "largest_patch_size": 20}),
      dup({
        "smallest_patch_size": 20,
        "largest_patch_size": 20}),
      dup({
        "smallest_patch_size": 40,
        "largest_patch_size": 40})
    )

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


    def _test_dump_diags(self, dim, **simInput):
        test_id = self.ddt_test_id()

        # configure simulation dim sized values
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(dim)]
        b0 = [[10 for i in range(dim)], [19 for i in range(dim)]]
        simInput["refinement_boxes"] = {"L0": {"B0": b0}}

        for interp in range(1, 4):
            print("_test_dump_diags dim/interp:{}/{}".format(dim, interp))

            local_out = f"{out}_dim{dim}_interp{interp}_mpi_n_{cpp.mpi_size()}_id{test_id}"
            simInput["diag_options"]["options"]["dir"] = local_out

            simulation = ph.Simulation(**simInput)
            self.assertTrue(len(simulation.cells) == dim)

            dump_all_diags(setup_model().populations)
            self.simulator = Simulator(simulation).initialize().advance()

            for diagInfo in ph.global_vars.sim.diagnostics:
                # diagInfo.quantity starts with a / this interferes with os.path.join, hence   [1:]
                h5_filename = os.path.join(local_out, (diagInfo.quantity + ".h5").replace('/', '_')[1:])
                print("h5_filename", h5_filename)

                h5_file = h5py.File(h5_filename, "r")
                self.assertTrue("t0.000000" in h5_file) #    init dump

                self.assertTrue("t0.000100" in h5_file)
                self.assertTrue("pl1" in h5_file["t0.000100"])
                self.assertFalse("pl0" in h5_file["t0.000100"])

                self.assertTrue("t0.001000" in h5_file) # advance dump


                # SEE https://github.com/PHAREHUB/PHARE/issues/275
                if dim == 1: # REMOVE WHEN PHARESEE SUPPORTS 2D
                    self.assertTrue(os.path.exists(h5_filename))
                    hier = hierarchy_from(h5_filename=h5_filename)
                    if h5_filename.endswith("domain.h5"):
                        for patch in hier.level(0).patches:
                            for qty_name, pd in patch.patch_datas.items():
                                splits = pd.dataset.split(ph.global_vars.sim)
                                self.assertTrue(splits.size() == pd.dataset.size() * 2)
                                print("splits.iCell", splits.iCells)
                                print("splits.delta", splits.deltas)
                                print("splits.weight", splits.weights)
                                print("splits.charge", splits.charges)
                                print("splits.v", splits.v)

            self.simulator = None
            ph.global_vars.sim = None


if __name__ == "__main__":
    unittest.main()
