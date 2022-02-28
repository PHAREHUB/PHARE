#!/usr/bin/env python3

from pyphare.cpp import cpp_lib
cpp = cpp_lib()

import os
import h5py
import unittest
import numpy as np
from ddt import ddt, data

from tests.diagnostic import dump_all_diags
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharein.simulation import supported_dimensions
from pyphare.pharesee.hierarchy import hierarchy_from, h5_filename_from, h5_time_grp_key
import pyphare.pharein as ph


def setup_model(ppc=100):

    def density(x):
        return 1.

    def S(x,x0,l):
        return 0.5*(1+np.tanh((x-x0)/l))

    def bx(x):
        return 0.

    def by(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()[0]
        v1=-1
        v2=1.
        return v1 + (v2-v1)*(S(x,L*0.25,1) -S(x, L*0.75, 1))

    def bz(x):
        return 0.5

    def b2(x):
        return bx(x)**2 + by(x)**2 + bz(x)**2

    def T(x):
        K = 1
        return 1/density(x)*(K - b2(x)*0.5)

    def vx(x):
        return 2.

    def vy(x):
        return 0.

    def vz(x):
        return 0.

    def vthx(x):
        return T(x)

    def vthy(x):
        return T(x)

    def vthz(x):
        return T(x)

    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz
    }

    model = ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"mass":1, "charge": 1, "density": density, **vvv, "nbr_part_per_cell":ppc, "init": {"seed": 1337}},
        alpha={"mass":4, "charge": 1, "density": density, **vvv, "nbr_part_per_cell":ppc, "init": {"seed": 2334}},
    )
    ElectronModel(closure="isothermal", Te=0.12)
    return model


out = "phare_outputs/tagging_test/"
simArgs = {
  "time_step_nbr":30000,
  "final_time":30.,
  "boundary_types":"periodic",
  "cells":200,
  "dl":0.3,
  "refinement":"tagging",
  "max_nbr_levels": 3,
  "diag_options": {"format": "phareh5", "options": {"dir": out, "mode":"overwrite", "fine_dump_lvl_max": 10}}
}

def dup(dic):
    dic.update(simArgs.copy())
    return dic


@ddt
class TaggingTest(unittest.TestCase):

    _test_cases = (
      dup({
        "smallest_patch_size": 10,
        "largest_patch_size": 20}),
      dup({
        "smallest_patch_size": 20,
        "largest_patch_size": 20}),
      dup({
        "smallest_patch_size": 20,
        "largest_patch_size": 40})
    )

    def __init__(self, *args, **kwargs):
        super(TaggingTest, self).__init__(*args, **kwargs)
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


    @data(*_test_cases)
    def test_tagging(self, simInput):
        # UPDATE when 2d tagging is finished
        for ndim in [1]: #supported_dimensions():
            self._test_dump_diags(ndim, **simInput)

    def _test_dump_diags(self, dim, **simInput):
        test_id = self.ddt_test_id()
        for key in ["cells", "dl", "boundary_types"]:
            simInput[key] = [simInput[key] for d in range(dim)]

        for interp in range(1, 4):
            local_out = f"{out}_dim{dim}_interp{interp}_mpi_n_{cpp.mpi_size()}_id{test_id}"
            simInput["diag_options"]["options"]["dir"] = local_out

            simulation = ph.Simulation(**simInput)
            self.assertTrue(len(simulation.cells) == dim)

            dump_all_diags(setup_model().populations)
            self.simulator = Simulator(simulation).initialize().advance().reset()

            self.assertTrue(any([diagInfo.quantity.endswith("tags") for diagInfo in ph.global_vars.sim.diagnostics]))

            checks = 0
            found = 0
            for diagInfo in ph.global_vars.sim.diagnostics:
                h5_filepath = os.path.join(local_out, h5_filename_from(diagInfo))
                self.assertTrue(os.path.exists(h5_filepath))

                h5_file = h5py.File(h5_filepath, "r")
                self.assertTrue("0.0000000000" in h5_file[h5_time_grp_key]) # init dump
                n_patches = len(list(h5_file[h5_time_grp_key]["0.0000000000"]["pl0"].keys()))

                if h5_filepath.endswith("tags.h5"):
                    found = 1
                    hier = hierarchy_from(h5_filename=h5_filepath)
                    patches = hier.level(0).patches
                    tag_found = 0
                    for patch in patches:
                        self.assertTrue(len(patch.patch_datas.items()))
                        for qty_name, pd in patch.patch_datas.items():
                            self.assertTrue((pd.dataset[:] >= 0).all())
                            self.assertTrue((pd.dataset[:] <  2).all())
                            tag_found |= (pd.dataset[:] == 1).any()
                        checks += 1

            self.assertEqual(found, 1)
            self.assertEqual(tag_found, 1)
            self.assertEqual(checks, n_patches)
            self.simulator = None
            ph.global_vars.sim = None


if __name__ == "__main__":
    unittest.main()
