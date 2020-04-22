#!/usr/bin/env python3
#
# formatted with black

from phare import cpp

from phare.data.wrangler import DataWrangler
from phare.core.gridlayout import yee_element_is_primal
from phare.pp.diagnostics import Diagnostics
from phare.pp.diagnostic.patch import aggregate_level0_patch_domain

from tests.simulator.py import InitValueValidation
from tests.simulator.py import create_simulator
from tests.diagnostic import dump_all_diags

import unittest
import numpy as np


def _as_np_float32(array):
    return np.array(array, dtype=np.float32)

diag_out_dir = "phare_outputs/data_wrangler"

class DataWranglerTest(unittest.TestCase):
    def tearDown(self):
        cpp.reset()  # test failure = hang forever without this

    def test_1d(self):
        dim = 1
        for interp in range(1, 4):
            simInput = InitValueValidation.diag_options(diag_out_dir)
            simInput.update(
                {"diags_fn": lambda model: dump_all_diags(model.populations)}
            )

            self.dman, self.sim, self.hier = create_simulator(1, interp, **simInput)
            self.dman.dump(timestamp=0, timestep=1)
            del self.dman  # flush diagnostics

            diags = Diagnostics(diag_out_dir)
            dw = DataWrangler(self.sim, self.hier)

            self.check_field(dw.lvl0IonDensity(), diags.ionDensity(), "density")
            self.check_fields(dw.lvl0PopDensity(), diags.popDensities(), "density")

            self.check_vec_field(
                dw.lvl0BulkVelocity(), diags.ionVelocity(), "ions_bulkVel"
            )
            self.check_vec_fields(dw.lvl0PopFluxs(), diags.popFluxes(), "flux")

            dw_EB = dw.lvl0EM()
            self.checkEM(dw, dw_EB, diags.E(), "EM_E")
            self.checkEM(dw, dw_EB, diags.B(), "EM_B")

            diags.close()
            del (dw, self.sim, self.hier)
            cpp.reset()

    def check_field(self, dw_field, diag, data_name):
        patch_level0 = diag.levels[0]
        dw_field = _as_np_float32(dw_field)
        diag_field = aggregate_level0_patch_domain(
            patch_level0, shared_patch_border=True
        )[data_name]
        self.assertTrue(np.array_equiv(dw_field, diag_field))

    def check_fields(self, dw_fields, pop_diags, data_name):
        for pop_diag in pop_diags:
            patch_level0 = pop_diag.levels[0]
            pop_name = patch_level0.patches_list()[0].patch_data.pop_name
            self.check_field(dw_fields[pop_name], pop_diag, data_name)

    def checkEM(self, dw, dw_EB, diags_EM, data_name):
        dw_vecf = dw_EB[data_name]
        patch_level0 = diags_EM.levels[0]
        for component_name in patch_level0.data_names():
            dataset_name = data_name + "_" + component_name
            dw_vecf_xyz = np.array(dw_vecf[dataset_name], dtype=np.float32)
            is_primal = yee_element_is_primal(data_name[-1] + component_name)
            diag_vecf = aggregate_level0_patch_domain(
                patch_level0, shared_patch_border=is_primal, ds_names=[component_name]
            )
            self.assertTrue(np.array_equiv(dw_vecf_xyz, diag_vecf[component_name]))

    def check_vec_field(self, dw_vecf, diag, data_name):
        patch_level0 = diag.levels[0]
        diag_vecf = aggregate_level0_patch_domain(
            patch_level0, shared_patch_border=True
        )
        for component_name in patch_level0.data_names():
            dw_vecf_xyz = np.array(
                dw_vecf[data_name + "_" + component_name], dtype=np.float32
            )
            self.assertTrue(np.array_equiv(dw_vecf_xyz, diag_vecf[component_name]))

    def check_vec_fields(self, dw_vecfs, pop_diags, pop_data_name):
        for pop_diag in pop_diags:
            patch_level0 = pop_diag.levels[0]
            pop_name = patch_level0.patches_list()[0].patch_data.pop_name
            self.check_vec_field(
                dw_vecfs[pop_name], pop_diag, pop_name + "_" + pop_data_name
            )


if __name__ == "__main__":
    unittest.main()
