#!/usr/bin/env python3
#
# formatted with black

from tests.simulator.py import InitValueValidation

import numpy as np

from phare.core.gridlayout import yee_element_is_primal


def _as_np_float32(array):
    return np.array(array, dtype=np.float32)


class EMInitValidation(InitValueValidation):
    def test_1d(self):
        from phare.pp.diagnostics import _EMPatchData
        from tests.diagnostic import dump_all_diags

        diag_out_dir = "phare_outputs/em_init"
        simInput = InitValueValidation.diag_options(diag_out_dir)
        simInput.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
        diags = self.runAndDump(dim=1, interp=1, input=simInput)
        self.checkEMDataEqualsUserFunctions(diags[_EMPatchData.__name__])

    def checkEMDataEqualsUserFunctions(self, emsDiags):
        for diag in emsDiags:
            patch_level0 = diag.levels[0]
            for patch in patch_level0.patches_list():
                for xyz in ["x", "y", "z"]:
                    hdf5_data = patch.data(xyz)
                    nGhosts = patch.patch_data.nGhosts(xyz)
                    fn_name = patch.patch_data.quantity_key.lower() + xyz
                    fn = self.getSimulation().model.model_dict[fn_name]
                    xorigin = patch.origin[0]
                    cell_width = patch.patch_level.cell_width("x")
                    is_primal = yee_element_is_primal(fn_name, "x")

                    x_start = xorigin if is_primal else xorigin + cell_width / 2
                    physical_hdf5_dataset = hdf5_data[nGhosts:-nGhosts]

                    x = np.arange(len(physical_hdf5_dataset)) * cell_width + x_start
                    expected_values = _as_np_float32(fn(x))

                    self.assertTrue(
                        np.allclose(expected_values, physical_hdf5_dataset, atol=1e-05)
                    )


if __name__ == "__main__":
    import unittest

    unittest.main()
