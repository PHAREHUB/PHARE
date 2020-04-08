#!/usr/bin/env python3
#
# formatted with black

from tests.simulator.py import InitValueValidation

import numpy as np

from phare.data.wrangler import primal_datasets

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

    def checkEMDataEqualsUserFunctions(self, ems):
        for diag in ems:
            patch_level0 = diag.levels[0]
            for patch in list(patch_level0.patches.values()):
                for xyz in ["x", "y", "z"]:
                    hdf5_data = patch.patch_data.data(xyz)
                    nGhosts = patch.patch_data.nGhosts(xyz)
                    fn_name = patch.patch_data.quantity_key.lower() + xyz
                    fn = self.getSimulation().model.model_dict[fn_name]
                    xorigin = patch.origin[0]
                    cell_width = patch.patch_level.cell_width("x")
                    is_primal = primal_datasets[fn_name]

                    x_start = xorigin if is_primal else xorigin + cell_width / 2
                    physical_hdf5_dataset = hdf5_data[nGhosts:-nGhosts]

                    fn_dataset = _as_np_float32([
                        fn(x_start + (cell_width * i))
                        for i in range(len(physical_hdf5_dataset))
                    ])

                    self.assertTrue(np.allclose(fn_dataset, physical_hdf5_dataset, atol=1e-05))


if __name__ == "__main__":
    import unittest

    unittest.main()
