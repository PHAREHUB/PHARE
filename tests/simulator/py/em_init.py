#!/usr/bin/env python3
#
# formatted with black

from tests.simulator.py import InitValueValidation

import numpy as np


class EMInitValidation(InitValueValidation):
    is_primal = {"bx": 1, "by": 0, "bz": 0, "ex": 0, "ey": 1, "ez": 1}

    def test_1d(self):
        from phare.pp.diagnostics import _EM
        from tests.diagnostic import dump_all_diags

        diag_out_dir = "phare_outputs/em_init"
        dic = InitValueValidation.diag_options(diag_out_dir)
        dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
        diags = self._simulate_diagnostics(dim=1, interp=1, input=dic)
        self._checkEMFromFN(diags[_EM.__name__])

    def _checkEMFromFN(self, ems):
        truncation, tolerance = 5, 1e-4  # the nearness is not great

        for diag in ems:
            patch_level0 = diag.levels[0]
            for patch in patch_level0.patches:
                for xyz in ["x", "y", "z"]:
                    hdf5_data = patch.patch_data.get()[xyz]
                    nGhosts = patch.patch_data.nGhosts(xyz)
                    fn_name = patch.patch_data.key.lower() + xyz
                    fn = diag.sim.model.model_dict[fn_name]
                    xorigin = patch.origin[0]
                    cell_width = patch.patch_level.cell_width("x")
                    is_primal = EMInitValidation.is_primal[fn_name]
                    x_pos = xorigin if is_primal else xorigin + cell_width / 2
                    physical_dataset = [
                        self.truncate(x, truncation)
                        for x in hdf5_data[nGhosts : int(nGhosts * -1)]
                    ]
                    fn_x = [
                        self.truncate(fn(x_pos + ((cell_width) * i)), truncation)
                        for i in range(len(physical_dataset))
                    ]
                    np.testing.assert_allclose(fn_x, physical_dataset, rtol=tolerance)


if __name__ == "__main__":
    import unittest

    unittest.main()
