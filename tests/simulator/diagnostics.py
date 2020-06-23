#!/usr/bin/env python3
#
# formatted with black


from pybindlibs import cpp
from tests.diagnostic import dump_all_diags
from tests.simulator import populate_simulation
from pyphare.simulator.simulator import Simulator
import unittest, os, shutil
from pyphare.pharesee.hierarchy import hierarchy_from


out = "phare_outputs/diagnostic_test"
diags = {"diag_options": {"format": "phareh5", "options": {"dir": out}}}


def dup(dic):
    dic.update(diags.copy())
    dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
    return dic

class DiagnosticsTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(DiagnosticsTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def test_dump_diags_with_killing_dman_1d(self):
        dim = 1
        for interp in range(1, 4):
            simput = dup({})
            # MPI rank only available after first "simulator.initialize()"
            # new files per output to verify each without deleting previous
            local_out = out + "_" + str(dim) + "_" + str(interp)
            simput["diag_options"]["options"]["dir"] = local_out
            self.simulator = Simulator(populate_simulation(dim, interp, **simput))
            self.simulator.initialize()

            self.simulator.diagnostics().dump(timestamp=0, timestep=1)

            em_b_file = os.path.join(local_out, "EM_B.h5")
            self.assertTrue(os.path.exists(os.path.join(local_out, "EM_B.h5")))

            print("hier", hierarchy_from(em_b_file))

            self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()

if __name__ == "__main__":
    unittest.main()
