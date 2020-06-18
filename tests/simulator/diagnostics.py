#!/usr/bin/env python3
#
# formatted with black


from pybindlibs import cpp
from tests.diagnostic import dump_all_diags
from tests.simulator import create_simulator
from pyphare.data.wrangler import DataWrangler
import unittest, os, shutil
from pyphare.pharesee.hierarchy import hierarchy_from


out = "phare_outputs/diagnostic_test"
diags = {"diag_options": {"format": "phareh5", "options": {"dir": out}}}
rank = cpp.mpi_rank()

def dup(dic):
    dic.update(diags.copy())
    dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
    return dic

class DiagnosticsTest(unittest.TestCase):

    def test_dump_diags_with_killing_dman_1d(self):

        for interp in range(1, 4):

            if rank == 0 and os.path.exists(out):
                shutil.rmtree(out)

            self.dman, self.sim, self.hier = create_simulator(1, interp, **dup({}))
            self.dman.dump(timestamp=0, timestep=1)
            em_b_file = os.path.join(out, "EM_B.h5")
            self.assertTrue(os.path.exists(os.path.join(out, "EM_B.h5")))
            hier = hierarchy_from(em_b_file)
            print("hier", hier)

            del (
                self.dman,
                self.sim,
                self.hier,
            )
            cpp.reset()

    def tearDown(self):
        for k in ["dw", "dman", "sim", "hier"]:
            if hasattr(self, k):
                v = getattr(self, k)
                del v  # blocks segfault on test failure, could be None
        cpp.reset()


if __name__ == "__main__":
    unittest.main()
