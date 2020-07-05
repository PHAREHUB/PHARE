#!/usr/bin/env python3
#
# formatted with black


from pybindlibs import cpp
from tests.diagnostic import dump_all_diags
from tests.simulator import create_simulator
from pyphare.pharesee.hierarchy import hierarchy_from
import pyphare.pharein as ph
import unittest, os, shutil
from ddt import ddt, data


out = "phare_outputs/diagnostic_test"
diags = {"diag_options": {"format": "phareh5", "options": {"dir": out}}}
rank = cpp.mpi_rank()


def dup(dic):
    dic.update(diags.copy())
    dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
    return dic

@ddt
class DiagnosticsTest(unittest.TestCase):

    @data(
      dup({
        "smallest_patch_size": 5,
        "largest_patch_size": 64}),
      dup({
        "smallest_patch_size": 65,
        "largest_patch_size": 65})
    )
    def test_dump_diags_with_killing_dman_1d(self, input):

        dim = 1
        for interp in range(1, 4):
            local_out = out + "_" + str(dim) + "_" + str(interp)
            input["diag_options"]["options"]["dir"] = local_out

            self.dman, self.sim, self.hier = create_simulator(1, interp, **input)
            self.dman.dump(timestamp=0, timestep=1)

            for diagInfo in ph.globals.sim.diagnostics:
                # diagInfo.quantity starts with a / this interferes with os.path.join, hence   [1:]
                h5_file = os.path.join(local_out, (diagInfo.quantity + ".h5").replace('/', '_')[1:])
                self.assertTrue(os.path.exists(h5_file))
                print("hier", hierarchy_from(h5_filename=h5_file))

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
