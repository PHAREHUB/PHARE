
#!/usr/bin/env python3
#
# formatted with black



from pybindlibs import cpp
from tests.diagnostic import dump_all_diags
from tests.simulator import populate_simulation
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from
import pyphare.pharein as ph
import unittest, os
from ddt import ddt, data



out = "phare_outputs/diagnostic_test"
diags = {"diag_options": {"format": "phareh5", "options": {"dir": out}}}


def dup(dic):
    dic.update(diags.copy())
    dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
    return dic

@ddt
class DiagnosticsTest(unittest.TestCase):

    _test_cases = (
      dup({
        "smallest_patch_size": 5,
        "largest_patch_size": 64}),
      dup({
        "smallest_patch_size": 65,
        "largest_patch_size": 65})
    )

    def __init__(self, *args, **kwargs):
        super(DiagnosticsTest, self).__init__(*args, **kwargs)
        self.simulator = None


    def _test_dump_diags_with_killing_dman_nd(self, dim, **input):

        for interp in range(1, 4):
            local_out = out + "_" + str(dim) + "_" + str(interp)
            input["diag_options"]["options"]["dir"] = local_out

            self.simulator = Simulator(populate_simulation(dim, interp, **input))
            self.simulator.initialize()
            self.simulator.diagnostics().dump(timestamp=0, timestep=1)

            for diagInfo in ph.globals.sim.diagnostics:
                # diagInfo.quantity starts with a / this interferes with os.path.join, hence   [1:]
                h5_file = os.path.join(local_out, (diagInfo.quantity + ".h5").replace('/', '_')[1:])
                self.assertTrue(os.path.exists(h5_file))
                print("hier", hierarchy_from(h5_filename=h5_file))

            self.simulator = None


    @data(*_test_cases)
    def test_dump_diags_with_killing_dman_1d(self, input):

        self._test_dump_diags_with_killing_dman_nd(1, **input)

    @data(*_test_cases)
    def test_dump_diags_with_killing_dman_2d(self, input):

        self._test_dump_diags_with_killing_dman_nd(2, **input)

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()


if __name__ == "__main__":
    unittest.main()
