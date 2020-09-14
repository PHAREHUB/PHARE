
#!/usr/bin/env python3


from pybindlibs import cpp
from tests.diagnostic import dump_all_diags
from tests.simulator import NoOverwriteDict, populate_simulation
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from
import pyphare.pharein as ph
import unittest, os, h5py, shutil
from ddt import ddt, data

from tests.simulator.config import project_root


out = "phare_outputs/diagnostic_test/"
diags = {"diag_options": {"format": "phareh5", "options": {"dir": out}}}


def dup(dic):
    dic = NoOverwriteDict(dic)
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
        Simulator.startMPI() # so we cna delete previous diags only on mpi_rank 0


    def _test_dump_diags_with_killing_dman_nd(self, dim, **simInput):
        print("_test_dump_diags_with_killing_dman_nd {}".format(dim))
        for interp in range(1, 2):
            local_out = out + str(dim) + "_" + str(interp)
            simInput["diag_options"]["options"]["dir"] = local_out

            # delete previous diags / can't truncate
            if cpp.mpi_rank() == 0 and os.path.exists(local_out):
                shutil.rmtree(local_out)

            self.simulator = Simulator(populate_simulation(dim, interp, **simInput))
            self.simulator.initialize()
            # self.simulator.advance() causes crash under MPI
            # EXCEPTION CAUGHT: The geometry on the patch is not set, please verify your configuration

            for diagInfo in ph.global_vars.sim.diagnostics:
                # diagInfo.quantity starts with a / this interferes with os.path.join, hence   [1:]
                h5_filename = os.path.join(local_out, (diagInfo.quantity + ".h5").replace('/', '_')[1:])
                print("h5_filename", h5_filename)

                h5_file = h5py.File(h5_filename, "r")
                self.assertTrue("t0.000000" in h5_file)
                # self.assertTrue("t0.001000" in h5_file) # fails under MPI runs

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


    @data(*_test_cases)
    def test_dump_diags_with_killing_dman_1d(self, simInput):
        self._test_dump_diags_with_killing_dman_nd(1, **simInput)


    @data(*_test_cases)
    def test_dump_diags_with_killing_dman_2d(self, simInput):
        self._test_dump_diags_with_killing_dman_nd(2, **simInput)


    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()


if __name__ == "__main__":
    unittest.main()
