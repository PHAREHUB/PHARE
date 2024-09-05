#
#

import unittest
import numpy as np
import pyphare.pharein as ph

from pyphare.simulator.simulator import Simulator

from copy import deepcopy
from tests.simulator import SimulatorTest

simArgs = dict(
    time_step_nbr=1,  # avoid regrid for refinement boxes https://github.com/LLNL/SAMRAI/issues/199
    time_step=0.001,
    boundary_types="periodic",
    cells=np.array([20]),
    dl=0.3,
)


class SimulatorValidation(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(SimulatorValidation, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(SimulatorValidation, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        ph.global_vars.sim = None

    def test_no_numpy_on_serialization(self):
        import dill
        import codecs

        simput = deepcopy(simArgs)
        sim = ph.Simulation(**deepcopy(simArgs))

        # check is np array before serialization
        assert isinstance(sim.cells, np.ndarray)
        hex = ph.simulation.serialize(sim)

        def check_numpify_recursive(owner, key):
            obj = getattr(owner, key)
            if obj is None:
                return 0
            if hasattr(obj, "__dict__"):
                for k in obj.__dict__:
                    check_numpify_recursive(obj, k)
            else:
                assert not isinstance(obj, np.ndarray)
                return 1
            return 0

        desim = dill.loads(codecs.decode(hex, "hex"))
        checks = 0
        for k in desim.__dict__:
            checks += check_numpify_recursive(desim, k)
        assert checks > 1  # we asserted something


if __name__ == "__main__":
    unittest.main()
