#!/usr/bin/env python3
#
# formatted with black


from pyphare.cpp import cpp_lib
cpp = cpp_lib()
from tests.simulator import populate_simulation
import numpy as np
from pyphare.simulator.simulator import Simulator
import unittest

# TODO - validate data from somewhere!

class DataWranglerTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(DataWranglerTest, self).__init__(*args, **kwargs)
        self.dw = None
        self.simulator = None

    def test_1d(self):

        for interp in range(1, 4):

            self.simulator = Simulator(populate_simulation(1, interp))
            self.simulator.initialize()
            self.dw = self.simulator.data_wrangler()

            print("\n", self.dw.lvl0IonDensity())
            print("\n", self.dw.lvl0BulkVelocity())
            print("\n", self.dw.lvl0PopDensity())
            print("\n", self.dw.lvl0PopFluxes())
            print("\n", self.dw.lvl0EM())

            for pop, particles in self.dw.getPatchLevel(0).getParticles().items():
                for key, patches in particles.items():
                    for patch in patches:
                        self.assertTrue(isinstance(patch.lower, np.ndarray))
                        self.assertTrue(isinstance(patch.upper, np.ndarray))

            self.simulator = None

    def tearDown(self):
        del self.dw
        if self.simulator is not None:
            self.simulator.reset()



if __name__ == "__main__":
    unittest.main()
