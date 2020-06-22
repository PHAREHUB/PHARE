#!/usr/bin/env python3
#
# formatted with black


from pybindlibs import cpp
from tests.simulator import create_simulator
from pyphare.data.wrangler import DataWrangler
import unittest, numpy as np

# TODO - validate data from somewhere!

class DataWranglerTest(unittest.TestCase):


    def test_1d(self):

        for interp in range(1, 4):

            self.dman, self.sim, self.hier = create_simulator(1, interp)
            self.dw = DataWrangler(self.sim, self.hier)

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

            del (
                self.dw,
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
