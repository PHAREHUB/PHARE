#!/usr/bin/env python3
#
# formatted with black

from phare import cpp
from tests.simulator.py import create_simulator
from phare.data.wrangler import DataWrangler

import unittest

# TODO - validate data from somewhere!


class DataWranglerTest(unittest.TestCase):
    def tearDown(self):
        cpp.reset()

    def test_1d(self):

        for interp in range(1, 4):

            self.dman, self.sim, self.hier = create_simulator(1, interp)
            dw = DataWrangler(self.sim, self.hier)

            print("\n", dw.lvl0IonDensity())
            print("\n", dw.lvl0BulkVelocity())
            print("\n", dw.lvl0PopDensity())
            print("\n", dw.lvl0PopFluxs())
            print("\n", dw.lvl0EM())

            del (
                dw,
                self.dman,
                self.sim,
                self.hier,
            )
            cpp.reset()


if __name__ == "__main__":
    unittest.main()
