#!/usr/bin/env python3
#
# formatted with black


from pybindlibs import cpp
from tests.simulator import populate_simulation
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
            print("\n", self.dw.lvl0PopFluxs())
            print("\n", self.dw.lvl0EM())

            self.simulator = None

    def tearDown(self):
        del self.dw
        if self.simulator is not None:
            self.simulator.reset()



if __name__ == "__main__":
    unittest.main()
