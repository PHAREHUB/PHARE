#!/usr/bin/env python3

import unittest
import numpy as np

from pyphare import cpp
from pyphare.simulator.simulator import Simulator
from tests.simulator import populate_simulation


from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim


class DataWranglerTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(DataWranglerTest, self).__init__(*args, **kwargs)
        self.dw = None
        self.simulator = None

    def test(self):
        ndim, interp = 2, 1
        self.simulator = Simulator(populate_simulation(ndim, interp))
        self.simulator.initialize()

        hier = None
        hier = hierarchy_from_sim(self.simulator, qty="EM_B_x", hier=hier, sync=True)
        hier = hierarchy_from_sim(self.simulator, qty="EM_B_y", hier=hier, sync=True)
        hier = hierarchy_from_sim(self.simulator, qty="density", hier=hier, sync=True)

        if cpp.mpi_rank() == 0:
            hier.plot(
                filename="data_wrangler.Bx.png",
                qty="EM_B_x",
                plot_patches=True,
                levels=(0,),
            )
            hier.plot(
                filename="data_wrangler.By.png",
                qty="EM_B_y",
                plot_patches=True,
                levels=(0,),
            )
            # hier.plot( # fails?
            #     filename="data_wrangler.Ni.png",
            #     qty="density",
            #     plot_patches=True,
            #     levels=(0,),
            # )

    def tearDown(self):
        del self.dw
        if self.simulator is not None:
            self.simulator.reset()


if __name__ == "__main__":
    unittest.main()
