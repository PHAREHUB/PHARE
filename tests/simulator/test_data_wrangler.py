#!/usr/bin/env python3
"""
the data wrangler and hierarchy_from_sim on a live simulation, the wrangler is also used
directly by test_patch_overlap.py and test_multi_species_momentum.py
"""

import unittest
from pathlib import Path

import numpy as np

from pyphare import cpp
from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim
from pyphare.simulator.simulator import Simulator
from tests.simulator import populate_simulation

ndim, interp = 2, 1
cells = 20  # basicSimulatorArgs
ppc = 100  # makeBasicModel / defaultPopulationSettings
out = Path("phare_outputs/test_data_wrangler")


class DataWranglerTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None

    def test_hierarchy_from_sim(self):
        self.simulator = Simulator(populate_simulation(ndim, interp)).initialize()

        for qty in ("EM_B_x", "density"):
            hier = hierarchy_from_sim(self.simulator, qty=qty)
            patches = hier.level(0).patches
            self.assertGreater(len(patches), 0)
            for patch in patches:
                pdata = patch.patch_datas[qty]
                self.assertEqual(pdata.dataset.shape, tuple(pdata.size))
                self.assertTrue(np.isfinite(pdata.dataset).all())

        out.mkdir(parents=True, exist_ok=True)
        hier.plot(filename=str(out / "density.png"), qty="density", levels=(0,))

        hier = hierarchy_from_sim(self.simulator, qty="particles", pop="protons")
        total = sum(
            p.patch_datas["protons_particles"].dataset.size()
            for p in hier.level(0).patches
        )
        expected = ppc * cells**ndim
        self.assertGreaterEqual(total, expected * 0.95 / cpp.mpi_size())

    def test_1d(self):
        for interp in range(1, 4):
            self.simulator = Simulator(populate_simulation(1, interp))
            self.simulator.initialize()
            dw = self.simulator.data_wrangler()

            print("\n", dw.lvl0IonDensity())
            print("\n", dw.lvl0BulkVelocity())
            print("\n", dw.lvl0PopDensity())
            print("\n", dw.lvl0PopFluxes())
            print("\n", dw.lvl0EM())

            for particles in dw.getPatchLevel(0).getParticles().values():
                for patches in particles.values():
                    for patch in patches:
                        self.assertTrue(isinstance(patch.lower, np.ndarray))
                        self.assertTrue(isinstance(patch.upper, np.ndarray))

            self.simulator.reset()
            self.simulator = None


if __name__ == "__main__":
    unittest.main()
