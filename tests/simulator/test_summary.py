#!/usr/bin/env python3
#
# Simulator.summary(): per level values reduced over all ranks onto rank 0.
#  Checked against what the layout fixes whatever the load balancer does:
#  a periodic 20x20 L0 and one refinement box of 10x10 coarse cells.

import unittest

from pyphare import cpp
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest, populate_simulation

cells = [20, 20]
refinement_box = [[3, 3], [12, 12]]  # coarse cells, L1 covers 10x10 of them at ratio 2


class SummaryTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None
        self.diag_dirs = []
        self.success = False

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        super().tearDown()

    def test_summary(self):
        simulation = populate_simulation(
            2, 1, cells=cells, refinement_boxes={"L0": {"B0": refinement_box}}
        )
        self.simulator = Simulator(simulation).initialize()
        summary = self.simulator.summary()  # collective

        if cpp.mpi_rank() > 0:
            self.assertIsNone(summary)
            self.success = True
            return

        L0, L1 = summary["levels"]
        size = cpp.mpi_size()

        # global values: a rank 0 only reduction loses the other ranks' cells
        self.assertEqual(L0["cells"], 20 * 20)
        self.assertEqual(L1["cells"], (10 * 2) ** 2)
        for lvl in (L0, L1):
            self.assertGreaterEqual(lvl["patches"], 1)
            self.assertGreater(lvl["particles"], 0)
            self.assertGreaterEqual(lvl["cell_imbalance"], 1)
            self.assertLessEqual(lvl["cell_imbalance"], size)
            self.assertGreaterEqual(lvl["particle_imbalance"], 1)
            self.assertLessEqual(lvl["remote_ghosts"], lvl["patch_ghosts"])
            self.assertLessEqual(lvl["remote_ranks_max"], size - 1)
            self.assertLessEqual(lvl["remote_ranks_avg"], lvl["remote_ranks_max"])

        # periodic L0: every ghost is covered by a patch or a periodic image
        self.assertEqual(L0["level_ghosts"], 0)
        self.assertGreater(L0["patch_ghosts"], 0)
        # L1 is inside the domain: its outer ghosts are level ghosts
        self.assertGreater(L1["level_ghosts"], 0)

        if size == 1:
            for lvl in (L0, L1):
                self.assertEqual(lvl["remote_ghosts"], 0)
                self.assertEqual(lvl["remote_ranks_max"], 0)
        elif L0["cell_imbalance"] < size:  # every rank owns L0 cells
            # periodic: some L0 patches of different ranks touch
            self.assertGreater(L0["remote_ghosts"], 0)
            self.assertGreaterEqual(L0["remote_ranks_max"], 1)

        self.success = True


if __name__ == "__main__":
    unittest.main()
