#!/usr/bin/env python3
"""
Coverage for hierarchy_from(simulator=..., qty=...).
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from, PatchHierarchy
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

diag_outputs = "phare_outputs/test_fromsim"
time_step = 0.001
cells = (20, 20)
dl = (0.2, 0.2)
nghosts = 2  # interp_order 1


class FromSimTest(unittest.TestCase):
    def diag_dir(self):
        return f"{diag_outputs}/{self._testMethodName}"

    def setUp(self):
        ph.global_vars.sim = None

        sim = ph.Simulation(
            time_step=time_step,
            time_step_nbr=1,
            cells=cells,
            dl=dl,
            refinement_boxes={"L0": {"B0": [(5, 5), (14, 14)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
            diag_options={
                "format": "phareh5",
                "options": {"dir": self.diag_dir(), "mode": "overwrite"},
            },
        )

        def density(x, y):
            return 1.0

        def bx(x, y):
            return 1.0

        def zero(x, y):
            return 0.0

        def vth(x, y):
            return 0.3

        ph.MaxwellianFluidModel(
            bx=bx,
            by=zero,
            bz=zero,
            protons={
                "charge": 1,
                "density": density,
                **{f"vbulk{c}": zero for c in "xyz"},
                **{f"vth{c}": vth for c in "xyz"},
                "nbr_part_per_cell": 10,
                "init": {"seed": 13337},
            },
        )
        ph.ElectronModel(closure="isothermal", Te=0.12)
        ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([0.0]))

        self.sim = sim
        self.simulator = Simulator(sim)
        self.simulator.initialize()

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_field_quantities_keep_their_yee_shape(self):
        """A flat wrangler buffer must come back shaped like the layout, as an
        HDF5 read would, and per component: Bx is primal in x and dual in y."""
        for qty, (primal_x, primal_y) in {
            "EM_B_x": (True, False),
            "EM_B_y": (False, True),
            "EM_B_z": (False, False),
            "EM_E_x": (False, True),
            "EM_E_y": (True, False),
            "EM_E_z": (True, True),
            "density": (True, True),
        }.items():
            hier = hierarchy_from(simulator=self.simulator, qty=qty)
            levels = hier.levels()
            self.assertGreaterEqual(len(levels), 2, f"{qty}: {len(levels)} level(s)")

            for ilvl, level in levels.items():
                self.assertTrue(level.patches, f"{qty} lvl={ilvl}: no patch")
                for patch in level.patches:
                    ghost_cells = patch.box.shape + 2 * nghosts
                    self.assertEqual(
                        tuple(patch.patch_datas[qty].dataset.shape),
                        (
                            ghost_cells[0] + int(primal_x),
                            ghost_cells[1] + int(primal_y),
                        ),
                        f"{qty} lvl={ilvl}",
                    )

    def test_hierarchy_carries_the_simulator_time(self):
        """not the PatchHierarchy default: the kwarg is `times`, not `time`"""
        hier = hierarchy_from(simulator=self.simulator, qty="EM_B_x")
        self.assertIsInstance(hier, PatchHierarchy)
        self.assertEqual(len(hier.times()), 1)
        self.assertAlmostEqual(float(hier.times()[0]), self.simulator.currentTime())

    def test_in_memory_fields_match_the_dumped_ones(self):
        """the strongest check available: same run, same patches, two read paths"""
        in_memory = {
            qty: hierarchy_from(simulator=self.simulator, qty=qty)
            for qty in ("EM_B_x", "EM_B_y", "EM_B_z")
        }
        self.simulator.dump(self.simulator.currentTime(), self.sim.time_step)
        dumped = Run(self.diag_dir()).GetB(0.0, all_primal=False)

        compared = 0
        for qty, hier in in_memory.items():
            component = qty[-1]
            for ilvl, level in hier.levels().items():
                for patch in level.patches:
                    match = [
                        p
                        for p in dumped.levels(0.0)[ilvl].patches
                        if p.box == patch.box
                    ]
                    self.assertTrue(match, f"{qty} lvl={ilvl}: no {patch.box} in dump")
                    np.testing.assert_array_equal(
                        patch.patch_datas[qty].dataset[:],
                        match[0].patch_datas[f"B{component}"].dataset[:],
                        err_msg=f"{qty} lvl={ilvl} {patch.box}",
                    )
                    compared += 1
        self.assertGreater(compared, 0)

    def test_particles_are_readable_in_memory(self):
        """the particles branch: flat buffers need reshaping and a per particle dl"""
        hier = hierarchy_from(simulator=self.simulator, qty="particles", pop="protons")
        levels = hier.levels()
        self.assertGreaterEqual(len(levels), 2)

        total = 0
        for ilvl, level in levels.items():
            self.assertTrue(level.patches, f"lvl={ilvl}: no patch")
            for patch in level.patches:
                particles = patch.patch_datas["protons_particles"].dataset
                self.assertEqual(particles.iCells.shape[1], len(cells))
                self.assertEqual(particles.deltas.shape, particles.iCells.shape)
                self.assertEqual(particles.dl.shape, particles.iCells.shape)
                self.assertEqual(particles.v.shape[1], 3)
                total += particles.size()
        self.assertGreater(total, 0)


if __name__ == "__main__":
    unittest.main()
