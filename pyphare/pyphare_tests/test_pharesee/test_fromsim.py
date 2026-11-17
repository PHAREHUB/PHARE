"""
The live data reader (pharesee.hierarchy.fromsim) must take the field ghost count the C++ layout
reports (PatchData.nGhosts, python3/patch_data.hpp), not the hybrid table of the interpolation
order: the fully kinetic layout has 4, 6, 6 field ghosts where the hybrid has 2, 4, 4 (P4). The
same for particles: 3, 3, 4 ghost cells where the hybrid has 1, 2, 2. Fake simulator and data wrangler objects stand for the C++ ones.
"""

import unittest

import numpy as np

from pyphare.core.box import Box
from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim, make_layout_for
from pyphare.pharesee.hierarchy.hierarchy_utils import field_qties
from pyphare.pharesee.hierarchy.patchdata import FieldData

CELLS = 10
GHOSTS = 4  # fully kinetic, order 1; the hybrid table says 2


class FakePatchData:
    """a PatchData of the data wrangler: Ex (dual along x) of one 1D patch, with ghosts"""

    def __init__(self, ghosts):
        self.lower = np.array([0])
        self.upper = np.array([CELLS - 1])
        self.origin = np.array([0.0])
        self.nGhosts = ghosts
        self.data = np.arange(CELLS + 2 * ghosts, dtype=float)


class FakeParticlePatchData:
    """a particle PatchData of the data wrangler: two domain particles of one 1D patch, flat
    arrays as ContiguousParticles packs them"""

    class Data:
        iCell = [2, 7]
        delta = [0.25, 0.5]
        weight = [0.1, 0.2]
        v = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

    def __init__(self, ghosts):
        self.lower = np.array([0])
        self.upper = np.array([CELLS - 1])
        self.origin = " 0.000000"
        self.nGhosts = ghosts
        self.data = self.Data()


class FakePatchLevel:
    def __init__(self, patch):
        self.patch = patch

    def __getattr__(self, getter):  # getEx, getBx, getDensity, ...
        return lambda *args: [self.patch]

    def getParticles(self, pop):
        return {pop: {"domain": [self.patch]}}


class FakeWrangler:
    def __init__(self, patch):
        self.patch = patch

    def getNumberOfLevels(self):
        return 1

    def getPatchLevel(self, ilvl):
        return FakePatchLevel(self.patch)


class FakeSimulator:
    def __init__(self, patch):
        self.patch = patch

    def data_wrangler(self):
        return FakeWrangler(self.patch)

    def cell_width(self):
        return np.array([0.1])

    def domain_box(self):
        return np.array([CELLS - 1])

    def interp_order(self):
        return 1

    def currentTime(self):
        return 0.0


class FromSimTest(unittest.TestCase):
    def check_physical(self, pdata, ghosts):
        physical = pdata[pdata.box]
        self.assertEqual(len(physical), CELLS)
        self.assertEqual(physical[0], ghosts)  # data[i] = i: the first physical node

    def test_layout_takes_the_reported_field_ghosts(self):
        patch = FakePatchData(GHOSTS)
        simulator = FakeSimulator(patch)
        layout = make_layout_for(simulator, patch, "EM_E_x", simulator.cell_width())
        pdata = FieldData(layout, field_qties["EM_E_x"], patch.data)
        self.assertEqual(list(pdata.ghosts_nbr), [GHOSTS])
        self.check_physical(pdata, GHOSTS)

    def test_hybrid_ghosts_unchanged(self):
        patch = FakePatchData(2)  # what the hybrid layout reports at order 1
        simulator = FakeSimulator(patch)
        layout = make_layout_for(simulator, patch, "EM_E_x", simulator.cell_width())
        self.assertEqual(list(layout.ghosts_nbr), [2])

    def test_2d_field_is_shaped_by_the_reported_ghosts(self):
        """Bz (dual, dual) of a 10 x 6 patch: the flat data takes the ghost box shape"""
        patch = FakePatchData(GHOSTS)
        patch.lower, patch.upper = np.array([0, 0]), np.array([CELLS - 1, 5])
        patch.origin = " 0.000000, 0.000000"  # as C++ sends it
        shape = (CELLS + 2 * GHOSTS, 6 + 2 * GHOSTS)
        patch.data = list(np.arange(shape[0] * shape[1], dtype=float))
        simulator = FakeSimulator(patch)
        simulator.cell_width = lambda: [0.1, 0.2]
        simulator.domain_box = lambda: np.array([CELLS - 1, 5])
        hier = hierarchy_from_sim(simulator, "EM_B_z")
        pdata = hier.level(0).patches[0].patch_datas["EM_B_z"]
        self.assertEqual(pdata.dataset.shape, shape)
        physical = pdata[pdata.box]
        self.assertEqual(physical.shape, (CELLS, 6))
        self.assertEqual(physical[0, 0], GHOSTS * shape[1] + GHOSTS)

    def test_particles_take_the_reported_ghosts(self):
        PARTICLE_GHOSTS = 3  # fully kinetic, order 1; the hybrid table says 1
        hier = hierarchy_from_sim(
            FakeSimulator(FakeParticlePatchData(PARTICLE_GHOSTS)), "particles", pop="protons"
        )
        pdata = hier.level(0).patches[0].patch_datas["protons_particles"]
        self.assertEqual(list(pdata.ghosts_nbr), [PARTICLE_GHOSTS])
        self.assertEqual(pdata.ghost_box, Box([-PARTICLE_GHOSTS], [CELLS - 1 + PARTICLE_GHOSTS]))
        particles = pdata.dataset
        np.testing.assert_array_equal(particles.iCells, [[2], [7]])
        np.testing.assert_allclose(particles.x, [0.225, 0.75])  # (iCell + delta) dx
        np.testing.assert_array_equal(particles.v, [[1, 2, 3], [4, 5, 6]])

    def test_hierarchy_from_sim_reads_the_physical_nodes(self):
        hier = hierarchy_from_sim(FakeSimulator(FakePatchData(GHOSTS)), "EM_E_x")
        patches = hier.level(0).patches
        self.assertEqual(len(patches), 1)
        self.assertEqual(patches[0].box, Box([0], [CELLS - 1]))
        self.check_physical(patches[0].patch_datas["EM_E_x"], GHOSTS)


if __name__ == "__main__":
    unittest.main()
