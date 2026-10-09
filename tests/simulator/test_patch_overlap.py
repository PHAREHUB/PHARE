#!/usr/bin/env python3
"""
Same-level patches may overlap: SAMRAI grows boxes smaller than the minimum patch size
without looking at the other new boxes. The configuration below produces such an overlap
on level 2 at initialization. Cells covered by several patches belong to the patch with
the smallest GlobalId, which alone holds their particles, while the fields on the shared
nodes stay identical on all patches.
"""

import unittest
import numpy as np

import pyphare.pharein as ph
from pyphare.core.box import Box
from pyphare.core.gridlayout import yee_centering
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest

ph.NO_GUI()

ndim = 2
interp = 2
cells = (64, 64)
dl = (0.4, 0.4)
time_step = 0.001
time_step_nbr = 3
ppc = 20

# Gaussian bumps of Bx tagged for refinement, positions chosen so that level 2 has a
# two-cell wide overlap between a large patch and a grown one
bump_centers = [
    (12.98, 19.72),
    (7.33, 19.69),
    (9.91, 11.62),
    (17.83, 11.41),
]
bump_widths = [1.37, 0.64, 1.65, 1.35]


def config():
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_nbr_levels=3,
        interp_order=interp,
        nesting_buffer=1,
        tagging_threshold=0.1,
        hyper_resistivity=0.002,
        resistivity=0.001,
        strict=True,
    )

    def bx(x, y):
        bumps = 0.0 * x
        for (cx, cy), w in zip(bump_centers, bump_widths):
            bumps = bumps + np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / w**2)
        return 1.0 + bumps

    def zero(x, y):
        return 0.0 * x

    def one(x, y):
        return 1.0 + 0.0 * x

    def vth(x, y):
        return 0.3 + 0.0 * x

    ph.MaxwellianFluidModel(
        bx=bx,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": one,
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    return sim


def global_id(patch):
    rank, local = patch.patchID.split("#")
    return int(rank), int(local)


def box_of(patch):
    return Box(patch.lower, patch.upper)


def overlapping_pairs(patches):
    pairs = []
    for i, p0 in enumerate(patches):
        for p1 in patches[i + 1 :]:
            if box_of(p0) * box_of(p1) is not None:
                pairs.append(tuple(sorted((p0, p1), key=global_id)))
    return pairs


def node_range(patch, centering, direction):
    # AMR index range of the nodes of `patch` interior to its domain
    lower, upper = patch.lower[direction], patch.upper[direction]
    return lower, upper + (1 if centering == "primal" else 0)


def shared_nodes(owner, other, qty):
    """values of `qty` on the nodes interior to both patch domains"""
    slices = ([], [])
    for d, direction in enumerate("xy"[:ndim]):
        centering = yee_centering[direction][qty]
        lo0, hi0 = node_range(owner, centering, d)
        lo1, hi1 = node_range(other, centering, d)
        lo, hi = max(lo0, lo1), min(hi0, hi1)
        for patch, s in zip((owner, other), slices):
            start = lo - patch.lower[d] + patch.nGhosts
            s.append(slice(start, start + hi - lo + 1))

    def field(patch):
        shape = [
            patch.upper[d]
            - patch.lower[d]
            + 1
            + 2 * patch.nGhosts
            + (1 if yee_centering[direction][qty] == "primal" else 0)
            for d, direction in enumerate("xy"[:ndim])
        ]
        return np.asarray(patch.data).reshape(shape)

    return field(owner)[tuple(slices[0])], field(other)[tuple(slices[1])]


class PatchOverlapTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        super().tearDown()

    def level_patches(self, ilvl, getter, *args):
        level = self.simulator.data_wrangler().getPatchLevel(ilvl)
        return getattr(level, getter)(*args)

    def check_particles_are_on_owner_only(self, ilvl, pairs):
        particles = self.level_patches(ilvl, "getParticles", "protons")
        particles = particles["protons"]["domain"]
        by_id = {p.patchID: p for p in particles}

        def nbr_in(patch, box):
            if patch.patchID not in by_id:  # no domain particles at all
                return 0
            icells = np.asarray(by_id[patch.patchID].data.iCell).reshape(-1, ndim)
            inside = np.all((icells >= box.lower) & (icells <= box.upper), axis=1)
            return np.count_nonzero(inside)

        for owner, other in pairs:
            overlap = box_of(owner) * box_of(other)
            self.assertGreater(nbr_in(owner, overlap), 0)
            self.assertEqual(nbr_in(other, overlap), 0)

    def check_shared_nodes_are_equal(self, ilvl, pairs):
        getters = {
            "Bx": "getBx", "By": "getBy", "Bz": "getBz",
            "Ex": "getEx", "Ey": "getEy", "Ez": "getEz",
            "rho": "getDensity",
        }  # fmt: skip
        for qty, getter in getters.items():
            by_id = {p.patchID: p for p in self.level_patches(ilvl, getter)}
            for owner, other in pairs:
                owner_values, other_values = shared_nodes(
                    by_id[owner.patchID], by_id[other.patchID], qty
                )
                self.assertGreater(owner_values.size, 0)
                np.testing.assert_array_equal(owner_values, other_values, err_msg=qty)

    def check_density_is_not_duplicated(self, ilvl, pairs):
        by_id = {p.patchID: p for p in self.level_patches(ilvl, "getDensity")}
        for owner, other in pairs:
            values, _ = shared_nodes(by_id[owner.patchID], by_id[other.patchID], "rho")
            # the initial density is 1, it doubles if both patches hold the particles
            self.assertLess(abs(values.mean() - 1.0), 0.1)

    def test_overlapping_patches_partition_particles(self):
        self.simulator = Simulator(config()).initialize()

        ilvl = 2
        pairs = overlapping_pairs(self.level_patches(ilvl, "getBx"))
        self.assertGreater(len(pairs), 0, "the configuration no longer has an overlap")

        for step in range(time_step_nbr + 1):
            if step > 0:
                self.simulator.advance()
            pairs = overlapping_pairs(self.level_patches(ilvl, "getBx"))
            self.check_particles_are_on_owner_only(ilvl, pairs)
            self.check_shared_nodes_are_equal(ilvl, pairs)
            self.check_density_is_not_duplicated(ilvl, pairs)


if __name__ == "__main__":
    unittest.main()
