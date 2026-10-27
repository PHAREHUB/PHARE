"""
hierarchy.finest(): on the domain, the finest level wins where it exists,
the coarser level elsewhere, for every Yee centering
"""

import unittest

import numpy as np
from ddt import data, ddt, unpack

from pyphare.pharesee.hierarchy.patchdata import local_domain_coords
from pyphare_tests.test_pharesee import build_hierarchy

quantities = ["Bx", "By", "Bz", "Ex", "Ey", "Ez"]
fine_offset = 100.0  # analytic fields are in [-1, 1]: tells fine from coarse values

cases = [
    (1, [40], {"L0": {"B0": [(10,), (19,)]}}),
    (2, [20, 20], {"L0": {"B0": [(5, 5), (9, 9)]}}),
    (2, [20, 20], {"L0": {"B0": [(2, 3), (6, 8)], "B1": [(12, 12), (16, 17)]}}),
]


def hierarchy(ndim, nbr_cells, refinement_boxes):
    domain_size = np.ones(ndim)
    return build_hierarchy(
        nbr_cells=np.asarray(nbr_cells),
        origin=np.zeros(ndim),
        interp_order=1,
        domain_size=domain_size,
        cell_width=domain_size / np.asarray(nbr_cells),
        refinement_boxes=refinement_boxes,
        quantities=list(quantities),
    )


def distance_to(coords, grid_coords, cell_width):
    """per grid point, the distance to the box spanned by coords, in cell widths"""
    mesh = np.meshgrid(*grid_coords, indexing="ij")
    dist = np.zeros(mesh[0].shape)
    for m, c, dl in zip(mesh, coords, cell_width):
        dist = np.maximum(dist, np.maximum(c.min() - m, m - c.max()).clip(0) / dl)
    return dist


def indices_of(coords, grid_coords):
    """indices of coords in grid_coords, which must contain them all"""
    idx = [np.searchsorted(g, c - 1e-9) for c, g in zip(coords, grid_coords)]
    for i, c, g in zip(idx, coords, grid_coords):
        np.testing.assert_allclose(g[i], c, atol=1e-12)
    return np.ix_(*idx)


@ddt
class FinestTest(unittest.TestCase):
    @data(*cases)
    @unpack
    def test_finest_level_wins(self, ndim, nbr_cells, refinement_boxes):
        hier = hierarchy(ndim, nbr_cells, refinement_boxes)
        time, fine = 0.0, hier.finest_level()
        coarse_dl = hier.level(0).cell_width
        self.assertEqual(fine, 1)
        for patch in hier.level(fine).patches:
            for qty in quantities:
                patch[qty].dataset[:] += fine_offset

        for qty in quantities:
            grid = hier.finest(time, qty)
            values = np.asarray(grid[grid.box])
            self.assertFalse(np.isnan(values).any(), qty)

            # every fine patch is reproduced on its own coordinates
            grid_coords = local_domain_coords(grid)
            near_fine = np.zeros(values.shape, dtype=bool)
            for patch in hier.level(fine).patches:
                pd = patch[qty]
                coords = local_domain_coords(pd)
                select = indices_of(coords, grid_coords)
                np.testing.assert_array_equal(values[select], pd[pd.box], err_msg=qty)
                # nearest interpolation may take a fine value up to half a coarse cell out
                near_fine |= distance_to(coords, grid_coords, coarse_dl) <= 0.5 + 1e-9

            # further out, values come from the coarse level
            self.assertTrue((np.abs(values[~near_fine]) <= 1).all(), qty)


if __name__ == "__main__":
    unittest.main()
