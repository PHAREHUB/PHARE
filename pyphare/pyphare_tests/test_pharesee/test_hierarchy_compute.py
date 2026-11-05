"""
compute_hier_from: callbacks take a patch and return the FieldDatas of the new patch
"""

import unittest

import numpy as np

from pyphare.pharesee.hierarchy import hierarchy_compute as hc
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from
from pyphare_tests.test_pharesee import build_hierarchy


def hierarchy():
    return build_hierarchy(
        nbr_cells=np.asarray([20, 20]),
        origin=np.zeros(2),
        interp_order=1,
        domain_size=np.ones(2),
        cell_width=np.ones(2) / 20,
        refinement_boxes={"L0": {"B0": [(5, 5), (9, 9)]}},
        quantities=["Bx", "By"],
    )


def all_patches(hier):
    return [p for lvl in hier.levels().values() for p in lvl.patches]


class HierarchyComputeTest(unittest.TestCase):
    def test_callback_gets_patch_returns_field_datas(self):
        hier = hierarchy()
        seen = []

        def twice(patch, **kwargs):
            seen.append(kwargs["factor"])
            pd = patch["Bx"]
            return (pd.copy_as(pd[:] * kwargs["factor"], name="twice"),)

        twice_hier = compute_hier_from(twice, hier, factor=2)

        # called once per patch, with the kwargs given to compute_hier_from
        self.assertEqual(len(seen), len(all_patches(hier)))
        self.assertEqual(twice_hier.levelNbr(), hier.levelNbr())
        for ref, new in zip(all_patches(hier), all_patches(twice_hier)):
            self.assertEqual(list(new.patch_datas), ["twice"])
            self.assertEqual(new.box, ref.box)
            pd, ref_pd = new["twice"], ref["Bx"]
            self.assertEqual(pd.ghosts_nbr.tolist(), ref_pd.ghosts_nbr.tolist())
            self.assertEqual(pd.centerings, ref_pd.centerings)
            np.testing.assert_array_equal(pd[:], 2 * ref_pd[:])

    def test_drop_ghosts_keeps_dataset_without_ghosts(self):
        """no ghosts, nothing to drop: the dataset is passed on unread, as is"""
        no_ghosts = compute_hier_from(hc.drop_ghosts, hierarchy())
        again = compute_hier_from(hc.drop_ghosts, no_ghosts)
        for patch, patch_again in zip(all_patches(no_ghosts), all_patches(again)):
            for name, pd in patch.patch_datas.items():
                self.assertEqual(pd.ghosts_nbr.tolist(), [0, 0])
                self.assertIs(patch_again[name].dataset, pd.dataset)

    def test_rename(self):
        renamed = hc.rename(hierarchy(), ["a", "b"])
        for patch in all_patches(renamed):
            self.assertEqual(list(patch.patch_datas), ["a", "b"])

        with self.assertRaises(ValueError) as error:
            hc.rename(hierarchy(), ["a"])  # two quantities, one name
        self.assertIn("count mismatch", str(error.exception))


if __name__ == "__main__":
    unittest.main()
