import unittest

from ddt import ddt, data, unpack
from pyphare.core.box import Box
from pyphare.pharesee.geometry import get_periodic_list, ghost_area_boxes
from pyphare.pharesee.geometry import (
    level_ghost_boxes,
    hierarchy_overlaps,
    touch_domain_border,
)

from pyphare.core.gridlayout import yee_element_is_primal

import numpy as np

from pyphare_tests.test_pharesee import build_hierarchy


@ddt
class GeometryTest(unittest.TestCase):
    def setup_hierarchy(self, dim, interp_order, nbr_cells, refinement_boxes, **kwargs):
        domain_size = np.asarray([1.0] * dim)
        return build_hierarchy(
            nbr_cells=nbr_cells,
            origin=np.asarray([0.0] * dim),
            interp_order=interp_order,
            domain_size=domain_size,
            cell_width=domain_size / nbr_cells,
            refinement_boxes=refinement_boxes,
            **kwargs
        )

    # used for tests without ddt hierarchy overrides
    def basic_hierarchy(self):
        dim, interp_order, nbr_cells = (1, 1, 65)
        refinement_boxes = {"L0": {"B0": [(5,), (29,)], "B1": [(32,), (55,)]}}
        return self.setup_hierarchy(dim, interp_order, nbr_cells, refinement_boxes)

    def test_overlaps(self):
        hierarchy = self.basic_hierarchy()

        expected = {
            0:
            # Middle overlap, for all quantities
            [
                {"box": Box(28, 38), "offset": (0, 0)},
                {"box": Box(28, 37), "offset": (0, 0)},
                {"box": Box(28, 37), "offset": (0, 0)},
                {"box": Box(28, 37), "offset": (0, 0)},
                {"box": Box(28, 38), "offset": (0, 0)},
                {"box": Box(28, 38), "offset": (0, 0)},
                {"box": Box(32, 33), "offset": (0, 0)},
                # left side overlap with periodicity, for all quantities
                {"box": Box(-5, 5), "offset": (0, -65)},
                # right side overlap with periodicity, for all quantities
                {"box": Box(60, 70), "offset": (65, 0)},
                {"box": Box(-5, 4), "offset": (0, -65)},
                {"box": Box(60, 69), "offset": (65, 0)},
                {"box": Box(-5, 4), "offset": (0, -65)},
                {"box": Box(60, 69), "offset": (65, 0)},
                {"box": Box(-5, 4), "offset": (0, -65)},
                {"box": Box(60, 69), "offset": (65, 0)},
                {"box": Box(-5, 5), "offset": (0, -65)},
                {"box": Box(60, 70), "offset": (65, 0)},
                {"box": Box(-5, 5), "offset": (0, -65)},
                {"box": Box(60, 70), "offset": (65, 0)},
                {"box": Box(-1, 0), "offset": (0, -65)},
                {"box": Box(64, 65), "offset": (65, 0)},
            ],
            1: [  # level 1
                {"box": Box(59, 65), "offset": (0, 0)},
                {"box": Box(59, 64), "offset": (0, 0)},
                {"box": Box(59, 64), "offset": (0, 0)},
                {"box": Box(59, 64), "offset": (0, 0)},
                {"box": Box(59, 65), "offset": (0, 0)},
                {"box": Box(59, 65), "offset": (0, 0)},
            ],
        }

        overlaps = hierarchy_overlaps(hierarchy)

        for ilvl, lvl in enumerate(hierarchy.patch_levels):
            self.assertEqual(len(expected[ilvl]), len(overlaps[ilvl]))

            for exp, actual in zip(expected[ilvl], overlaps[ilvl]):
                act_box = actual["box"]
                act_offset = actual["offset"]
                exp_box = exp["box"]
                exp_offset = exp["offset"]

                self.assertEqual(act_box, exp_box)
                self.assertEqual(act_offset, exp_offset)

    def test_touch_border(self):
        hierarchy = self.basic_hierarchy()

        self.assertFalse(
            touch_domain_border(Box(10, 20), hierarchy.domain_box, "upper")
        )
        self.assertFalse(
            touch_domain_border(Box(10, 20), hierarchy.domain_box, "lower")
        )
        self.assertTrue(touch_domain_border(Box(0, 20), hierarchy.domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box(-5, 20), hierarchy.domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box(-5, 70), hierarchy.domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box(-5, 70), hierarchy.domain_box, "upper"))
        self.assertTrue(touch_domain_border(Box(40, 70), hierarchy.domain_box, "upper"))
        self.assertTrue(touch_domain_border(Box(40, 64), hierarchy.domain_box, "upper"))

    def test_particle_ghost_area_boxes(self):
        hierarchy = self.basic_hierarchy()

        expected = {
            0: [Box(-1, -1), Box(33, 33), Box(32, 32), Box(65, 65)],
            1: [Box(9, 9), Box(60, 60), Box(63, 63), Box(112, 112)],
        }

        gaboxes = ghost_area_boxes(hierarchy, "particles")
        particles = "particles"

        # same number of levels
        self.assertEqual(len(expected), len(gaboxes))

        for ilvl, lvl in enumerate(hierarchy.patch_levels):
            qtyNbr = len(gaboxes[ilvl].keys())
            self.assertEqual(qtyNbr, 1)

            key = list(gaboxes[ilvl].keys())[0]

            level_ghost_area_boxes = sum(  # aggregate to single list
                [actual["boxes"] for actual in gaboxes[ilvl][key]], []
            )

            self.assertEqual(expected[ilvl], level_ghost_area_boxes)

    def test_level_ghost_boxes(self):
        hierarchy = self.basic_hierarchy()

        expected = {
            1: [
                {
                    "boxes": [
                        Box(9, 9),
                    ]
                },
                {
                    "boxes": [
                        Box(60, 60),
                    ]
                },
                {
                    "boxes": [
                        Box(63, 63),
                    ]
                },
                {
                    "boxes": [
                        Box(112, 112),
                    ]
                },
            ]
        }

        lvl_gaboxes = level_ghost_boxes(hierarchy, "particles")
        for ilvl in range(1, len(hierarchy.patch_levels)):
            qtyNbr = len(lvl_gaboxes[ilvl].keys())
            self.assertEqual(qtyNbr, 1)

            key = list(lvl_gaboxes[ilvl].keys())[0]

            for actual, exp in zip(lvl_gaboxes[ilvl][key], expected[ilvl]):
                act_boxes = actual["boxes"]
                exp_boxes = exp["boxes"]

                for act_box, exp_box in zip(act_boxes, exp_boxes):
                    self.assertEqual(act_box, exp_box)

    def test_level_ghost_boxes_do_not_overlap_patch_interiors(self):
        hierarchy = self.basic_hierarchy()

        lvl_gaboxes = level_ghost_boxes(hierarchy, "particles")

        for ilvl in range(1, len(hierarchy.patch_levels)):
            qtyNbr = len(lvl_gaboxes[ilvl].keys())
            self.assertEqual(qtyNbr, 1)

            key = list(lvl_gaboxes[ilvl].keys())[0]

            for pdatainfo in lvl_gaboxes[ilvl][key]:
                for box in pdatainfo["boxes"]:
                    for patch in hierarchy.patch_levels[ilvl].patches:
                        self.assertIsNone(patch.box * box)

    @data(
        (
            {
                "L0": [
                    Box(0, 4),
                    Box(15, 19),
                ]
            },
            {
                1: [
                    {"box": Box(-5, 5), "offset": (0, -40)},
                    {"box": Box(35, 45), "offset": (40, 0)},
                ]
            },
        ),
        (
            {
                "L0": [
                    Box(1, 5),
                    Box(14, 18),
                ]
            },
            {
                1: [
                    {"box": Box(-3, 3), "offset": (0, -40)},
                    {"box": Box(37, 43), "offset": (40, 0)},
                ]
            },
        ),
    )
    @unpack
    def test_periodic_overlaps(self, refinement_boxes, expected):
        dim, interp_order, nbr_cells = (1, 1, 20)
        hierarchy = self.setup_hierarchy(
            dim, interp_order, nbr_cells, refinement_boxes, quantities="Bx"
        )

        overlaps = hierarchy_overlaps(hierarchy)
        for ilvl, lvl in enumerate(hierarchy.patch_levels):
            if ilvl not in expected:
                continue
            self.assertEqual(len(expected[ilvl]), len(overlaps[ilvl]))
            for exp, actual in zip(expected[ilvl], overlaps[ilvl]):
                self.assertEqual(actual["box"], exp["box"])
                self.assertEqual(actual["offset"], exp["offset"])

    @data(
        (
            {
                "L0": [
                    Box(0, 4),
                    Box(15, 19),
                ]
            },
            {
                1: [
                    Box(-10, -1),
                    Box(0, 9),
                    Box(30, 39),
                    Box(40, 49),
                ]
            },
        ),
        (
            {
                "L0": [
                    Box(1, 5),
                    Box(14, 18),
                ]
            },
            {
                1: [
                    Box(-12, -3),
                    Box(2, 11),
                    Box(28, 37),
                    Box(42, 51),
                ]
            },
        ),
    )
    @unpack
    def test_periodic_list(self, refinement_boxes, expected):
        dim, interp_order, nbr_cells = (1, 1, 20)
        hierarchy = self.setup_hierarchy(
            dim, interp_order, nbr_cells, refinement_boxes, quantities="Bx"
        )

        for ilvl in range(1, len(hierarchy.patch_levels)):
            refined_domain_box = hierarchy.refined_domain_box(ilvl)

            n_ghosts = (
                hierarchy.patch_levels[ilvl].patches[0].patch_datas["Bx"].ghosts_nbr
            )
            patches = get_periodic_list(
                hierarchy.patch_levels[ilvl].patches, refined_domain_box, n_ghosts
            )

            periodic_boxes = [patch.box for patch in patches]

            for ref_box_i, ref_box in enumerate(periodic_boxes):
                for cmp_box in periodic_boxes[ref_box_i + 1 :]:
                    self.assertTrue(ref_box * cmp_box is None)

            self.assertEqual(expected[ilvl], periodic_boxes)

    @data(
        (
            {
                "L0": [
                    Box(0, 4),
                    Box(15, 19),
                ]
            },
            {
                1: [
                    Box(10, 14),
                    Box(25, 29),
                ]
            },
        ),
        (
            {
                "L0": [
                    Box(1, 5),
                    Box(14, 18),
                ]
            },
            {
                1: [
                    Box(-2, 1),
                    Box(12, 16),
                    Box(23, 27),
                    Box(38, 41),
                ]
            },
        ),
        (
            {
                "L0": [
                    Box(5, 9),
                    Box(10, 14),
                ]
            },
            {
                1: [
                    Box(5, 9),
                    Box(30, 34),
                ]
            },
        ),
    )
    @unpack
    def test_level_ghostboxes(self, refinement_boxes, expected):
        dim, interp_order, nbr_cells = (1, 1, 20)
        hierarchy = self.setup_hierarchy(
            dim, interp_order, nbr_cells, refinement_boxes, quantities="Bx"
        )

        lvl_gaboxes = level_ghost_boxes(hierarchy, "Bx")
        for ilvl in range(1, len(hierarchy.patch_levels)):
            qtyNbr = len(lvl_gaboxes[ilvl].keys())
            self.assertEqual(qtyNbr, 1)

            key = list(lvl_gaboxes[ilvl].keys())[0]

            ghost_area_boxes = sum(  # aggregate to single list
                [actual["boxes"] for actual in lvl_gaboxes[ilvl][key]], []
            )

            self.assertEqual(expected[ilvl], ghost_area_boxes)

    def test_field_data_select(self):
        dim, interp_order, nbr_cells = (1, 1, 20)

        for qty in ["Bx", "By"]:
            hierarchy = self.setup_hierarchy(
                dim, interp_order, nbr_cells, {}, quantities=qty
            )

            pdata = hierarchy.level(0).patches[0].patch_datas[qty]
            lower_gb, upper_gb = sorted(
                pdata.ghost_box - pdata.box, key=lambda box: box.lower.all()
            )
            lower_ds, upper_ds = pdata[lower_gb], pdata[upper_gb]

            qty_is_primal = yee_element_is_primal(qty)
            assert lower_ds.shape[0] == pdata.ghosts_nbr[0] + qty_is_primal
            assert upper_ds.shape[0] == pdata.ghosts_nbr[0] + qty_is_primal

            np.testing.assert_array_equal(lower_ds, pdata.dataset[: lower_ds.shape[0]])
            np.testing.assert_array_equal(upper_ds, pdata.dataset[-upper_ds.shape[0] :])


if __name__ == "__main__":
    unittest.main()
