import unittest
import numpy as np
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box2D
from pyphare.pharesee.geometry import (
    level_ghost_boxes,
    hierarchy_overlaps,
    touch_domain_border,
    ghost_area_boxes,
    get_periodic_list,
)

from pyphare_tests.test_pharesee import build_hierarchy


class AGeometryTest(unittest.TestCase):
    def ddt_test_id(self):
        """
        ddt test functions end with "_${ID}"
        """
        return self._testMethodName.split("_")[-1]

    def setup_hierarchy(self, dim, interp_order, nbr_cells, refinement_boxes, **kwargs):
        domain_size = np.asarray([1.0] * dim)
        return build_hierarchy(
            nbr_cells=nbr_cells,
            origin=np.asarray([0.0] * dim),
            interp_order=interp_order,
            domain_size=domain_size,
            cell_width=domain_size / nbr_cells,
            refinement_boxes=refinement_boxes,
            **kwargs,
        )

    # used for tests without ddt hierarchy overrides
    def basic_hierarchy(self):
        dim, interp_order, nbr_cells = (2, 1, [20] * 2)
        refinement_boxes = {"L0": {"B0": Box2D(5, 9), "B1": Box2D(14, 19)}}
        return self.setup_hierarchy(dim, interp_order, nbr_cells, refinement_boxes)


@ddt
class GeometryTest(AGeometryTest):
    @data(
        (
            {
                # no level 1
            },
            {
                0: [
                    {"box": Box([5, 5], [15, 14]), "offset": ([0, 0], [0, 0])},
                    {"box": Box([5, -5], [15, 4]), "offset": ([0, 0], [0, -20])},
                    {"box": Box([5, 15], [15, 24]), "offset": ([0, 20], [0, 0])},
                    {"box": Box([-5, 5], [5, 14]), "offset": ([0, 0], [-20, 0])},
                    {"box": Box([15, 5], [25, 14]), "offset": ([20, 0], [0, 0])},
                    {"box": Box([-5, -5], [5, 4]), "offset": ([0, 0], [-20, -20])},
                    {"box": Box([15, 15], [25, 24]), "offset": ([20, 20], [0, 0])},
                    {"box": Box([-5, 5], [15, 14]), "offset": ([0, 0], [0, 0])},
                    {"box": Box([-5, -5], [15, 4]), "offset": ([0, 0], [0, -20])},
                    {"box": Box([-5, 15], [15, 24]), "offset": ([0, 20], [0, 0])},
                    {"box": Box([5, -5], [15, 14]), "offset": ([0, 0], [0, 0])},
                    {"box": Box([-5, -5], [5, 14]), "offset": ([0, 0], [-20, 0])},
                    {"box": Box([15, -5], [25, 14]), "offset": ([20, 0], [0, 0])},
                    {"box": Box([5, 5], [15, 24]), "offset": ([0, 0], [0, 0])},
                    {"box": Box([15, 5], [25, 24]), "offset": ([0, 0], [20, 0])},
                    {"box": Box([-5, 5], [5, 24]), "offset": ([-20, 0], [0, 0])},
                    {"box": Box([5, 5], [25, 14]), "offset": ([0, 0], [0, 0])},
                    {"box": Box([5, 15], [25, 24]), "offset": ([0, 0], [0, 20])},
                    {"box": Box([5, -5], [25, 4]), "offset": ([0, -20], [0, 0])},
                    {"box": Box([5, 5], [15, 14]), "offset": ([0, 0], [0, 0])},
                    {"box": Box([5, 15], [15, 24]), "offset": ([0, 0], [0, 20])},
                    {"box": Box([5, -5], [15, 4]), "offset": ([0, -20], [0, 0])},
                    {"box": Box([-5, 15], [5, 24]), "offset": ([0, 0], [-20, 20])},
                    {"box": Box([15, -5], [25, 4]), "offset": ([20, -20], [0, 0])},
                    {"box": Box([-5, 5], [5, 14]), "offset": ([0, 0], [-20, 0])},
                    {"box": Box([15, 5], [25, 14]), "offset": ([20, 0], [0, 0])},
                ]
            },
        ),
        (
            {
                "L0": [
                    Box2D(0, 4),
                    Box([0, 15], [4, 19]),
                    Box([15, 0], [19, 4]),
                    Box2D(15, 19),
                ]
            },
            {
                1: [  # level 0 same as previous
                    {"box": Box([-5, -5], [15, 4]), "offset": ([0, 0], [0, -40])},
                    {"box": Box([-5, 35], [15, 44]), "offset": ([0, 40], [0, 0])},
                    {"box": Box([-5, -5], [5, 14]), "offset": ([0, 0], [-40, 0])},
                    {"box": Box([35, -5], [45, 14]), "offset": ([40, 0], [0, 0])},
                    {"box": Box([-5, -5], [5, 4]), "offset": ([0, 0], [-40, -40])},
                    {"box": Box([35, 35], [45, 44]), "offset": ([40, 40], [0, 0])},
                    {"box": Box([-5, 35], [5, 44]), "offset": ([0, 0], [-40, 40])},
                    {"box": Box([35, -5], [45, 4]), "offset": ([40, -40], [0, 0])},
                    {"box": Box([-5, 25], [5, 44]), "offset": ([0, 0], [-40, 0])},
                    {"box": Box([35, 25], [45, 44]), "offset": ([40, 0], [0, 0])},
                    {"box": Box([25, -5], [45, 4]), "offset": ([0, 0], [0, -40])},
                    {"box": Box([25, 35], [45, 44]), "offset": ([0, 40], [0, 0])},
                ]
            },
        ),
    )
    @unpack
    def test_overlaps(self, refinement_boxes, expected):
        print("GeometryTest.test_overlaps")
        test_id = self._testMethodName.split("_")[-1]
        dim, interp_order, nbr_cells = (2, 1, [20] * 2)
        hierarchy = self.setup_hierarchy(
            dim, interp_order, nbr_cells, refinement_boxes, quantities="Bx"
        )

        level_overlaps = hierarchy_overlaps(hierarchy)

        for ilvl, lvl in enumerate(hierarchy.levels().items()):
            if ilvl not in expected:
                continue
            self.assertEqual(len(expected[ilvl]), len(level_overlaps[ilvl]))
            for exp, actual in zip(expected[ilvl], level_overlaps[ilvl]):
                self.assertEqual(actual["box"], exp["box"])
                self.assertTrue(
                    (np.asarray(actual["offset"]) == np.asarray(exp["offset"])).all()
                )

        if 1 in level_overlaps:
            fig = hierarchy.plot_2d_patches(
                1, [p.box for p in hierarchy.level(1).patches]
            )
            fig.savefig(
                "hierarchy_2d_lvl_" + str(1) + "_simple_test_" + str(test_id) + ".png"
            )

    @data(
        [
            {"box": Box([15, 15], [26, 25]), "offset": ([0, 0], [20, 20])},
            {"box": Box([-5, -5], [6, 5]), "offset": ([-20, -20], [0, 0])},
            {"box": Box([15, -5], [26, 25]), "offset": ([0, 0], [20, 0])},
            {"box": Box([-5, -5], [6, 25]), "offset": ([-20, 0], [0, 0])},
            {"box": Box([15, -5], [26, 5]), "offset": ([0, 0], [20, -20])},
            {"box": Box([-5, 15], [6, 25]), "offset": ([-20, 20], [0, 0])},
            {"box": Box([-5, 15], [26, 25]), "offset": ([0, 0], [0, 20])},
            {"box": Box([-5, -5], [26, 5]), "offset": ([0, -20], [0, 0])},
            {"box": Box([-5, -5], [26, 5]), "offset": ([0, 0], [0, -20])},
            {"box": Box([-5, 15], [26, 25]), "offset": ([0, 20], [0, 0])},
            {"box": Box([-5, 15], [6, 25]), "offset": ([0, 0], [-20, 20])},
            {"box": Box([15, -5], [26, 5]), "offset": ([20, -20], [0, 0])},
            {"box": Box([-5, -5], [6, 25]), "offset": ([0, 0], [-20, 0])},
            {"box": Box([15, -5], [26, 25]), "offset": ([20, 0], [0, 0])},
            {"box": Box([-5, -5], [6, 5]), "offset": ([0, 0], [-20, -20])},
            {"box": Box([15, 15], [26, 25]), "offset": ([20, 20], [0, 0])},
        ]
    )
    def test_large_patchoverlaps(self, expected):
        print(f"GeometryTest.{self._testMethodName}")
        test_id = self._testMethodName.split("_")[-1]
        dim, interp_order, nbr_cells = (2, 1, [20] * 2)
        hierarchy = self.setup_hierarchy(
            dim, interp_order, nbr_cells, {}, quantities="Bx", largest_patch_size=20
        )

        level_overlaps = hierarchy_overlaps(hierarchy)
        ilvl = 0
        overlap_boxes = []

        for exp, actual in zip(expected, level_overlaps[ilvl]):
            self.assertEqual(actual["box"], exp["box"])
            self.assertTrue(
                (np.asarray(actual["offset"]) == np.asarray(exp["offset"])).all()
            )
            overlap_boxes += [actual["box"]]

        fig = hierarchy.plot_2d_patches(
            ilvl,
            collections=[
                {
                    "boxes": overlap_boxes,
                    "facecolor": "yellow",
                },
                {
                    "boxes": [p.box for p in hierarchy.level(ilvl).patches],
                    "facecolor": "grey",
                },
            ],
        )
        fig.savefig(f"hierarchy_2d_lvl_0_large_patch_test_{test_id}.png")

    def test_touch_border(self):
        print("GeometryTest.test_touch_border")
        hierarchy = super().basic_hierarchy()

        domain_box = hierarchy.domain_box
        self.assertFalse(touch_domain_border(Box2D(10, 14), domain_box, "upper"))
        self.assertFalse(touch_domain_border(Box2D(10, 14), domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box2D(0, 14), domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box2D(-5, 14), domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box2D(-5, 10), domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box2D(-5, 30), domain_box, "upper"))
        self.assertTrue(touch_domain_border(Box2D(10, 30), domain_box, "upper"))
        self.assertTrue(touch_domain_border(Box2D(10, 24), domain_box, "upper"))

    def test_particle_ghost_area_boxes(self):
        print("GeometryTest.test_particle_ghost_area_boxes")
        hierarchy = super().basic_hierarchy()

        expected = {
            0: [
                {
                    "boxes": [
                        Box([-1, -1], [-1, 10]),
                        Box([10, -1], [10, 10]),
                        Box([0, -1], [9, -1]),
                        Box([0, 10], [9, 10]),
                    ]
                },
                {
                    "boxes": [
                        Box([9, 9], [9, 20]),
                        Box([20, 9], [20, 20]),
                        Box([10, 9], [19, 9]),
                        Box([10, 20], [19, 20]),
                    ]
                },
                {
                    "boxes": [
                        Box([-1, 9], [-1, 20]),
                        Box([10, 9], [10, 20]),
                        Box([0, 9], [9, 9]),
                        Box([0, 20], [9, 20]),
                    ]
                },
                {
                    "boxes": [
                        Box([9, -1], [9, 10]),
                        Box([20, -1], [20, 10]),
                        Box([10, -1], [19, -1]),
                        Box([10, 10], [19, 10]),
                    ]
                },
            ],
            1: [
                {
                    "boxes": [
                        Box([9, 9], [9, 20]),
                        Box([20, 9], [20, 20]),
                        Box([10, 9], [19, 9]),
                        Box([10, 20], [19, 20]),
                    ]
                },
                {
                    "boxes": [
                        Box([27, 27], [27, 40]),
                        Box([40, 27], [40, 40]),
                        Box([28, 27], [39, 27]),
                        Box([28, 40], [39, 40]),
                    ]
                },
            ],
        }

        gaboxes = ghost_area_boxes(hierarchy, "particles")
        particles = "particles"

        self.assertEqual(len(expected), len(gaboxes))

        for ilvl, lvl in enumerate(hierarchy.levels().items()):
            self.assertEqual(len(gaboxes[ilvl][particles]), len(expected[ilvl]))
            for act_pdata, exp_pdata in zip(gaboxes[ilvl][particles], expected[ilvl]):
                self.assertEqual(len(exp_pdata["boxes"]), len(act_pdata["boxes"]))
                for exp_box in exp_pdata["boxes"]:
                    self.assertTrue(exp_box in act_pdata["boxes"])

    def test_particle_level_ghost_boxes_do_not_overlap_patch_interiors(self):
        print(
            "GeometryTest.test_particle_level_ghost_boxes_do_not_overlap_patch_interiors"
        )
        hierarchy = super().basic_hierarchy()

        lvl_gboxes = level_ghost_boxes(hierarchy, "particles")

        assert len(lvl_gboxes) > 0
        for ilvl, pdatainfos in lvl_gboxes.items():
            assert len(pdatainfos) > 0
            for particles_id, gaboxes_list in pdatainfos.items():
                assert len(gaboxes_list) > 0
                for pdatainfo in gaboxes_list:
                    for box in pdatainfo["boxes"]:
                        for patch in hierarchy.level(ilvl).patches:
                            self.assertIsNone(patch.box * box)


@ddt
class ParticleLevelGhostGeometryTest(AGeometryTest):
    @data(
        (  # no patch ghost on level 1
            {"L0": [Box2D(5, 9), Box2D(14, 19)]},
            {
                1: [
                    Box([9, 9], [9, 20]),
                    Box([20, 9], [20, 20]),
                    Box([10, 9], [19, 9]),
                    Box([10, 20], [19, 20]),
                    Box([27, 27], [27, 40]),
                    Box([40, 27], [40, 40]),
                    Box([28, 27], [39, 27]),
                    Box([28, 40], [39, 40]),
                ]
            },
        ),
        (  # top right of gabox0 overlaps bottom left gabox1
            {"L0": [Box2D(5, 9), Box2D(10, 11)]},
            {
                1: [
                    Box([9, 9], [9, 20]),
                    Box([20, 9], [20, 19]),
                    Box([10, 9], [19, 9]),
                    Box([10, 20], [19, 20]),
                    Box([19, 20], [19, 24]),
                    Box([24, 19], [24, 24]),
                    Box([20, 19], [23, 19]),
                    Box([20, 24], [23, 24]),
                ]
            },
        ),
        (  # right side of gabox0 overlaps left side of gabox1
            {"L0": [Box2D(2, 9), Box([10, 1], [14, 10])]},
            {
                1: [
                    Box([3, 3], [3, 20]),
                    Box([4, 3], [19, 3]),
                    Box([4, 20], [19, 20]),
                    Box([19, 1], [19, 3]),
                    Box([19, 20], [19, 22]),
                    Box([30, 1], [30, 22]),
                    Box([20, 1], [29, 1]),
                    Box([20, 22], [29, 22]),
                ]
            },
        ),
        (
            {  # right side of gabox0 overlaps left of gabox1/gabox2
                # left side of gabox3 overlaps right of gabox1/gabox2
                "L0": [
                    Box([0, 0], [4, 19]),
                    Box([10, 0], [14, 19]),
                    Box([5, 2], [9, 6]),
                    Box([5, 12], [9, 16]),
                ]
            },
            {
                1: [
                    Box([-1, 1], [-1, 38]),
                    Box([10, 14], [10, 23]),
                    Box([10, 1], [10, 3]),
                    Box([10, 34], [10, 38]),
                    Box([19, 14], [19, 23]),
                    Box([19, 1], [19, 3]),
                    Box([19, 34], [19, 38]),
                    Box([30, 1], [30, 38]),
                    Box([10, 3], [19, 3]),
                    Box([10, 14], [19, 14]),
                    Box([10, 23], [19, 23]),
                    Box([10, 34], [19, 34]),
                ]
            },
        ),
    )
    @unpack
    def test_level_ghost_boxes(self, refinement_boxes, expected):
        print("ParticleLevelGhostGeometryTest.test_level_ghost_boxes")
        dim, interp_order, nbr_cells = (2, 1, [20] * 2)
        hierarchy = self.setup_hierarchy(
            dim, interp_order, nbr_cells, refinement_boxes, quantities="particles"
        )
        lvl_gaboxes = level_ghost_boxes(hierarchy, "particles")

        for ilvl in range(1, len(hierarchy.levels())):
            self.assertEqual(len(lvl_gaboxes[ilvl].keys()), 1)

            key = list(lvl_gaboxes[ilvl].keys())[0]

            ghost_area_box_list = sum(  # aggregate to single list
                [actual["boxes"] for actual in lvl_gaboxes[ilvl][key]], []
            )

            fig = hierarchy.plot_2d_patches(
                ilvl,
                collections=[
                    {
                        "boxes": ghost_area_box_list,
                        "facecolor": "yellow",
                    },
                    {
                        "boxes": [p.box for p in hierarchy.level(ilvl).patches],
                        "facecolor": "grey",
                    },
                ],
                title="".join(
                    [str(box) + "," for box in refinement_boxes["L" + str(ilvl - 1)]]
                ),
            )
            fig.savefig(f"{type(self).__name__}_lvl_{ilvl}_{self._testMethodName}.png")
            self.assertEqual(expected[ilvl], ghost_area_box_list)

    @data(
        (
            {
                # no level 1
            },
            {
                0: [
                    Box([0, 0], [9, 9]),
                    Box([0, 10], [9, 19]),
                    Box([10, 0], [19, 9]),
                    Box([10, 10], [19, 19]),
                    Box([-1, 19], [10, 30]),
                    Box([19, -1], [30, 10]),
                    Box([19, 19], [30, 30]),
                    Box([-11, -11], [0, 0]),
                    Box([-11, 9], [0, 20]),
                    Box([9, -11], [20, 0]),
                    Box([-1, -11], [10, 0]),
                    Box([19, -11], [30, 0]),
                    Box([19, 9], [30, 20]),
                    Box([-11, -1], [0, 10]),
                    Box([-11, 19], [0, 30]),
                    Box([9, 19], [20, 30]),
                ]
            },
        ),
        (
            {
                "L0": [
                    Box([0, 0], [4, 4]),
                    Box([15, 15], [19, 19]),
                    Box([15, 0], [19, 4]),
                    Box([0, 15], [4, 19]),
                ]
            },
            {
                1: [
                    Box([0, 0], [9, 9]),
                    Box([30, 0], [39, 9]),
                    Box([0, 30], [9, 39]),
                    Box([30, 30], [39, 39]),
                    Box([-1, 39], [10, 50]),
                    Box([39, -1], [50, 10]),
                    Box([39, 39], [50, 50]),
                    Box([-11, -11], [0, 0]),
                    Box([-11, 29], [0, 40]),
                    Box([29, -11], [40, 0]),
                    Box([-11, -1], [0, 10]),
                    Box([-11, 39], [0, 50]),
                    Box([29, 39], [40, 50]),
                    Box([-1, -11], [10, 0]),
                    Box([39, -11], [50, 0]),
                    Box([39, 29], [50, 40]),
                ]
            },
        ),
        (
            {
                "L0": [
                    Box([1, 1], [5, 5]),
                    Box([14, 14], [18, 18]),
                    Box([14, 1], [18, 5]),
                    Box([1, 14], [5, 18]),
                ]
            },
            {
                # NO lvl 1 PERIODIC PARTICLE GHOST BOX
            },
        ),
    )
    @unpack
    def test_patch_periodicity_copy(self, refinement_boxes, expected):
        print("ParticleLevelGhostGeometryTest.test_patch_periodicity_copy")
        dim, interp_order, nbr_cells = (2, 1, [20] * 2)
        hierarchy = self.setup_hierarchy(
            dim, interp_order, nbr_cells, refinement_boxes, quantities="particles"
        )

        for ilvl, lvl in hierarchy.levels().items():
            domain_box = hierarchy.level_domain_box(ilvl)

            qtyNbr = len(lvl.patches[0].patch_datas.keys())
            self.assertEqual(qtyNbr, 1)
            pop_name = list(lvl.patches[0].patch_datas.keys())[0]
            n_ghosts = lvl.patches[0].patch_datas[pop_name].ghosts_nbr

            periodic_list = get_periodic_list(lvl.patches, domain_box, n_ghosts)

            box_list = [p.box for p in periodic_list]

            if ilvl in expected:
                self.assertEqual(expected[ilvl], box_list)

            fig = hierarchy.plot_2d_patches(
                ilvl,
                collections=[
                    {
                        "boxes": [p.box for p in periodic_list],
                        "facecolor": "grey",
                    },
                ],
            )
            fig.savefig(f"{type(self).__name__}_lvl_{ilvl}_{self._testMethodName}.png")


if __name__ == "__main__":
    import matplotlib as mpl

    mpl.use("Agg")  # for cli only environments
    unittest.main()
