import unittest
from ddt import ddt, data, unpack
from pyphare.core import box as boxm
from pyphare.core.box import Box, Box1D, Box2D, Box3D


@ddt
class BoxTesT(unittest.TestCase):
    @data(
        (Box(10, 20), Box(10, 20)),
        (Box(-5, 5), Box(-5, 5)),
        (Box(-5, -2), Box(-5, -2)),
        (Box([10, 10], [20, 20]), Box([10, 10], [20, 20])),
        (Box([-5, -5], [5, 5]), Box([-5, -5], [5, 5])),
        (Box([-5, -5], [-2, -2]), Box([-5, -5], [-2, -2])),
    )
    @unpack
    def test_box_equality(self, box1, box2):
        self.assertTrue(box1 == box2)

    @data(
        (Box(10, 20), Box(20, 41)),
        (Box(-5, 11), Box(-10, 23)),
        (Box(-5, -2), Box(-10, -3)),
        (Box(1, 1), Box(2, 3)),
        (Box([10, 10], [20, 20]), Box([20, 20], [41, 41])),
        (Box([-5, -5], [11, 11]), Box([-10, -10], [23, 23])),
        (Box([-5, -5], [-2, -2]), Box([-10, -10], [-3, -3])),
        (Box([1, 1], [1, 1]), Box([2, 2], [3, 3])),
    )
    @unpack
    def test_box_refinement(self, coarse, exp_refined):
        actual_refined = boxm.refine(coarse, 2)
        self.assertEqual(actual_refined, exp_refined)

    @data(
        (Box(10, 20), Box(15, 30), Box(15, 20)),  # upper intersection
        (Box(-5, 20), Box(5, 5), Box(5, 5)),  # box2 within and 1 cell
        (Box(-5, -5), Box(-5, -5), Box(-5, -5)),  # all negative & 1 cell
        (Box(-5, 10), Box(-10, -5), Box(-5, -5)),  # negative and box2 lower inter.
        (Box(10, 11), Box(11, 12), Box(11, 11)),  # positive upper cell inter
        (Box(10, 11), Box(14, 15), None),  # no inter
        # upper intersection
        (Box([10, 10], [20, 20]), Box([15, 15], [30, 30]), Box([15, 15], [20, 20])),
        # box2 within and 1 cell
        (Box([-5, -5], [20, 20]), Box([5, 5], [5, 5]), Box([5, 5], [5, 5])),
        # all negative & 1 cell
        (Box([-5, -5], [-5, -5]), Box([-5, -5], [-5, -5]), Box([-5, -5], [-5, -5])),
        # negative and box2 lower inter.
        (Box([-5, -5], [10, 10]), Box([-10, -10], [-5, -5]), Box([-5, -5], [-5, -5])),
        # positive upper cell inter
        (Box([10, 10], [11, 11]), Box([11, 11], [12, 12]), Box([11, 11], [11, 11])),
        # no inter
        (Box1D(10, 11), Box1D(14, 15), None),
        (Box2D(10, 11), Box2D(14, 15), None),
        (Box2D(5, 24), Box2D(23, 44), Box2D(23, 24)),
        (Box3D(10, 11), Box3D(14, 15), None),
    )
    @unpack
    def test_intersection(self, box1, box2, exp_inter):
        self.assertEqual(box1 * box2, exp_inter)

    @data(
        (Box(10, 20), 5, Box(15, 25)),
        (Box(10, 20), -5, Box(5, 15)),
        (Box(10, 20), -15, Box(-5, 5)),
        (Box(-5, 20), 10, Box(5, 30)),
        (Box(-5, -5), -2, Box(-7, -7)),
        (Box([10, 10], [20, 20]), 5, Box([15, 15], [25, 25])),
        (Box([10, 10], [20, 20]), [5, 0], Box([15, 10], [25, 20])),
        (Box([10, 10], [20, 20]), [0, 5], Box([10, 15], [20, 25])),
        (Box([10, 10], [20, 20]), -5, Box([5, 5], [15, 15])),
        (Box([10, 10], [20, 20]), -15, Box([-5, -5], [5, 5])),
        (Box([-5, -5], [20, 20]), 10, Box([5, 5], [30, 30])),
        (Box([-5, -5], [-5, -5]), -2, Box([-7, -7], [-7, -7])),
    )
    @unpack
    def test_shift(self, box, shift, expected):
        self.assertEqual(boxm.shift(box, shift), expected)

    @data(
        (Box(10, 20), 5, Box(5, 25)),
        (Box(-5, 20), 10, Box(-15, 30)),
        (Box(-5, -5), 2, Box(-7, -3)),
        (Box([10, 10], [20, 20]), [5, 5], Box([5, 5], [25, 25])),
        (Box([-5, -5], [20, 20]), [10, 10], Box([-15, -15], [30, 30])),
        (Box([-5, -5], [-5, -5]), [2, 2], Box([-7, -7], [-3, -3])),
    )
    @unpack
    def test_grow(self, box, size, expected):
        self.assertEqual(boxm.grow(box, [size] * box.ndim), expected)

    @data(
        (Box(10, 12)),
        (Box([10, 10], [12, 12])),
    )
    def test_grow_neg_size_raises(self, box):
        with self.assertRaises(ValueError):
            boxm.grow(box, [-1] * box.ndim)

    @data(
        (Box1D(10, 29), Box1D(10, 25), [Box1D(26, 29)]),
        (Box1D(10, 29), Box1D(25, 29), [Box1D(10, 24)]),
        (Box1D(10, 29), Box1D(15, 19), [Box1D(10, 14), Box1D(20, 29)]),
        (Box1D(10, 29), Box1D(20, 29), [Box1D(10, 19)]),
        (Box(10, 30), Box(15, 25), [Box(10, 14), Box(26, 30)]),  # remove middle part
        (
            Box(10, 30),
            Box(5, 20),
            [
                Box(21, 30),
            ],
        ),  # remove lower part
        (
            Box(10, 30),
            Box(15, 35),
            [
                Box(10, 14),
            ],
        ),  # remove upper part
        (Box(10, 30), Box(5, 35), []),  # remove all
        (Box(-10, 30), Box(-5, 25), [Box(-10, -6), Box(26, 30)]),  # remove middle part
        (
            Box(-10, 30),
            Box(-15, 20),
            [
                Box(21, 30),
            ],
        ),  # remove lower part
        (
            Box(-10, 30),
            Box(15, 35),
            [
                Box(-10, 14),
            ],
        ),  # remove upper part
        (Box(-10, 30), Box(-15, 35), []),  # remove all
        (Box2D(10, 29), Box([10, 10], [25, 35]), [Box([26, 10], [29, 29])]),
        (Box2D(10, 29), Box([10, 10], [35, 25]), [Box([10, 26], [29, 29])]),
        (
            Box2D(10, 29),
            Box2D(10, 25),
            [  # remove down left
                Box([26, 10], [29, 29]),  # right
                Box([10, 26], [25, 29]),  # up
            ],
        ),
        (
            Box2D(10, 29),
            Box2D(20, 29),
            [  # remove up right
                Box([10, 10], [19, 29]),  # left
                Box([20, 10], [29, 19]),  # down
            ],
        ),
        (
            Box2D(10, 29),
            Box2D(15, 25),
            [  # remove middle
                Box([10, 10], [14, 29]),  # left
                Box([26, 10], [29, 29]),  # right
                Box([15, 10], [25, 14]),  # down
                Box([15, 26], [25, 29]),  # up
            ],
        ),
        (
            Box2D(10, 29),
            Box2D(20, 25),
            [  # remove middle
                Box([10, 10], [19, 29]),  # left
                Box([26, 10], [29, 29]),  # right
                Box([20, 10], [25, 19]),  # down
                Box([20, 26], [25, 29]),  # up
            ],
        ),
        (
            Box2D(10, 29),
            Box([15, 10], [15, 25]),
            [  # enters one side
                Box([10, 10], [14, 29]),  # left
                Box([16, 10], [29, 29]),  # right
                Box([15, 26], [15, 29]),  # up
            ],
        ),
        (
            Box3D(10, 29),
            Box3D(10, 15),
            [  # remove down left back
                Box([16, 10, 10], [29, 29, 29]),  # right
                Box([10, 16, 10], [15, 29, 29]),  # up
                Box([10, 10, 16], [15, 15, 29]),  # front
            ],
        ),
        (
            Box3D(10, 29),
            Box3D(20, 25),
            [  # remove middle
                Box([10, 10, 10], [19, 29, 29]),  # left
                Box([26, 10, 10], [29, 29, 29]),  # right
                Box([20, 10, 10], [25, 19, 29]),  # down
                Box([20, 26, 10], [25, 29, 29]),  # up
                Box([20, 20, 10], [25, 25, 19]),  # back
                Box([20, 20, 26], [25, 25, 29]),  # front
            ],
        ),
        (
            Box3D(10, 29),
            Box([20, 10, 20], [20, 29, 20]),
            [  # remove middle Y block
                Box([10, 10, 10], [19, 29, 29]),  # left
                Box([21, 10, 10], [29, 29, 29]),  # right
                Box([20, 10, 10], [20, 29, 19]),  # back
                Box([20, 10, 21], [20, 29, 29]),  # front
            ],
        ),
        (
            Box3D(10, 29),
            Box([10, 10, 10], [10, 29, 29]),
            [  # remove left
                Box([11, 10, 10], [29, 29, 29]),  # right
            ],
        ),
        (
            Box3D(10, 29),
            Box([29, 10, 10], [29, 29, 29]),
            [  # remove right
                Box([10, 10, 10], [28, 29, 29]),  # left
            ],
        ),
        (
            Box3D(10, 29),
            Box([10, 10, 10], [29, 29, 10]),
            [  # remove back
                Box([10, 10, 11], [29, 29, 29]),  # front
            ],
        ),
        (
            Box3D(10, 29),
            Box([10, 10, 29], [29, 29, 29]),
            [  # remove front
                Box([10, 10, 10], [29, 29, 28]),  # back
            ],
        ),
    )
    @unpack
    def test_remove(self, box, to_remove, expected):
        remaining = boxm.remove(box, to_remove)
        self.assertEqual(expected, remaining)

        for i, k0 in enumerate(remaining):
            self.assertTrue(k0 in box)
            self.assertTrue(k0 not in to_remove)
            for k1 in remaining[i + 1 :]:
                self.assertTrue(k0 not in k1)

        expectedCells, remainingCells = 0, 0
        for exp, keep in zip(expected, remaining):
            expectedCells += exp.nCells()
            remainingCells += keep.nCells()

        self.assertTrue(expectedCells == remainingCells)
        self.assertTrue(box.nCells() == remainingCells + (box * to_remove).nCells())

    @data(
        (Box(10, 20), 10, True),
        (Box(10, 20), 11, True),
        (Box(10, 20), 20, True),
        (Box(-20, -10), -10, True),
        (Box(-20, -10), -20, True),
        (Box(-20, -10), -11, True),
        (Box(-20, -10), -9, False),
        (Box(10, 20), Box(10, 11), True),
        (Box(10, 20), Box(10, 10), True),
        (Box(10, 20), Box(15, 20), True),
        (Box(10, 20), Box(20, 20), True),
        (Box(10, 20), Box(20, 21), False),
        (Box(10, 20), Box(20, 21), False),
        (Box([10, 10], [20, 20]), [10, 10], True),
        (Box([10, 10], [20, 20]), [11, 11], True),
        (Box([10, 10], [20, 20]), [20, 20], True),
        (Box([-20, -20], [-10, -10]), [-10, -10], True),
        (Box([-20, -20], [-10, -10]), [-20, -20], True),
        (Box([-20, -20], [-10, -10]), [-11, -11], True),
        (Box([-20, -20], [-10, -10]), [-9, -9], False),
        (Box([10, 10], [20, 20]), Box([10, 10], [11, 11]), True),
        (Box([10, 10], [20, 20]), Box([10, 10], [10, 10]), True),
        (Box([10, 10], [20, 20]), Box([15, 15], [20, 20]), True),
        (Box([10, 10], [20, 20]), Box([20, 20], [20, 20]), True),
        (Box2D(10, 20), Box2D(21, 21), False),
        (Box3D(10, 20), Box3D(20, 20), True),
        (Box3D(10, 20), Box3D(21, 21), False),
    )
    @unpack
    def test_in(self, box, element, expected):
        self.assertEqual(element in box, expected)

    @data(
        (Box(33, 64), Box(-5, 69), Box(38, 69)),
        (Box([33, 33], [64, 64]), Box([-5, -5], [69, 69]), Box([38, 38], [69, 69])),
    )
    @unpack
    def test_amr_to_local(self, box, ref_box, expected):
        self.assertEqual(boxm.amr_to_local(box, ref_box), expected)


if __name__ == "__main__":
    unittest.main()
