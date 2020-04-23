

import unittest
from ddt import ddt, data, unpack
from pyphare.pharesee import box as boxm
from pyphare.pharesee.box import Box


@ddt
class BoxTesT(unittest.TestCase):


    @data((Box(10,20),Box(10,20)) ,
           (Box(-5, 5), Box(-5, 5)),
           (Box(-5, -2), Box(-5, -2)))
    @unpack
    def test_box_equality(self, box1, box2):
        self.assertTrue(box1 == box2)



    @data( (Box(10,20),Box(20,41)) ,
           (Box(-5,11), Box(-10,23)),
           (Box(-5, -2), Box(-10, -3)),
           (Box(1, 1), Box(2, 3)),
           )
    @unpack
    def test_box_refinement(self, coarse, exp_refined):
        actual_refined = boxm.refine(coarse, 2)
        self.assertEqual(actual_refined, exp_refined)



    @data((Box(10,20), Box(15,30), Box(15,20)),        # upper intersection
          (Box(-5, 20), Box(5, 5), Box(5, 5)),         # box2 within and 1 cell
          (Box(-5, -5), Box(-5, -5), Box(-5, -5)),     # all negative & 1 cell
          (Box(-5, 10), Box(-10, -5), Box(-5, -5)),    # negative and box2 lower inter.
          (Box(10, 11), Box(11, 12), Box(11, 11)),     # positive upper cell inter
          (Box(10, 11), Box(14, 15), None),            # no inter
          )
    @unpack
    def test_intersection(self, box1, box2, exp_inter):
        self.assertEqual(box1*box2, exp_inter)



    @data( (Box(10, 20),  5 , Box(15, 25)),
           (Box(10, 20), -5 , Box( 5, 15)),
           (Box(10, 20), -15, Box(-5,  5)),
           (Box(-5, 20),  10, Box( 5, 30)),
           (Box(-5, -5), -2,  Box(-7,-7))
           )
    @unpack
    def test_shift(self, box, shift, expected):
        self.assertEqual(boxm.shift(box, shift), expected)



    @data( (Box(10, 20),  5 , Box( 5, 25)),
           (Box(-5, 20),  10, Box(-15, 30)),
           (Box(-5, -5),  2,  Box(-7,-3)),
           )
    @unpack
    def test_grow(self, box, size, expected):
        self.assertEqual(boxm.grow(box, size), expected)



    def test_grow_neg_size_raises(self):
        with self.assertRaises(ValueError):
            boxm.grow(Box(10,12),-1)





    @data((Box(10, 30), Box(15, 25), [Box(10,14),Box(26,30)]),   # remove middle part
          (Box(10, 30), Box(5, 20), [Box(21, 30),]),             # remove lower part
          (Box(10, 30), Box(15, 35), [Box(10, 14), ]),           # remove upper part
          (Box(10, 30), Box(5, 35), []),                         # remove all
          (Box(-10, 30), Box(-5, 25), [Box(-10, -6), Box(26, 30)]),  # remove middle part
          (Box(-10, 30), Box(-15, 20), [Box(21, 30), ]),             # remove lower part
          (Box(-10, 30), Box(15, 35), [Box(-10, 14), ]),             # remove upper part
          (Box(-10, 30), Box(-5, 35), []),                           # remove all

          )
    @unpack
    def test_remove(self, box, to_remove, expected):


        keeps = boxm.remove(box, to_remove)

        if len(keeps) == 0:
            self.assertEqual(len(expected), 0)

        for exp, keep in zip(expected, keeps):
            self.assertEqual(keep, exp)



    @data( (Box(10, 20), 10, True),
           (Box(10, 20), 11, True),
           (Box(10, 20), 20, True),
           (Box(-20, -10), -10, True),
           (Box(-20, -10), -20, True),
           (Box(-20, -10), -11, True),
           (Box(-20, -10), -9, False),
           (Box(10, 20), Box(10,11), True),
           (Box(10, 20), Box(10, 10), True),
           (Box(10, 20), Box(15, 20), True),
           (Box(10, 20), Box(20, 20), True),
           (Box(10, 20), Box(20, 21), False),
            )
    @unpack
    def test_in(self, box, element, expected):
        self.assertEqual(element in box, expected)



    def test_amr_to_local(self):
        self.assertEqual(boxm.amr_to_local(Box(33,64),Box(-5,69)), Box(38,69))