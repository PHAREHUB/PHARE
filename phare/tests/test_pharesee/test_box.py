

import unittest
from ddt import ddt, data, unpack
from pharesee import box as boxm
from pharesee.box import Box


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



    @data((Box(10,20), Box(15,30), Box(15,20)),
          (Box(10, 20), Box(15, 30), Box(15, 20)),
          (Box(-5, 20), Box(5, 5), Box(5, 5))
          )
    @unpack
    def test_intersection(self, box1, box2, exp_inter):
        print("assert = {}*{} = {} and get {}".format(box1, box2, exp_inter, box1*box2))
        self.assertEqual(box1*box2, exp_inter)




