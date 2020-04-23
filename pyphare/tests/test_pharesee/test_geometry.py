import unittest

import pyphare.pharesee.box as boxm
from pyphare.pharesee.box import Box
from pyphare.pharesee.particles import Particles
from pyphare.pharesee.hierarchy import FieldData
from pyphare.pharesee.hierarchy import ParticleData
from pyphare.pharesee.hierarchy import PatchHierarchy
from pyphare.pharesee.hierarchy import Patch, PatchLevel
from pyphare.pharesee.geometry import particle_ghost_area_boxes
from pyphare.pharesee.geometry import level_ghost_boxes, hierarchy_overlaps, touch_domain_border
from pyphare.core.gridlayout import GridLayout

import numpy as np


def bx(ghost_box, dx, Lx, origin):
    x = origin + np.arange(ghost_box.size() + 1) * dx - 5 * dx
    return np.sin(2 * np.pi / Lx * x)


def by(ghost_box, dx, Lx, origin):
    x = origin + np.arange(ghost_box.size()) * dx - 5 * dx
    return np.cos(2 * np.pi / Lx * x)


def bz(ghost_box, dx, Lx, origin):
    x = origin + np.arange(ghost_box.size()) * dx - 5 * dx
    return np.sin(4 * np.pi / Lx * x)


def ex(ghost_box, dx, Lx, origin):
    x = origin + np.arange(ghost_box.size()) * dx - 5 * dx
    return np.sin(2 * np.pi / Lx * x)


def ey(ghost_box, dx, Lx, origin):
    x = origin + np.arange(ghost_box.size() + 1) * dx - 5 * dx
    return np.cos(2 * np.pi / Lx * x)


def ez(ghost_box, dx, Lx, origin):
    x = origin + np.arange(ghost_box.size() + 1) * dx - 5 * dx
    return np.sin(4 * np.pi / Lx * x)





def build_hierarchy(nbr_cells, origin, interp_order, domain_size, cell_width, refinement_ratio, refinement_boxes):

    domain_box = boxm.Box(0, nbr_cells-1)
    domain_layout = GridLayout(domain_box, origin, cell_width, interp_order)

    coarse_particles = Particles(box=domain_box)

    # copy domain particles an put them in ghost cells
    particle_ghost_nbr = domain_layout.particleGhostNbr(interp_order)
    box_extend         =  particle_ghost_nbr-1

    upper_slct_box = Box(domain_box.upper - box_extend, domain_box.upper)
    lower_slct_box = Box(domain_box.lower, domain_box.lower+box_extend)

    upper_cell_particles = coarse_particles.select(upper_slct_box)
    lower_cell_particles = coarse_particles.select(lower_slct_box)

    coarse_particles.add(upper_cell_particles.shift_icell(-domain_box.size()))
    coarse_particles.add(lower_cell_particles.shift_icell(domain_box.size()))


    boxes = {}
    for ilvl, boxes_data in refinement_boxes.items():

        level_number = ilvl+1

        if level_number not in boxes:
            boxes[level_number] = []


        for boxname, lower_upper in boxes_data.items():

            refinement_box = Box(lower_upper[0], lower_upper[1])
            refined_box    = boxm.refine(refinement_box, refinement_ratio)
            boxes[level_number].append(refined_box)


    # coarse level boxes are arbitrarily divided in 2 patches in the middle
    middle_cell = round(domain_box.upper/2)
    lower_box = Box(0,middle_cell)
    upper_box = Box(middle_cell+1, domain_box.upper)
    boxes[0] = [lower_box, upper_box]


    patch_datas = {}

    for ilvl, lvl_box in boxes.items():

        lvl_cell_width = cell_width/(refinement_ratio**ilvl)

        if ilvl not in patch_datas:
            patch_datas[ilvl] = []

        for box in lvl_box:

            ghost_box = boxm.grow(box, 5)
            origin = box.lower * lvl_cell_width
            layout = GridLayout(box, origin, lvl_cell_width, interp_order)

            datas = {"Bx"        : bx(ghost_box, cell_width, domain_size, origin),
                     "By"        : by(ghost_box, cell_width, domain_size, origin),
                     "Bz"        : bz(ghost_box, cell_width, domain_size, origin),
                     "Ex"        : ex(ghost_box, cell_width, domain_size, origin),
                     "Ey"        : ey(ghost_box, cell_width, domain_size, origin),
                     "Ez"        : ez(ghost_box, cell_width, domain_size, origin),
                     "particles" : coarse_particles.select(boxm.grow(box, layout.particleGhostNbr(interp_order)))
                     }

            boxed_patch_datas = {}
            for qty_name, data in datas.items():
                if qty_name == 'particles':
                    pdata = ParticleData(layout, data)
                else:
                    pdata = FieldData(layout, qty_name, data)

                boxed_patch_datas[qty_name] = pdata

            patch_datas[ilvl].append(boxed_patch_datas)



    patches = {}
    for ilvl, lvl_patch_datas in patch_datas.items():
        if ilvl not in patches :
            patches[ilvl] = []

        for patch_datas in lvl_patch_datas:
            patches[ilvl].append(Patch(patch_datas))



    patch_levels = {}
    for ilvl, lvl_patches in patches.items():
        patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

    sorted_levels_numbers = sorted(patch_levels)
    patch_levels = {ilvl: patch_levels[ilvl] for ilvl in sorted_levels_numbers}
    return PatchHierarchy(list(patch_levels.values()), domain_box, refinement_ratio)




class GeometryTest(unittest.TestCase):



    def setUp(self):




        nbr_cells = 65
        Lx = 1.
        dx = Lx / nbr_cells
        xtot = np.arange(nbr_cells) * dx
        ratio = 2
        domain_box = boxm.Box(0, 64)

        self.coarse_particles = Particles(box=domain_box)
        upper_cell_particles = self.coarse_particles.select(Box(domain_box.upper, domain_box.upper))
        lower_cell_particles = self.coarse_particles.select(Box(domain_box.lower, domain_box.lower))
        self.coarse_particles.add(upper_cell_particles.shift_icell(-domain_box.size()))
        self.coarse_particles.add(lower_cell_particles.shift_icell(domain_box.size()))

        self.fine_particles = Particles(box=boxm.refine(domain_box, ratio))  # should be split of coarse of course

        L0B1 = Box(0, 32)
        L0B2 = Box(33, 64)
        L1B1 = boxm.refine(Box(5, 29), ratio)
        L1B2 = boxm.refine(Box(32, 55), ratio)



        L0GB1 = boxm.grow(L0B1, 5)
        origin = 0.
        bx_dataL0P1 = bx(L0GB1, dx, Lx, origin)
        by_dataL0P1 = by(L0GB1, dx, Lx, origin)
        bz_dataL0P1 = bz(L0GB1, dx, Lx, origin)
        ex_dataL0P1 = ex(L0GB1, dx, Lx, origin)
        ey_dataL0P1 = ey(L0GB1, dx, Lx, origin)
        ez_dataL0P1 = ez(L0GB1, dx, Lx, origin)

        L0GB2 = boxm.grow(L0B2, 5)
        origin = xtot[33]
        bx_dataL0P2 = bx(L0GB2, dx, Lx, origin)
        by_dataL0P2 = by(L0GB2, dx, Lx, origin)
        bz_dataL0P2 = bz(L0GB2, dx, Lx, origin)
        ex_dataL0P2 = ex(L0GB2, dx, Lx, origin)
        ey_dataL0P2 = ey(L0GB2, dx, Lx, origin)
        ez_dataL0P2 = ez(L0GB2, dx, Lx, origin)

        L1GB1 = boxm.grow(L1B1, 5)
        origin = 5 * dx
        bx_dataL1P1 = bx(L1GB1, dx / ratio, Lx, origin)
        by_dataL1P1 = by(L1GB1, dx / ratio, Lx, origin)
        bz_dataL1P1 = bz(L1GB1, dx / ratio, Lx, origin)
        ex_dataL1P1 = ex(L1GB1, dx / ratio, Lx, origin)
        ey_dataL1P1 = ey(L1GB1, dx / ratio, Lx, origin)
        ez_dataL1P1 = ez(L1GB1, dx / ratio, Lx, origin)

        L1GB2 = boxm.grow(L1B2, 5)
        origin = 32 * dx
        bx_dataL1P2 = bx(L1GB2, dx / ratio, Lx, origin)
        by_dataL1P2 = by(L1GB2, dx / ratio, Lx, origin)
        bz_dataL1P2 = bz(L1GB2, dx / ratio, Lx, origin)
        ex_dataL1P2 = ex(L1GB2, dx / ratio, Lx, origin)
        ey_dataL1P2 = ey(L1GB2, dx / ratio, Lx, origin)
        ez_dataL1P2 = ez(L1GB2, dx / ratio, Lx, origin)

        interp_order = 1
        layout = GridLayout(L0B1, 0., dx, interp_order)

        self.L0P1_datas = {"EM_B_x": FieldData(layout, "Bx", bx_dataL0P1),
                      "EM_B_y": FieldData(layout, "By", by_dataL0P1),
                      "EM_B_z": FieldData(layout, "Bz", bz_dataL0P1),
                      "EM_E_x": FieldData(layout, "Ex", ex_dataL0P1),
                      "EM_E_y": FieldData(layout, "Ey", ey_dataL0P1),
                      "EM_E_z": FieldData(layout, "Ez", ez_dataL0P1),
                      "particles": ParticleData(layout, self.coarse_particles.select(boxm.grow(layout.box,1)))
                      }

        layout = GridLayout(L0B2, xtot[33], dx, interp_order)
        self.L0P2_datas = {"EM_B_x": FieldData(layout, "Bx", bx_dataL0P2),
                      "EM_B_y": FieldData(layout, "By", by_dataL0P2),
                      "EM_B_z": FieldData(layout, "Bz", bz_dataL0P2),
                      "EM_E_x": FieldData(layout, "Ex", ex_dataL0P2),
                      "EM_E_y": FieldData(layout, "Ey", ey_dataL0P2),
                      "EM_E_z": FieldData(layout, "Ez", ez_dataL0P2),
                      "particles": ParticleData(layout, self.coarse_particles.select(boxm.grow(layout.box, 1)))
                      }

        layout = GridLayout(L1B1, 5 * dx, dx / ratio, interp_order)
        self.L1P1_datas = {"EM_B_x": FieldData(layout, "Bx", bx_dataL1P1),
                      "EM_B_y": FieldData(layout, "By", by_dataL1P1),
                      "EM_B_z": FieldData(layout, "Bz", bz_dataL1P1),
                      "EM_E_x": FieldData(layout, "Ex", ex_dataL1P1),
                      "EM_E_y": FieldData(layout, "Ey", ey_dataL1P1),
                      "EM_E_z": FieldData(layout, "Ez", ez_dataL1P1),
                      "particles": ParticleData(layout, self.fine_particles.select(boxm.grow(layout.box, 1)))
                      }

        layout = GridLayout(L1B2, 32 * dx, dx / ratio, interp_order)
        self.L1P2_datas = {"EM_B_x": FieldData(layout, "Bx", bx_dataL1P2),
                      "EM_B_y": FieldData(layout, "By", by_dataL1P2),
                      "EM_B_z": FieldData(layout, "Bz", bz_dataL1P2),
                      "EM_E_x": FieldData(layout, "Ex", ex_dataL1P2),
                      "EM_E_y": FieldData(layout, "Ey", ey_dataL1P2),
                      "EM_E_z": FieldData(layout, "Ez", ez_dataL1P2),
                      "particles": ParticleData(layout, self.fine_particles.select(boxm.grow(layout.box, 1)))
                      }

        self.L0P1 = Patch(self.L0P1_datas)
        self.L0P2 = Patch(self.L0P2_datas)

        self.L1P1 = Patch(self.L1P1_datas)
        self.L1P2 = Patch(self.L1P2_datas)

        self.L0 = PatchLevel(0, [self.L0P1, self.L0P2])
        self.L1 = PatchLevel(1, [self.L1P1, self.L1P2])
        self.hierarchy = PatchHierarchy([self.L0, self.L1], domain_box, ratio)



    def test_overlaps(self):
        expected = {0 :
                        # Middle overlap, for all quantities
                        [  {
                            "pdatas":[self.L0P1_datas["EM_B_x"], self.L0P2_datas["EM_B_x"]],
                            "box":Box(28, 38),
                            'offset':(0,0)
                          },

                          {
                             "pdatas": [self.L0P1_datas["EM_B_y"], self.L0P2_datas["EM_B_y"]],
                            "box": Box(28, 37),
                              "offset":(0,0)
                          },
                        {
                            "pdatas": [self.L0P1_datas["EM_B_z"], self.L0P2_datas["EM_B_z"]],
                            "box": Box(28, 37),
                            "offset":(0,0)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_E_x"], self.L0P2_datas["EM_E_x"]],
                            "box": Box(28, 37),
                            'offset': (0,0)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_E_y"], self.L0P2_datas["EM_E_y"]],
                            "box": Box(28, 38),
                            "offset":(0,0)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_z"], self.L0P2_datas["EM_E_z"]],
                            "box": Box(28, 38),
                            "offset":(0,0)
                        },

                        {
                            "pdatas": [self.L0P1_datas["particles"], self.L0P2_datas["particles"]],
                            "box": Box(32,33),
                            "offset":(0,0)
                        },

                        # left side overlap with periodicity, for all quantities
                          {
                            "pdatas": [self.L0P1_datas["EM_B_x"], self.L0P2_datas["EM_B_x"]],
                            "box": Box(-5, 5),
                              "offset":(0,-65)
                          },
                        # right side overlap with periodicity, for all quantities
                        {
                            "pdatas": [self.L0P1_datas["EM_B_x"], self.L0P2_datas["EM_B_x"]],
                            "box": Box(60, 70),
                            "offset": (65,0)
                        },


                        {
                            "pdatas": [self.L0P1_datas["EM_B_y"], self.L0P2_datas["EM_B_y"]],
                            "box": Box(-5, 4),
                            "offset":(0,-65)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_B_y"], self.L0P2_datas["EM_B_y"]],
                            "box": Box(60, 69),
                            "offset": (65,0 )
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_B_z"], self.L0P2_datas["EM_B_z"]],
                            "box": Box(-5, 4),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_B_z"], self.L0P2_datas["EM_B_z"]],
                            "box": Box(60, 69),
                            "offset": (65,0)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_E_x"], self.L0P2_datas["EM_E_x"]],
                            "box": Box(-5, 4),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_x"], self.L0P2_datas["EM_E_x"]],
                            "box": Box(60, 69),
                            "offset":(65, 0)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_E_y"], self.L0P2_datas["EM_E_y"]],
                            "box": Box(-5, 5),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_y"], self.L0P2_datas["EM_E_y"]],
                            "box": Box(60, 70),
                            "offset":(65, 0)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_z"], self.L0P2_datas["EM_E_z"]],
                            "box": Box(-5, 5),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_z"], self.L0P2_datas["EM_E_z"]],
                            "box": Box(60, 70),
                            "offset":(65, 0)
                        },
                        {
                            "pdatas": [self.L0P1_datas["particles"], self.L0P2_datas["particles"]],
                            "box": Box(-1,0),
                            "offset":(0,-65)
                        },

                        {
                            "pdatas": [self.L0P1_datas["particles"], self.L0P2_datas["particles"]],
                            "box": Box(64, 65),
                            "offset":(65 ,0)
                        }
                    ],

            1:  # level 1
                [
                    {
                        "pdatas": [self.L0P1_datas["EM_B_x"], self.L0P2_datas["EM_B_x"]],
                        "box": Box(59, 65),
                        'offset': (0, 0)
                    },

                    {
                        "pdatas": [self.L0P1_datas["EM_B_y"], self.L0P2_datas["EM_B_y"]],
                        "box": Box(59, 64),
                        "offset": (0, 0)
                    },
                    {
                        "pdatas": [self.L0P1_datas["EM_B_z"], self.L0P2_datas["EM_B_z"]],
                        "box": Box(59, 64),
                        "offset": (0, 0)
                    },

                    {
                        "pdatas": [self.L0P1_datas["EM_E_x"], self.L0P2_datas["EM_E_x"]],
                        "box": Box(59, 64),
                        'offset': (0, 0)
                    },

                    {
                        "pdatas": [self.L0P1_datas["EM_E_y"], self.L0P2_datas["EM_E_y"]],
                        "box": Box(59, 65),
                        "offset": (0, 0)
                    },
                    {
                        "pdatas": [self.L0P1_datas["EM_E_z"], self.L0P2_datas["EM_E_z"]],
                        "box": Box(59, 65),
                        "offset": (0, 0)
                    }

                ]


        }


        overlaps = hierarchy_overlaps(self.hierarchy)

        for ilvl, lvl in enumerate(self.hierarchy.patch_levels):
            self.assertEqual(len(expected[ilvl]), len(overlaps[ilvl]))

            for exp, actual in zip(expected[ilvl], overlaps[ilvl]):

                act_box    = actual["box"]
                act_pdatas = actual["pdatas"]
                act_offset = actual["offset"]

                exp_box    = exp["box"]
                exp_pdatas = exp["pdatas"]
                exp_offset = exp["offset"]


                self.assertEqual(act_box, exp_box)
                self.assertEqual(act_offset, exp_offset)
                #self.assertEqual(act_box, exp_box)
                #self.assertEqual(act_box, exp_box)




    def test_touch_border(self):

        self.assertFalse(touch_domain_border(Box(10,20), self.hierarchy.domain_box, "upper"))
        self.assertFalse(touch_domain_border(Box(10, 20), self.hierarchy.domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box(0, 20), self.hierarchy.domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box(-5, 20), self.hierarchy.domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box(-5, 70), self.hierarchy.domain_box, "lower"))
        self.assertTrue(touch_domain_border(Box(-5, 70), self.hierarchy.domain_box, "upper"))
        self.assertTrue(touch_domain_border(Box(40, 70), self.hierarchy.domain_box, "upper"))
        self.assertTrue(touch_domain_border(Box(40, 64), self.hierarchy.domain_box, "upper"))




    def test_particle_ghost_area_boxes(self):

        expected = {

            0: [
                {
                    "pdatas":self.L0P1_datas["particles"],
                    "boxes":[Box(33, 33), Box(-1,-1)]
                },
                {
                    "pdatas": self.L0P2_datas["particles"],
                    "boxes": [Box(32, 32), Box(65, 65)]
                }
            ],


            1: [
                {
                    "pdatas": self.L1P1_datas["particles"],
                    "boxes": [Box(9, 9), Box(60, 60)]
                },
                {
                    "pdatas": self.L1P2_datas["particles"],
                    "boxes": [Box(63, 63), Box(112, 112)]
                }

            ]
        }

        gaboxes = particle_ghost_area_boxes(self.hierarchy)

        # same number of levels
        self.assertEqual(len(expected), len(gaboxes))

        for ilvl, lvl in enumerate(self.hierarchy.patch_levels):

            # same number of PatchDatas
            self.assertEqual(len(gaboxes[ilvl]), len(expected[ilvl]))

            for act_pdata, exp_pdata in zip(gaboxes[ilvl], expected[ilvl]):

                    self.assertEqual(len(exp_pdata["boxes"]), len(act_pdata["boxes"]))

                    for exp_box in exp_pdata["boxes"]:
                        self.assertTrue(exp_box in act_pdata["boxes"])





    def test_level_ghost_boxes(self):

        expected = {


            1: [
                {
                    "pdatas": self.L1P1_datas["particles"],
                    "boxes": [Box(9, 9),]
                },
                {
                    "pdatas": self.L1P1_datas["particles"],
                    "boxes": [Box(60, 60),]
                },
                {
                    "pdatas": self.L1P2_datas["particles"],
                    "boxes": [Box(63, 63),]
                },
                {
                    "pdatas": self.L1P2_datas["particles"],
                    "boxes": [Box(112, 112),]
                }
            ]
        }

        lvl_gaboxes = level_ghost_boxes(self.hierarchy)
        for ilvl  in range(1, len(self.hierarchy.patch_levels)):
            for actual, exp in zip(lvl_gaboxes[ilvl], expected[ilvl]):

                act_boxes = actual["boxes"]
                exp_boxes = exp["boxes"]

                for act_box, exp_box in zip(act_boxes, exp_boxes):
                    self.assertEqual(act_box, exp_box)





    def test_overlapped_fields_are_equal(self):

        overlaps = hierarchy_overlaps(self.hierarchy)

        for ilvl, lvl in enumerate(self.hierarchy.patch_levels):

            for overlap in overlaps[ilvl]:

                pd1, pd2 = overlap["pdatas"]
                box      = overlap["box"]
                offsets  = overlap["offset"]

                self.assertEqual(pd1.quantity, pd2.quantity)

                if pd1.quantity == 'field':

                    # we need to transform the AMR overlap box, which is thus
                    # (because AMR) common to both pd1 and pd2 into local index
                    # boxes that will allow to slice the data

                    # the patchData ghost box that serves as a reference box
                    # to transfrom AMR to local indexes first needs to be
                    # shifted by the overlap offset associated to it
                    # this is because the overlap box has been calculated from
                    # the intersection of possibly shifted patch data ghost boxes

                    loc_b1 = boxm.amr_to_local(box, boxm.shift(pd1.ghost_box, offsets[0]))
                    loc_b2 = boxm.amr_to_local(box, boxm.shift(pd2.ghost_box, offsets[1]))

                    data1 = pd1.dataset
                    data2 = pd2.dataset

                    slice1 = data1[loc_b1.lower:loc_b1.upper+1]
                    slice2 = data2[loc_b2.lower:loc_b2.upper+ 1]

                    self.assertTrue(np.allclose(slice1, slice2, atol=1e-12))







    def test_patch_particle_overlap(self):

        from copy import copy
        overlaps = hierarchy_overlaps(self.hierarchy)

        for ilvl, lvl in enumerate(self.hierarchy.patch_levels):

            for overlap in overlaps[ilvl]:

                pd1, pd2 = overlap["pdatas"]
                box      = overlap["box"]
                offsets  = overlap["offset"]

                self.assertEqual(pd1.quantity, pd2.quantity)

                if pd1.quantity == "particles":

                    # the following uses 'offset', we need to remember that offset
                    # is the quantity by which a patch has been moved to detect
                    # overlap with the other one.
                    # so shift by +offset when evaluating patch data in overlap box
                    # index space, and by -offset when we want to shift box indexes
                    # to the associated patch index space.

                    # overlap box must be shifted by -offset to select data in the patches
                    part1 = copy(pd1.dataset.select(boxm.shift(box, -offsets[0])))
                    part2 = copy(pd2.dataset.select(boxm.shift(box, -offsets[1])))

                    idx1 = np.argsort(part1.iCells + part1.deltas)
                    idx2 = np.argsort(part2.iCells + part2.deltas)


                    # particle iCells are in their patch AMR space
                    # so we need to shift them by +offset to move them to the box space
                    np.testing.assert_array_equal(part1.iCells[idx1]+offsets[0], part2.iCells[idx2]+offsets[1])

                    self.assertTrue(np.allclose(part1.deltas[idx1], part2.deltas[idx2], atol=1e-12))
                    self.assertTrue(np.allclose(part1.vx[idx1], part2.vx[idx2], atol=1e-12))
                    self.assertTrue(np.allclose(part1.vy[idx1], part2.vy[idx2], atol=1e-12))
                    self.assertTrue(np.allclose(part1.vz[idx1], part2.vz[idx2], atol=1e-12))






    def test_hier(self):

        nbr_cells = 65
        origin = 0.
        interp_order = 1
        domain_size = 1.
        cell_width = domain_size / nbr_cells
        refinement_ratio = 2
        refinement_boxes = {0: {"B0": (5, 29), "B1": (32, 55)}}

        hier = build_hierarchy(nbr_cells,
                               origin,
                               interp_order,
                               domain_size,
                               cell_width,
                               refinement_ratio,
                               refinement_boxes)

        print(hier)








