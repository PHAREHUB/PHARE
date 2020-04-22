import unittest

import pyphare.pharesee.box as boxm
from pyphare.pharesee.box import Box
from pyphare.pharesee.particles import Particles
from pyphare.pharesee.hierarchy import FieldData
from pyphare.pharesee.hierarchy import ParticleData
from pyphare.pharesee.hierarchy import PatchHierarchy
from pyphare.pharesee.hierarchy import Patch, PatchLevel
from pyphare.pharesee.geometry import compute_overlaps, toFieldBox
from pyphare.pharesee.geometry import particle_ghost_area_boxes
from pyphare.pharesee.geometry import level_ghost_boxes, hierarchy_overlaps, touch_domain_border
from pyphare.core.gridlayout import GridLayout

import numpy as np



class GeometryTest(unittest.TestCase):

    def setUp(self):


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

        def select_particles(particles, box):
            idx = np.where((particles.iCells >= box.lower) & (particles.iCells <= box.upper))[0]
            return Particles(icell=particles.iCells[idx],
                             deltas=particles.deltas[idx, :],
                             vx=particles.vx[idx, :],
                             vy=particles.vy[idx, :],
                             vz=particles.vz[idx, :])

        nbr_cells = 65
        Lx = 1.
        dx = Lx / nbr_cells
        xtot = np.arange(nbr_cells) * dx
        ratio = 2
        domain_box = boxm.Box(0, 64)

        coarse_particles = Particles(box=domain_box)
        fine_particles = Particles(box=boxm.refine(domain_box, ratio))  # should be split of coarse of course

        L0B1 = Box(0, 32)
        L0B2 = Box(33, 64)
        L1B1 = boxm.refine(Box(5, 29), ratio)
        L1B2 = boxm.refine(Box(32, 55), ratio)

        print("box 1 and 2 to refine are {0} and {1}".format(L0B1, L0B2))
        print("refined boxes are : {0} and  {1}".format(L1B1, L1B2))
        print("")

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
        layout = GridLayout(0., dx, interp_order)

        self.L0P1_datas = {"EM_B_x": FieldData(L0B1, layout, "Bx", bx_dataL0P1),
                      "EM_B_y": FieldData(L0B1, layout, "By", by_dataL0P1),
                      "EM_B_z": FieldData(L0B1, layout, "Bz", bz_dataL0P1),
                      "EM_E_x": FieldData(L0B1, layout, "Ex", ex_dataL0P1),
                      "EM_E_y": FieldData(L0B1, layout, "Ey", ey_dataL0P1),
                      "EM_E_z": FieldData(L0B1, layout, "Ez", ez_dataL0P1),
                      "particles": ParticleData(L0B1, layout, select_particles(coarse_particles, L0B1))
                      }

        layout = GridLayout(xtot[33], dx, interp_order)
        self.L0P2_datas = {"EM_B_x": FieldData(L0B2, layout, "Bx", bx_dataL0P2),
                      "EM_B_y": FieldData(L0B2, layout, "By", by_dataL0P2),
                      "EM_B_z": FieldData(L0B2, layout, "Bz", bz_dataL0P2),
                      "EM_E_x": FieldData(L0B2, layout, "Ex", ex_dataL0P2),
                      "EM_E_y": FieldData(L0B2, layout, "Ey", ey_dataL0P2),
                      "EM_E_z": FieldData(L0B2, layout, "Ez", ez_dataL0P2),
                      "particles": ParticleData(L0B2, layout, select_particles(coarse_particles, L0B2))
                      }

        layout = GridLayout(5 * dx, dx / ratio, interp_order)
        self.L1P1_datas = {"EM_B_x": FieldData(L1B1, layout, "Bx", bx_dataL1P1),
                      "EM_B_y": FieldData(L1B1, layout, "By", by_dataL1P1),
                      "EM_B_z": FieldData(L1B1, layout, "Bz", bz_dataL1P1),
                      "EM_E_x": FieldData(L1B1, layout, "Ex", ex_dataL1P1),
                      "EM_E_y": FieldData(L1B1, layout, "Ey", ey_dataL1P1),
                      "EM_E_z": FieldData(L1B1, layout, "Ez", ez_dataL1P1),
                      "particles": ParticleData(L1B1, layout, select_particles(fine_particles, L1B1))
                      }

        layout = GridLayout(32 * dx, dx / ratio, interp_order)
        self.L1P2_datas = {"EM_B_x": FieldData(L1B2, layout, "Bx", bx_dataL1P2),
                      "EM_B_y": FieldData(L1B2, layout, "By", by_dataL1P2),
                      "EM_B_z": FieldData(L1B2, layout, "Bz", bz_dataL1P2),
                      "EM_E_x": FieldData(L1B2, layout, "Ex", ex_dataL1P2),
                      "EM_E_y": FieldData(L1B2, layout, "Ey", ey_dataL1P2),
                      "EM_E_z": FieldData(L1B2, layout, "Ez", ez_dataL1P2),
                      "particles": ParticleData(L1B2, layout, select_particles(fine_particles, L1B2))
                      }

        self.L0P1 = Patch(L0B1, GridLayout(0., dx, interp_order), patch_datas = self.L0P1_datas)
        self.L0P2 = Patch(L0B2, GridLayout(xtot[33], dx, interp_order), patch_datas = self.L0P2_datas)

        self.L1P1 = Patch(L1B1, GridLayout(5 * dx, dx / ratio, interp_order), patch_datas = self.L1P1_datas)
        self.L1P2 = Patch(L1B2, GridLayout(32 * dx, dx / ratio, interp_order), patch_datas = self.L1P2_datas)

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
                            "offset": (0,65)
                        },


                        {
                            "pdatas": [self.L0P1_datas["EM_B_y"], self.L0P2_datas["EM_B_y"]],
                            "box": Box(-5, 4),
                            "offset":(0,-65)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_B_y"], self.L0P2_datas["EM_B_y"]],
                            "box": Box(60, 69),
                            "offset": (0,65)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_B_z"], self.L0P2_datas["EM_B_z"]],
                            "box": Box(-5, 4),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_B_z"], self.L0P2_datas["EM_B_z"]],
                            "box": Box(60, 69),
                            "offset": (0,65)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_E_x"], self.L0P2_datas["EM_E_x"]],
                            "box": Box(-5, 4),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_x"], self.L0P2_datas["EM_E_x"]],
                            "box": Box(60, 69),
                            "offset":(0,65)
                        },

                        {
                            "pdatas": [self.L0P1_datas["EM_E_y"], self.L0P2_datas["EM_E_y"]],
                            "box": Box(-5, 5),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_y"], self.L0P2_datas["EM_E_y"]],
                            "box": Box(60, 70),
                            "offset":(0,65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_z"], self.L0P2_datas["EM_E_z"]],
                            "box": Box(-5, 5),
                            "offset":(0,-65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["EM_E_z"], self.L0P2_datas["EM_E_z"]],
                            "box": Box(60, 70),
                            "offset":(0,65)
                        },
                        {
                            "pdatas": [self.L0P1_datas["particles"], self.L0P2_datas["particles"]],
                            "box": Box(-1,0),
                            "offset":(0,-65)
                        },

                        {
                            "pdatas": [self.L0P1_datas["particles"], self.L0P2_datas["particles"]],
                            "box": Box(64, 65),
                            "offset":(0,65)
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
            refined_domain_box = self.hierarchy.refined_domain_box(ilvl)


            overlaps = compute_overlaps(lvl.patches, refined_domain_box)

            print("level {} overlap : ".format(ilvl))

            self.assertEqual(len(expected[ilvl]), len(overlaps))

            for exp, actual in zip(expected[ilvl], overlaps):

                act_box    = actual["box"]
                act_pdatas = actual["pdatas"]
                act_offset = actual["offset"]

                exp_box    = exp["box"]
                exp_pdatas = exp["pdatas"]
                exp_offset = exp["offset"]

                if exp_pdatas[0].quantity == 'field':
                    print("    overlap for {} : {}".format(exp_pdatas[0].field_name, exp_box))
                else:
                    print("    overlap for particles : {}".format(exp_box))


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
