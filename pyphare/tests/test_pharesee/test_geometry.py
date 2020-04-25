import unittest

from ddt import ddt, data
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


def build_hierarchy(**kwargs):
    """accepted keywords:
     - simulation : a simulation Object

     or

     - nbr_cells
     - origin
     - interp_order
     - domain_size
     - cell_width
     - refinement_ratio
     - refinement_boxes
     """
    if "simulation" in kwargs:
        for k in kwargs:
            if k != "simulation":
                print("warning: 'simulation' given, {} discarded".format(k))

        sim = kwargs["simulation"]
        nbr_cells = sim.cells[0]
        origin = sim.origin[0]
        interp_order = sim.interp_order
        domain_size = sim.simulation_domain()[0]
        cell_width = sim.dl[0]
        refinement_ratio = 2
        refinement_boxes = sim.refinement_boxes

    else:
        nbr_cells = kwargs["nbr_cells"]
        origin = kwargs["origin"]
        interp_order = kwargs["interp_order"]
        domain_size = kwargs["domain_size"]
        cell_width = kwargs["cell_width"]
        refinement_ratio = kwargs.get("refinement_ratio", 1)
        refinement_boxes = kwargs.get("refinement_boxes", {})

    domain_box = boxm.Box(0, nbr_cells - 1)
    domain_layout = GridLayout(domain_box, origin, cell_width, interp_order)

    coarse_particles = Particles(box=domain_box)

    # copy domain particles an put them in ghost cells
    particle_ghost_nbr = domain_layout.particleGhostNbr(interp_order)
    box_extend = particle_ghost_nbr - 1

    upper_slct_box = Box(domain_box.upper - box_extend, domain_box.upper)
    lower_slct_box = Box(domain_box.lower, domain_box.lower + box_extend)

    upper_cell_particles = coarse_particles.select(upper_slct_box)
    lower_cell_particles = coarse_particles.select(lower_slct_box)

    coarse_particles.add(upper_cell_particles.shift_icell(-domain_box.size()))
    coarse_particles.add(lower_cell_particles.shift_icell(domain_box.size()))

    boxes = {}
    for ilvl, boxes_data in refinement_boxes.items():

        level_number = int(ilvl.strip("L")) + 1

        if level_number not in boxes:
            boxes[level_number] = []

        for boxname, lower_upper in boxes_data.items():
            refinement_box = Box(lower_upper[0][0], lower_upper[1][0])
            refined_box = boxm.refine(refinement_box, refinement_ratio)
            boxes[level_number].append(refined_box)

    # coarse level boxes are arbitrarily divided in 2 patches in the middle
    middle_cell = round(domain_box.upper / 2)
    lower_box = Box(0, middle_cell)
    upper_box = Box(middle_cell + 1, domain_box.upper)
    boxes[0] = [lower_box, upper_box]

    patch_datas = {}

    for ilvl, lvl_box in boxes.items():

        lvl_cell_width = cell_width / (refinement_ratio ** ilvl)

        if ilvl == 0:
            lvl_particles = coarse_particles
        else:
            level_domain_box = boxm.refine(domain_box, refinement_ratio)
            lvl_ghost_domain_box = boxm.grow(level_domain_box, domain_layout.particleGhostNbr(interp_order))
            lvl_particles = Particles(box = lvl_ghost_domain_box)

        if ilvl not in patch_datas:
            patch_datas[ilvl] = []

        for box in lvl_box:

            ghost_box = boxm.grow(box, 5)
            origin = box.lower * lvl_cell_width
            layout = GridLayout(box, origin, lvl_cell_width, interp_order)

            datas = {"Bx": bx(ghost_box, lvl_cell_width, domain_size, origin),
                     "By": by(ghost_box, lvl_cell_width, domain_size, origin),
                     "Bz": bz(ghost_box, lvl_cell_width, domain_size, origin),
                     "Ex": ex(ghost_box, lvl_cell_width, domain_size, origin),
                     "Ey": ey(ghost_box, lvl_cell_width, domain_size, origin),
                     "Ez": ez(ghost_box, lvl_cell_width, domain_size, origin),
                     "particles": lvl_particles.select(ghost_box)
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
        if ilvl not in patches:
            patches[ilvl] = []

        for patch_datas in lvl_patch_datas:
            patches[ilvl].append(Patch(patch_datas))

    patch_levels = {}
    for ilvl, lvl_patches in patches.items():
        patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

    sorted_levels_numbers = sorted(patch_levels)
    patch_levels = {ilvl: patch_levels[ilvl] for ilvl in sorted_levels_numbers}
    return PatchHierarchy(list(patch_levels.values()), domain_box, refinement_ratio)




@ddt
class GeometryTest(unittest.TestCase):

    def setUp(self):

        nbr_cells = 65
        origin = 0.
        interp_order = 1
        domain_size = 1.
        cell_width = domain_size / nbr_cells
        refinement_ratio = 2
        refinement_boxes = {"L0": {"B0": [(5,), (29,)], "B1": [(32,), (55,)]}}

        self.hierarchy = build_hierarchy(nbr_cells=nbr_cells,
                               origin=origin,
                               interp_order=interp_order,
                               domain_size=domain_size,
                               cell_width=cell_width,
                               refinement_ratio=refinement_ratio,
                               refinement_boxes=refinement_boxes)




    def test_overlaps(self):
        expected = {0 :
                        # Middle overlap, for all quantities
                        [  {"box":Box(28, 38),'offset':(0,0)},

                        {"box": Box(28, 37),"offset":(0,0)},
                        {"box": Box(28, 37),"offset":(0,0)},
                        {"box": Box(28, 37),'offset': (0,0)},
                        {"box": Box(28, 38),"offset":(0,0)},
                        {"box": Box(28, 38),"offset":(0,0)},
                        {"box": Box(32,33),"offset":(0,0)},

                        # left side overlap with periodicity, for all quantities
                          {"box": Box(-5, 5),"offset":(0,-65)},
                        # right side overlap with periodicity, for all quantities
                        {"box": Box(60, 70),"offset": (65,0)},
                        {"box": Box(-5, 4),"offset":(0,-65)},
                        {"box": Box(60, 69),"offset": (65,0 )},
                        {"box": Box(-5, 4),"offset":(0,-65)},
                        {"box": Box(60, 69),"offset": (65,0)},
                        {"box": Box(-5, 4),"offset":(0,-65)},
                        {"box": Box(60, 69),"offset":(65, 0)},
                        {"box": Box(-5, 5),"offset":(0,-65)},
                        {"box": Box(60, 70),"offset":(65, 0)},
                        {"box": Box(-5, 5),"offset":(0,-65)},
                        {"box": Box(60, 70),"offset":(65, 0)},
                        {"box": Box(-1,0),"offset":(0,-65)},
                        {"box": Box(64, 65),"offset":(65 ,0)}
                    ],

            1:  # level 1
                [
                    {"box": Box(59, 65),'offset': (0, 0)},
                    {"box": Box(59, 64),"offset": (0, 0)},
                    {"box": Box(59, 64),"offset": (0, 0)},
                    {"box": Box(59, 64),'offset': (0, 0)},
                    {"box": Box(59, 65),"offset": (0, 0)},
                    {"box": Box(59, 65),"offset": (0, 0)}
                ]


        }


        overlaps = hierarchy_overlaps(self.hierarchy)

        for ilvl, lvl in enumerate(self.hierarchy.patch_levels):
            self.assertEqual(len(expected[ilvl]), len(overlaps[ilvl]))

            for exp, actual in zip(expected[ilvl], overlaps[ilvl]):

                act_box    = actual["box"]
                act_offset = actual["offset"]
                exp_box    = exp["box"]
                exp_offset = exp["offset"]

                self.assertEqual(act_box, exp_box)
                self.assertEqual(act_offset, exp_offset)




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

            0: [ {"boxes":[Box(33, 33), Box(-1,-1)]},
                 {"boxes": [Box(32, 32), Box(65, 65)]} ],


            1: [ {"boxes": [Box(9, 9), Box(60, 60)]},
                 {"boxes": [Box(63, 63), Box(112, 112)]} ]
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
                {"boxes": [Box(9, 9),]},
                {"boxes": [Box(60, 60),] },
                {"boxes": [Box(63, 63),]},
                {"boxes": [Box(112, 112),]}
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







    def test_overlapped_particledatas_have_identical_particles(self):

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




    def test_level_ghost_boxes_do_not_overlap_patch_interiors(self):

        lvl_gboxes = level_ghost_boxes(self.hierarchy)

        for ilvl, pdatainfos in lvl_gboxes.items():
            for pdatainfo in pdatainfos:
                for box in pdatainfo["boxes"]:
                    for patch in self.hierarchy.patch_levels[ilvl].patches:
                            self.assertIsNone(patch.box * box)




    @data(1,2,3)
    def test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(self, interp_order):

        nbr_cells = 65
        origin = 0.
        domain_size = 1.
        cell_width = domain_size / nbr_cells
        refinement_ratio = 2

        # one refined cell between the refined patches
        refinement_boxes = {"L0": {"B0": [(5,), (30,)], "B1": [(31,), (55,)]}}

        hier = build_hierarchy(nbr_cells=nbr_cells,
                               origin=origin,
                               interp_order=interp_order,
                               domain_size=domain_size,
                               cell_width=cell_width,
                               refinement_ratio=refinement_ratio,
                               refinement_boxes=refinement_boxes)

        overlaps = hierarchy_overlaps(hier)
        for ilvl, lvl_overlaps in overlaps.items():
            for overlap in lvl_overlaps:

                if overlap["pdatas"][0].quantity == "particles":
                    ref_pd, cmp_pd = overlap["pdatas"]

                    if ilvl == 0:
                        box = overlap["box"]
                        offsets = overlap["offset"]

                        # first let's shift the overlap box over the AMR
                        # indices of the patchdata. The box has been created
                        # by shifting the patdata ghost box by 'offset' so here
                        # the box is shifted by -offset to get over patchdata
                        shift_refbox, shift_cmpbox = [boxm.shift(box, -off) for off in offsets]

                        # the overlap box overlaps both ghost and domain cells
                        # we need to extract the domain ones to later select domain
                        # particles
                        ovlped_refdom = ref_pd.box * shift_refbox
                        ovlped_cmpdom = cmp_pd.box * shift_cmpbox

                        # on lvl 0 patches are adjacent
                        # therefore the overlap box must overlap the
                        # patchData box. 1 cell in interporder1, 2 cells for higher
                        self.assertIsNotNone(ovlped_refdom)
                        self.assertIsNotNone(ovlped_cmpdom)

                        refdomain = ref_pd.dataset.select(ovlped_refdom)
                        cmpdomain = cmp_pd.dataset.select(ovlped_cmpdom)

                        # now get the ghost cells of each patch data overlaped by
                        # the overlap box. To do this we need to intersect the shifted
                        # overlap box with the patchdata ghost box, and remove interior cells
                        # note that in 1D we don't expect remove to return more than 1 box, hence [0]
                        ovlped_refghost = boxm.remove(ref_pd.ghost_box * shift_refbox, ref_pd.box)[0]
                        ovlped_cmpghost = boxm.remove(cmp_pd.ghost_box * shift_cmpbox, cmp_pd.box)[0]

                        refghost  = ref_pd.dataset.select(ovlped_refghost)
                        cmpghost  = cmp_pd.dataset.select(ovlped_cmpghost)


                        # before comparing the particles we need to be sure particles of both patchdatas
                        # are sorted in the same order. We do that by sorting by x position
                        sort_refdomain_idx = np.argsort(refdomain.iCells + refdomain.deltas)
                        sort_cmpdomain_idx = np.argsort(cmpdomain.iCells + cmpdomain.deltas)
                        sort_refghost_idx = np.argsort(refghost.iCells + refghost.deltas)
                        sort_cmpghost_idx = np.argsort(cmpghost.iCells + cmpghost.deltas)

                        np.testing.assert_allclose(refdomain.deltas[sort_refdomain_idx], cmpghost.deltas[sort_cmpghost_idx], atol=1e-12)
                        np.testing.assert_allclose(cmpdomain.deltas[sort_cmpdomain_idx], refghost.deltas[sort_refghost_idx], atol=1e-12)









    def test_hier(self):

        # some draft function

        nbr_cells = 65
        origin = 0.
        interp_order = 2
        domain_size = 1.
        cell_width = domain_size / nbr_cells
        refinement_ratio = 2
        refinement_boxes = {"L0": {"B0": [(5,), (29,)], "B1": [(32,), (55,)]}}

        hier = build_hierarchy(nbr_cells=nbr_cells,
                               origin=origin,
                               interp_order=interp_order,
                               domain_size=domain_size,
                               cell_width=cell_width,
                               refinement_ratio=refinement_ratio,
                               refinement_boxes=refinement_boxes)


        #print(hier)

        nbr_cells = 65
        origin = 0.
        interp_order = 2
        domain_size = 1.
        cell_width = domain_size / nbr_cells
        refinement_ratio = 2
        refinement_boxes = {"L0": {"B0": [(5,), (29,)],  "B1": [(32,), (55,)]},
                            "L1": {"B0": [(15,),(28,)],  "B1": [(33,), (47,)], "B2": [(66,),(101,)]},
                            "L2": {"B0": [(38,),(52,)],  "B1": [(68,), (80,)], 'B3': [(84,),(95,)], "B4": [(134,),(148,)]}}

        hier = build_hierarchy(nbr_cells=nbr_cells,
                               origin=origin,
                               interp_order=interp_order,
                               domain_size=domain_size,
                               cell_width=cell_width,
                               refinement_ratio=refinement_ratio,
                               refinement_boxes=refinement_boxes)


        #print(hier)

        nbr_cells = 65
        origin = 0.
        interp_order = 2
        domain_size = 1.
        cell_width = domain_size / nbr_cells
        refinement_ratio = 2
        refinement_boxes = {"L0": {"B0": [(5,), (29,)], "B1": [(32,), (55,)]}}

        hier = build_hierarchy(nbr_cells=nbr_cells,
                               origin=origin,
                               interp_order=interp_order,
                               domain_size=domain_size,
                               cell_width=cell_width,
                               refinement_ratio=refinement_ratio,
                               refinement_boxes=refinement_boxes)

