import numpy as np

import pharesee.box as boxm
from pharesee.box import Box
from pharesee.particles import Particles
from pharesee.hierarchy import FieldData
from pharesee.hierarchy import ParticleData
from pharesee.hierarchy import PatchHierarchy
from pharesee.hierarchy import Patch, PatchLevel
from pharesee.geometry import compute_overlaps, toFieldBox
from pharesee.geometry import particle_ghost_area_boxes
from pharesee.geometry import level_ghost_boxes, hierarchy_overlaps
from core.gridlayout import GridLayout


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

coarse_particles = Particles(box = domain_box)
fine_particles = Particles(box = boxm.refine(domain_box, ratio))  # should be split of coarse of course

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

L0P1_datas = {"EM_B_x": FieldData(L0B1, layout, "Bx", bx_dataL0P1),
              "EM_B_y": FieldData(L0B1, layout, "By", by_dataL0P1),
              "EM_B_z": FieldData(L0B1, layout, "Bz", bz_dataL0P1),
              "EM_E_x": FieldData(L0B1, layout, "Ex", ex_dataL0P1),
              "EM_E_y": FieldData(L0B1, layout, "Ey", ey_dataL0P1),
              "EM_E_z": FieldData(L0B1, layout, "Ez", ez_dataL0P1),
              "particles": ParticleData(L0B1, layout, select_particles(coarse_particles, L0B1))
              }

layout =  GridLayout(xtot[33], dx, interp_order)
L0P2_datas = {"EM_B_x": FieldData(L0B2, layout, "Bx", bx_dataL0P2),
              "EM_B_y": FieldData(L0B2, layout, "By", by_dataL0P2),
              "EM_B_z": FieldData(L0B2, layout, "Bz", bz_dataL0P2),
              "EM_E_x": FieldData(L0B2, layout, "Ex", ex_dataL0P2),
              "EM_E_y": FieldData(L0B2, layout, "Ey", ey_dataL0P2),
              "EM_E_z": FieldData(L0B2, layout, "Ez", ez_dataL0P2),
              "particles": ParticleData(L0B2, layout, select_particles(coarse_particles, L0B2))
              }

layout = GridLayout(5 * dx, dx / ratio, interp_order)
L1P1_datas = {"EM_B_x": FieldData(L1B1, layout, "Bx", bx_dataL1P1),
              "EM_B_y": FieldData(L1B1, layout, "By", by_dataL1P1),
              "EM_B_z": FieldData(L1B1, layout, "Bz", bz_dataL1P1),
              "EM_E_x": FieldData(L1B1, layout, "Ex", ex_dataL1P1),
              "EM_E_y": FieldData(L1B1, layout, "Ey", ey_dataL1P1),
              "EM_E_z": FieldData(L1B1, layout, "Ez", ez_dataL1P1),
              "particles": ParticleData(L1B1, layout, select_particles(fine_particles, L1B1))
              }

layout = GridLayout(32 * dx, dx / ratio, interp_order)
L1P2_datas = {"EM_B_x": FieldData(L1B2, layout, "Bx", bx_dataL1P2),
              "EM_B_y": FieldData(L1B2, layout, "By", by_dataL1P2),
              "EM_B_z": FieldData(L1B2, layout, "Bz", bz_dataL1P2),
              "EM_E_x": FieldData(L1B2, layout, "Ex", ex_dataL1P2),
              "EM_E_y": FieldData(L1B2, layout, "Ey", ey_dataL1P2),
              "EM_E_z": FieldData(L1B2, layout, "Ez", ez_dataL1P2),
              "particles": ParticleData(L1B2, layout, select_particles(fine_particles, L1B2))
              }


L0P1 = Patch(L0B1, GridLayout(0., dx, interp_order), patch_datas=L0P1_datas)
L0P2 = Patch(L0B2, GridLayout(xtot[33], dx, interp_order), patch_datas=L0P2_datas)

L1P1 = Patch(L1B1, GridLayout(5 * dx, dx / ratio, interp_order), patch_datas=L1P1_datas)
L1P2 = Patch(L1B2, GridLayout(32 * dx, dx / ratio, interp_order), patch_datas=L1P2_datas)

L0 = PatchLevel(0, [L0P1, L0P2])
L1 = PatchLevel(1, [L1P1, L1P2])
hierarchy = PatchHierarchy([L0, L1], domain_box, ratio)

print(hierarchy)

print("")
print("ghost area boxes per level")

gaboxes = particle_ghost_area_boxes(hierarchy)

for ilvl, infos in gaboxes.items():
    for info in infos:
        boxes = info["boxes"]
        for b in boxes:
            print("level {} has box {}".format(ilvl, b))

lvl_gab = level_ghost_boxes(hierarchy)
print("\n")
print("level ghost area boxes per level")

for ilvl, infos in lvl_gab.items():
    for info in infos:
        boxes = info["boxes"]
        for b in boxes:
            print("level {} has box {}".format(ilvl, b))

print("")


def print_overlaps(hierarchy):
    for ilvl, lvl in enumerate(hierarchy.patch_levels):
        refined_domain_box = boxm.refine(domain_box, ratio ** (ilvl))
        # print("refined domain box {}**{}: {}".format(ratio, ilvl,refined_domain_box))
        overlaps = compute_overlaps(lvl.patches, refined_domain_box)
        print("level {} overlap : ".format(ilvl))
        for overlap in overlaps:
            box = overlap["box"]
            pdatas = overlap["pdatas"]
            if pdatas[0].quantity == "field":
                print("    overlap for {} : {}".format(pdatas[0].field_name, toFieldBox(box, pdatas[0])))
            else:
                print("    overlap for particles : {}".format(box))


print_overlaps(hierarchy)
print("")

overlaps = hierarchy_overlaps(hierarchy)

for ilvl, lvl in enumerate(hierarchy.patch_levels):
    print("level {} overlap data comparison".format(ilvl))

    for overlap in overlaps[ilvl]:

        pd1, pd2 = overlap["pdatas"]
        box = overlap["box"]

        if pd1.quantity == 'field':
            pd_boxes = [toFieldBox(box, pd) for pd in (pd1, pd2)]
        else:
            pd_boxes = [box, box]

        offsets = overlap["offsets"]
        local_boxes = [boxm.amr_to_local(b, boxm.shift(pd.ghost_box, offset)) for b, pd, offset in
                       zip(pd_boxes, (pd1, pd2), offsets)]

        print("pd_box = {} and local boxes = {} and {}".format(pd_boxes[0],
                                                               local_boxes[0], local_boxes[1]))

        if pd1.quantity == "field":
            datas = (pd.dataset for pd in (pd1, pd2))
            slices = [data[box.lower:box.upper + 1] for data, box in zip(datas, local_boxes)]

            print("    {} points on patch data {}".format(len(slices[0]), pd1.field_name))
            print("    data are not zero: {} ".format(np.sum(np.abs(slices[0])) > 1e-8))
            print("    data are the same down to 1e-12: {}\n".format(np.allclose(slices[0], slices[1], atol=1e-12)))
        else:
            print("we don't compare particles for now...\n")

print("")





