import numpy as np

import pyphare.core.box as boxm
from pyphare.core.box import Box
from pyphare.core.phare_utilities import listify
from pyphare.core.gridlayout import GridLayout, yee_element_is_primal

from pyphare.pharesee.particles import Particles

from pyphare.pharesee.hierarchy import FieldData
from pyphare.pharesee.hierarchy import ParticleData
from pyphare.pharesee.hierarchy import PatchHierarchy
from pyphare.pharesee.hierarchy import Patch, PatchLevel


def init(ghost_box, layout, L, qty, fn):

    assert layout.impl == "yee"

    ndim = ghost_box.dim()

    dl, origin = layout.dl, layout.origin
    xyz, directions = [], ["x", "y", "z"][: ndim]
    for dimdex, direction in enumerate(directions):
        primal = yee_element_is_primal(qty, direction)
        nbrGhosts = layout.nbrGhosts(primal, layout.interp_order)
        xyz.append(
            origin[dimdex]
            + np.arange(ghost_box.shape()[dimdex] + primal) * dl[dimdex]
            - nbrGhosts * dl[dimdex]
        )

    if ndim > 1:
        xx, yy = np.meshgrid(*xyz, indexing="ij")

        if ndim == 2:
            return fn(0, xx) * fn(1, yy)

        raise ValueError("3d unsupported")

    return fn(0, xyz[0])


def Bx(ghost_box, layout, L):
    return init(ghost_box, layout, L, "Bx", lambda i, v : np.sin(2 * np.pi / L[i] * v))


def By(ghost_box, layout, L):
    return init(ghost_box, layout, L, "By", lambda i, v : np.cos(2 * np.pi / L[i] * v))


def Bz(ghost_box, layout, L):
    return init(ghost_box, layout, L, "Bz", lambda i, v : np.sin(4 * np.pi / L[i] * v))


def Ex(ghost_box, layout, L):
    return init(ghost_box, layout, L, "Ex", lambda i, v : np.sin(2 * np.pi / L[i] * v))


def Ey(ghost_box, layout, L):
    return init(ghost_box, layout, L, "Ey", lambda i, v : np.cos(2 * np.pi / L[i] * v))


def Ez(ghost_box, layout, L):
    return init(ghost_box, layout, L, "Ez", lambda i, v : np.sin(4 * np.pi / L[i] * v))


def build_refinementboxes(domain_box, **kwargs):

    refinement_ratio = kwargs["refinement_ratio"]

    boxes = {}
    for ilvl, boxes_data in kwargs["refinement_boxes"].items():

        level_number = int(ilvl.strip("L")) + 1

        if level_number not in boxes:
            boxes[level_number] = []

        assert isinstance(boxes_data, (list, dict))

        if isinstance(boxes_data, list):
            for _, refinement_box in enumerate(boxes_data):
                refined_box = boxm.refine(refinement_box, refinement_ratio)
                boxes[level_number].append(refined_box)

        else:
            for boxname, lower_upper in boxes_data.items():
                refinement_box = Box(lower_upper[0][0], lower_upper[1][0])
                refined_box = boxm.refine(refinement_box, refinement_ratio)
                boxes[level_number].append(refined_box)

    # coarse level boxes are arbitrarily divided in 2 patches in the middle
    middle_cell = np.round(domain_box.upper / 2)
    lower_box = Box(0, middle_cell)
    upper_box = Box(middle_cell + 1, domain_box.upper)
    boxes[0] = [lower_box, upper_box]

    return boxes


def build_patch_datas(domain_box, boxes, **kwargs):
    quantities = kwargs["quantities"]
    origin = kwargs["origin"]
    interp_order = kwargs["interp_order"]
    domain_size = kwargs["domain_size"]
    cell_width = kwargs["cell_width"]
    refinement_ratio = kwargs["refinement_ratio"]

    skip_particles = False
    if "particles" in quantities:
        del quantities[quantities.index("particles")]
    else:
        skip_particles = True

    domain_layout = GridLayout(domain_box, origin, cell_width, interp_order)

    coarse_particles = Particles(box=domain_box)

    # copy domain particles and put them in ghost cells
    particle_ghost_nbr = domain_layout.particleGhostNbr(interp_order)
    box_extend = particle_ghost_nbr - 1

    upper_slct_box = Box(domain_box.upper - box_extend, domain_box.upper)
    lower_slct_box = Box(domain_box.lower, domain_box.lower + box_extend)

    if not skip_particles:
        coarse_particles = Particles(box=domain_box)
        upper_cell_particles = coarse_particles.select(upper_slct_box)
        lower_cell_particles = coarse_particles.select(lower_slct_box)

        coarse_particles.add(upper_cell_particles.shift_icell(-domain_box.shape()))
        coarse_particles.add(lower_cell_particles.shift_icell(domain_box.shape()))

    patch_datas = {}

    for ilvl, lvl_box in boxes.items():

        lvl_cell_width = cell_width / (refinement_ratio ** ilvl)

        if not skip_particles:
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

            datas = {qty: globals()[qty](ghost_box, layout, domain_size) for qty in quantities}

            if not skip_particles:
                datas["particles"] = lvl_particles.select(ghost_box)

            boxed_patch_datas = {}
            for qty_name, data in datas.items():
                if qty_name == 'particles':
                    pdata = ParticleData(layout, data, "pop_name")
                else:
                    pdata = FieldData(layout, qty_name, data)

                boxed_patch_datas[qty_name] = pdata

            patch_datas[ilvl].append(boxed_patch_datas)

    return patch_datas


def build_kwargs(**kwargs):

    quantities = ["Bx", "By", "Bz", "Ex", "Ey", "Ez", "particles"]

    if "simulation" in kwargs:
        for k in kwargs:
            if k != "simulation":
                print("warning: 'simulation' given, {} discarded".format(k))

        sim = kwargs["simulation"]
        kwargs["nbr_cells"] = sim.cells
        kwargs["origin"] = sim.origin
        kwargs["interp_order"] = sim.interp_order
        kwargs["domain_size"] = sim.simulation_domain()
        kwargs["cell_width"] = sim.dl
        kwargs["refinement_boxes"] = sim.refinement_boxes

    kwargs["quantities"] = listify(kwargs.get("quantities", quantities))
    kwargs["refinement_ratio"] = 2

    return kwargs


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
    kwargs = build_kwargs(**kwargs)

    nbr_cells = kwargs["nbr_cells"]
    origin = kwargs["origin"]
    interp_order = kwargs["interp_order"]
    domain_size = kwargs["domain_size"]
    cell_width = kwargs["cell_width"]
    refinement_ratio = kwargs["refinement_ratio"]

    domain_box = boxm.Box(0, nbr_cells - 1)
    boxes = build_refinementboxes(domain_box, **kwargs)
    patch_datas = build_patch_datas(domain_box, boxes, **kwargs)

    patches = {ilvl: [] for ilvl in list(patch_datas.keys())}
    for ilvl, lvl_patch_datas in patch_datas.items():
        for patch_datas in lvl_patch_datas:
            patches[ilvl].append(Patch(patch_datas))

    patch_levels = {}
    for ilvl, lvl_patches in patches.items():
        patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

    sorted_levels_numbers = sorted(patch_levels)
    patch_levels = {ilvl: patch_levels[ilvl] for ilvl in sorted_levels_numbers}
    return PatchHierarchy(patch_levels, domain_box, refinement_ratio)

