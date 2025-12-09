import numpy as np

from ..core import box as boxm
from pyphare.core.box import Box
from .hierarchy.patchdata import FieldData
from .hierarchy.hierarchy_utils import is_root_lvl

from pyphare.core.phare_utilities import listify, is_scalar


def toFieldBox(box, patch_data):
    """
    grows the box by one cell if patch_data is primal
    because box.upper+1 allows to get the rightmost primal node
    otherwise missing
    """
    assert isinstance(patch_data, FieldData)

    from copy import copy

    box = copy(box)

    directions = ["X", "Y", "Z"][: box.ndim]  # drop unused directions
    for i, direction in enumerate(directions):
        if patch_data.layout.centering[direction][patch_data.field_name] == "primal":
            box.upper[i] = box.upper[i] + 1

    return box


def shift_patch(patch, offset):
    patch.box = boxm.shift(patch.box, offset)
    for pdata in patch.patch_datas.values():
        pdata.box = boxm.shift(pdata.box, offset)
        pdata.ghost_box = boxm.shift(pdata.ghost_box, offset)


def domain_border_ghost_boxes(domain_box, patches):
    max_ghosts = max(
        [
            pd.ghosts_nbr.max()
            for pd in sum([list(patch.patch_datas.values()) for patch in patches], [])
        ]
    )
    ghost_box_width = max_ghosts - 1

    if domain_box.ndim == 1:
        upper_x = domain_box.upper
        return {
            "left": Box(0, ghost_box_width),
            "right": Box(upper_x - ghost_box_width, upper_x),
        }

    elif domain_box.ndim == 2:
        upper_x, upper_y = domain_box.upper
        return {
            "bottom": Box(
                (
                    0,
                    0,
                ),
                (upper_x, ghost_box_width),
            ),
            "top": Box((0, upper_y - ghost_box_width), (upper_x, upper_y)),
            "left": Box((0, 0), (ghost_box_width, upper_y)),
            "right": Box((upper_x - ghost_box_width, 0), (upper_x, upper_y)),
        }

    raise ValueError("Unhandeled dimension")


def touch_domain_border(box, domain_box, border):
    if border == "upper":
        return domain_box.upper in box
    elif border == "lower":
        return domain_box.lower in box
    else:
        raise RuntimeError("invalid border")


def periodicity_shifts(domain_box):
    if domain_box.ndim == 1:
        shape_x = domain_box.shape
        return {
            "left": shape_x,
            "right": -shape_x,
        }

    if domain_box.ndim == 2:
        shape_x, shape_y = domain_box.shape
        shifts = {
            "left": [(shape_x, 0)],
            "right": [(-shape_x, 0)],
            "bottom": [(0, shape_y)],
            "top": [(0, -shape_y)],
        }
        shifts.update(
            {
                "bottomleft": [*shifts["left"], *shifts["bottom"], (shape_x, shape_y)],
                "bottomright": [
                    *shifts["right"],
                    *shifts["bottom"],
                    (-shape_x, shape_y),
                ],
                "topleft": [*shifts["left"], *shifts["top"], (shape_x, -shape_y)],
                "topright": [*shifts["right"], *shifts["top"], (-shape_x, -shape_y)],
                "bottomtop": [*shifts["bottom"], *shifts["top"]],
                "leftright": [*shifts["left"], *shifts["right"]],
            }
        )
        shifts.update(
            {
                "bottomtopleft": [
                    *shifts["bottomtop"],
                    *shifts["left"],
                    shifts["bottomleft"][-1],
                    shifts["topleft"][-1],
                ],
                "bottomtopright": [
                    *shifts["bottomtop"],
                    *shifts["right"],
                    shifts["bottomright"][-1],
                    shifts["topright"][-1],
                ],
                "bottomleftright": [
                    *shifts["leftright"],
                    *shifts["bottom"],
                    shifts["bottomleft"][-1],
                    shifts["bottomright"][-1],
                ],
                "topleftright": [
                    *shifts["leftright"],
                    *shifts["top"],
                    shifts["topleft"][-1],
                    shifts["topright"][-1],
                ],
                "bottomtopleftright": [  # one patch covers domain
                    *shifts["bottomleft"],
                    *shifts["topright"],
                    shifts["bottomright"][-1],
                    shifts["topleft"][-1],
                ],
            }
        )

    if domain_box.ndim == 3:
        raise ValueError("Unhandeled dimension")

    return shifts


def compute_overlaps(patches, domain_box):
    """
    returns a list of overlaps for all patch datas in given patches
    and for a domain box. An overlap is defined as an intersection of
    two PatchData ghost box.

    if  - = domain cell, o = ghost cell

    PatchData1  : - - - - - - - - - o o o o o
    PatchData2  :                     o o o o o - - - - - - - -
    overlap box :                     ^ ^ ^ ^

    for now the function is purely 1D and assumes the direction is
    periodic. Therefore, PatchData close to the lower (resp. upper)
    domain boundary (hence the need of domain_box) may overlap another
    one close to the upper (resp. lower) end of the domain.

    an overlap is a dictionnary with 3 keys:
    - pdatas : value is a tuple of the two overlapped PatchDatas.
               this is useful when one wants to get overlapped data
    - box    :  is the intersection box between the two PatchData ghost boxes
    - offset : tuple of two offsets by which the patchData ghost box is shifted
               to compute the overlap box.
    """

    # sorting the patches per origin allows to
    # consider first and last as possible periodic overlaps
    dim = domain_box.ndim

    if dim == 1:
        patches = sorted(patches, key=lambda p: p.origin.all())

    zero_offset = [0] * dim if dim > 1 else 0

    overlaps = []

    # first deal with intra domain overlaps
    for ip, refPatch in enumerate(patches):
        for cmpPatch in patches[ip + 1 :]:
            # for two patches, compare patch_datas of the same quantity
            for ref_pdname, ref_pd in refPatch.patch_datas.items():
                cmp_pd = cmpPatch.patch_datas[ref_pdname]

                gb1 = ref_pd.ghost_box
                gb2 = cmp_pd.ghost_box
                overlap = gb1 * gb2

                if overlap is not None:
                    # boxes indexes represent cells
                    # therefore for fields, we need to
                    # adjust the box. This essentially
                    # add 1 to upper in case field is on corners
                    # because upper corner can only be grabbed
                    # if box extends to upper+1 cell
                    # particles don't need that as they are contained
                    # in cells.
                    if ref_pd.quantity == "field":
                        overlap = toFieldBox(overlap, ref_pd)

                    overlaps.append(
                        {
                            "name": ref_pdname,
                            "pdatas": (ref_pd, cmp_pd),
                            "patches": (refPatch, cmpPatch),
                            "box": overlap,
                            "offset": (zero_offset, zero_offset),
                        }
                    )

    def append(name, ref_pd, cmp_pd, refPatch, cmpPatch, overlap, offset_tuple):
        overlaps.append(
            {
                "name": name,
                "pdatas": (ref_pd, cmp_pd),
                "patches": (refPatch, cmpPatch),
                "box": overlap,
                "offset": offset_tuple,
            }
        )

    def borders_per(patch):
        return "".join(
            [key for key, side in sides.items() if patch.box * side is not None]
        )

    sides = domain_border_ghost_boxes(domain_box, patches)
    shifts = periodicity_shifts(domain_box)

    # now dealing with border patches to see their patchdata overlap
    if dim == 1 and len(patches) > 0:
        patches = [patches[0], patches[-1]]

    # filter out patches not near a border
    borders_per_patch = {p: borders_per(p) for p in patches}
    border_patches = [
        p for p, in_sides in borders_per_patch.items() if len(in_sides) > 0
    ]

    for patch_i, ref_patch in enumerate(border_patches):
        in_sides = borders_per_patch[ref_patch]
        assert in_sides in shifts

        for ref_pdname, ref_pd in ref_patch.patch_datas.items():
            for shift in shifts[in_sides]:
                for cmp_patch in border_patches[
                    patch_i:
                ]:  # patches can overlap with themselves
                    for cmp_pdname, cmp_pd in cmp_patch.patch_datas.items():
                        if cmp_pdname == ref_pdname:
                            gb1 = ref_pd.ghost_box
                            gb2 = cmp_pd.ghost_box

                            offset = np.asarray(shift)
                            overlap = gb1 * boxm.shift(gb2, -offset)

                            if overlap is not None:
                                other_ovrlp = boxm.shift(gb1, offset) * gb2
                                assert other_ovrlp is not None

                                if ref_pd.quantity == "field":
                                    overlap = toFieldBox(overlap, ref_pd)
                                    other_ovrlp = toFieldBox(other_ovrlp, ref_pd)

                                append(
                                    ref_pdname,
                                    ref_pd,
                                    cmp_pd,
                                    ref_patch,
                                    cmp_patch,
                                    overlap,
                                    (zero_offset, (-offset).tolist()),
                                )
                                append(
                                    ref_pdname,
                                    ref_pd,
                                    cmp_pd,
                                    ref_patch,
                                    cmp_patch,
                                    other_ovrlp,
                                    (offset.tolist(), zero_offset),
                                )

    return overlaps


def hierarchy_overlaps(hierarchy, time=0):
    """
    returns all overlaps for the given hierarchy
    """
    overlaps = {}

    for ilvl, lvl in hierarchy.levels(time).items():
        overlaps[ilvl] = compute_overlaps(
            lvl.patches, hierarchy.refined_domain_box(ilvl)
        )
    return overlaps


def get_periodic_list(patches, domain_box, n_ghosts):
    """
    given a list of patches and a domain box the function
    returns a list of patches sorted by origin where the
    first (resp. last) patch is the last (resp. first) patch
    of the sorted list, if that patch touches the upper (resp. lower)
    border of the domain
    """
    assert len(patches) > 0
    dim = patches[-1].box.ndim
    assert all([p.box.ndim == dim for p in patches])

    from copy import copy

    sorted_patches = sorted(patches, key=lambda p: p.origin.all())

    if dim == 1:  # only two possible border patches, [0] and [-1]
        # copy before check as the list is modified in-place
        last_patch = copy(sorted_patches[-1])
        first_patch = copy(sorted_patches[0])

        if touch_domain_border(
            boxm.grow(sorted_patches[-1].box, n_ghosts), domain_box, "upper"
        ):
            shift_patch(last_patch, -domain_box.shape)
            sorted_patches.insert(0, last_patch)

        if touch_domain_border(
            boxm.grow(sorted_patches[0].box, n_ghosts), domain_box, "lower"
        ):
            shift_patch(first_patch, domain_box.shape)
            sorted_patches.append(first_patch)

    if dim == 2:
        sides = {
            "bottom": Box([0, 0], [domain_box.upper[0], 0]),
            "top": Box(
                [0, domain_box.upper[1]], [domain_box.upper[0], domain_box.upper[1]]
            ),
            "left": Box([0, 0], [0, domain_box.upper[1]]),
            "right": Box(
                [domain_box.upper[0], 0], [domain_box.upper[0], domain_box.upper[1]]
            ),
        }

        shifts = periodicity_shifts(domain_box)

        def borders_per(box):
            return "".join(
                [key for key, side in sides.items() if box * side is not None]
            )

        for patch in patches:
            in_sides = borders_per(boxm.grow(patch.box, n_ghosts))

            if in_sides in shifts:  # in_sides might be empty, so no borders
                for shift in shifts[in_sides]:
                    patch_copy = copy(patch)
                    shift_patch(patch_copy, shift)
                    sorted_patches.append(patch_copy)

    if dim == 3:
        raise ValueError("not yet implemented")

    return sorted_patches


def ghost_area_boxes(hierarchy, quantities, levelNbrs=[], time=0):
    """
    this function returns boxes representing ghost cell boxes for all levels
    a ghost cell box is a box containing cells of contiguous AMR index not
    contained in the domain box.

    if  - = domain cell and o = ghost cell

    patchdata = o o o o o  - - - - - - - - - - o o o o o
    boxes =     ^-------^                      ^-------^

    return : {level_number : [{"pdata":patch_data1, "boxes":ghost_boxes},
                              {"pdata":patch_data2, "boxes":ghost_boxes}, ...]}
    """
    levelNbrs = listify(levelNbrs)
    if len(levelNbrs) == 0:
        levelNbrs = list(hierarchy.levels(time).keys())

    gaboxes = {}

    for ilvl in levelNbrs:
        lvl = hierarchy.level(ilvl, time)
        for patch in lvl.patches:
            for pd_key, pd in patch.patch_datas.items():
                skip = not any([pd_key.endswith(qty) for qty in quantities])

                if skip:
                    continue

                patch_data = patch.patch_datas[pd_key]
                gbox = patch_data.ghost_box
                box = patch.box

                if ilvl not in gaboxes:
                    gaboxes[ilvl] = {}

                if pd_key not in gaboxes[ilvl]:
                    gaboxes[ilvl][pd_key] = []

                gaboxes[ilvl][pd_key] += [
                    {"pdata": patch_data, "boxes": boxm.remove(gbox, box)}
                ]

    return gaboxes


def level_ghost_boxes(hierarchy, quantities, levelNbrs=[], time=None):
    """
    this function returns boxes representing level ghost cell boxes for all levels
    A level ghost cell box is a ghost cell box that does not overlap any cell contained
    in a patchData interior
    patchdata1           : o o o o o - - - - - - - - - - - - o o o o o
    patchdata2           :                               o o o o o - - - - - - - - -
    lvl ghost cell boxes :                                   ^---^
    returns a dictionnary which keys are level_number and value is a list of dict with :
     keys:value :
        - pdata : patch_data for which level ghost cell boxes are detected
        - boxes : level ghost cell boxes
    return : {level_number : [{"pdata":patch_data1, "boxes":lvl_ghost_boxes},
                              {"pdata":patch_data2, "boxes":lvl_ghost_boxes}, ...]}

    Other parameters
    ----------------
      levelNbrs : limit working set of hierarchy levels to those requested, if scalar, returns just that level
      time      : the simulation time to access the appropriate data for the requested time
    """
    quantities = listify(quantities)

    levelNbrs_is_scalar = is_scalar(levelNbrs)
    levelNbrs = listify(levelNbrs)
    if len(levelNbrs) == 0:
        levelNbrs = list(hierarchy.levels(time).keys())

    gaboxes = ghost_area_boxes(hierarchy, quantities, levelNbrs, time)
    lvl_gaboxes = {}

    for ilvl in levelNbrs:
        lvl = hierarchy.level(ilvl, time)

        if is_root_lvl(lvl):  # level ghost do not make sense for periodic root level
            continue

        for pd_key, info_list in gaboxes[ilvl].items():
            # if periodic, always true for now
            if True:  # lgtm [py/constant-conditional-expression]
                refined_domain_box = hierarchy.refined_domain_box(ilvl)
                n_ghosts = lvl.patches[0].patch_datas[pd_key].ghosts_nbr
                patches = get_periodic_list(lvl.patches, refined_domain_box, n_ghosts)

            for info in info_list:
                patch_data, ghostAreaBoxes = info["pdata"], info["boxes"]

                check_patches = [
                    p for p in patches if p.patch_datas[pd_key] is not patch_data
                ]

                if len(check_patches) == 0:
                    check_patches = patches

                for gabox in ghostAreaBoxes:
                    remaining = gabox - check_patches[0].box

                    for patch in check_patches[1:]:
                        tmp = []
                        remove = []
                        for i, rem in enumerate(remaining):
                            if rem * patch.box is not None:
                                remove.append(i)
                                tmp += rem - patch.box
                        for rm in reversed(remove):
                            del remaining[rm]
                        remaining += tmp

                    if ilvl not in lvl_gaboxes:
                        lvl_gaboxes[ilvl] = {}

                    if pd_key not in lvl_gaboxes[ilvl]:
                        lvl_gaboxes[ilvl][pd_key] = []

                    if len(remaining):
                        lvl_gaboxes[ilvl][pd_key] += [
                            {"pdata": patch_data, "boxes": remaining}
                        ]

    if levelNbrs_is_scalar:
        return lvl_gaboxes[levelNbrs[0]]
    return lvl_gaboxes
