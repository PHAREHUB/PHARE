#

import numpy as np

from ..core import box as boxm
from .hierarchy.patchdata import FieldData
from .hierarchy.hierarchy_utils import is_root_lvl

from pyphare.logger import getLogger
from pyphare.core.phare_utilities import listify, is_scalar

logger = getLogger(__name__)


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


def touch_domain_border(box, domain_box, border):
    if border == "upper":
        return domain_box.upper in box
    elif border == "lower":
        return domain_box.lower in box
    else:
        raise RuntimeError("invalid border")


def is_border_patch(patch, domain_box):
    gbox = boxm.grow(patch.box, [1] * domain_box.ndim)
    return gbox * domain_box != gbox


def periodic_shifts_for(domain_box):
    import itertools

    L = np.asarray(domain_box.shape)
    ndim = domain_box.ndim

    if ndim == 1:
        return [L[0], -L[0]]

    shifts = []
    for coeffs in itertools.product([-1, 0, 1], repeat=ndim):
        if all(c == 0 for c in coeffs):
            continue  # skip zero shift
        shifts.append(np.array(coeffs) * L)

    return shifts


def border_shifted_patches_for(patch, domain_box):
    from copy import deepcopy

    ndim = domain_box.ndim
    gbox = boxm.grow(patch.box, [1] * ndim)
    if gbox * domain_box == gbox:
        return []  # not a border patch

    shifted_patches = []
    for shift in periodic_shifts_for(domain_box):
        shifted_box = boxm.shift(gbox, shift)
        if shifted_box * domain_box is not None:
            shifted_patches.append(deepcopy(patch))
            shifted_patches[-1].box = shifted_box

    return shifted_patches


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

    def check_overlap(refPatch, cmpPatch):
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

    shifts = periodic_shifts_for(domain_box)

    def check_periodic_overlap(ref_patch, cmp_patch):
        for ref_pdname, ref_pd in ref_patch.patch_datas.items():
            cmp_pd = cmp_patch.patch_datas[ref_pdname]

            gb1 = ref_pd.ghost_box
            gb2 = cmp_pd.ghost_box

            for offset in shifts:
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

    def is_border(patch):
        return is_border_patch(patch, domain_box)

    def check_overlaps(refPatch, cmpPatch):
        # no interal overlap with self
        if refPatch.box != cmpPatch.box:
            check_overlap(refPatch, cmpPatch)

        # but possible periodic overlaps
        if is_border(refPatch) and is_border(cmpPatch):
            check_periodic_overlap(refPatch, cmpPatch)

    for ip, refPatch in enumerate(patches):
        for cmpPatch in patches[ip:]:  # patches can overlap with themselves
            check_overlaps(refPatch, cmpPatch)

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

        return sorted_patches

    for patch in patches:
        sorted_patches.extend(border_shifted_patches_for(patch, domain_box))

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
                    remaining = gabox - [p.box for p in check_patches]

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
