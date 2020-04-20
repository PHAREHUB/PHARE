

from . import box as boxm
from .hierarchy import FieldData, is_root_lvl



def toFieldBox(box, patch_data):
    from core.gridlayout import GridLayout
    """
    grows the box by one cell if patch_data is primal
    because box.upper+1 allows to get the rightmost primal node
    otherwise missing
    """
    assert (isinstance(patch_data, FieldData))
    if patch_data.layout.centering["X"][patch_data.field_name] == "primal":
        return boxm.Box(box.lower, box.upper + 1)
    else:
        return box





def shift_patch(patch, offset):
    patch.box = boxm.shift(patch.box, offset)
    for pdata in patch.patch_datas.values():
        pdata.ghost_box = boxm.shift(pdata.ghost_box, offset)




def touch_domain_border(box, domain_box, border):
    if border == "upper":
        return domain_box.upper in box
    elif border == "lower":
        return domain_box.lower in box
    else:
        raise RuntimeError("invalid border")



def compute_overlaps(patches, domain_box):
    # sorting the patches per origin allows to
    # consider first and last as possible periodic overlaps

    patches = sorted(patches, key=lambda p: p.origin)
    overlaps = []
    for ip, refPatch in enumerate(patches):
        for cmpPatch in patches[ip + 1:]:
            # for two patches, compare patch_datas of the same quantity
            for ref_pdname, ref_pd in refPatch.patch_datas.items():
                for cmp_pdname, cmp_pd in cmpPatch.patch_datas.items():
                    if cmp_pdname == ref_pdname:
                        gb1 = ref_pd.ghost_box
                        gb2 = cmp_pd.ghost_box
                        overlap = gb1 * gb2
                        if overlap is not None:
                            if ref_pd.quantity == 'field':
                                overlap = toFieldBox(overlap, ref_pd)

                            overlaps.append({"pdatas": (ref_pd, cmp_pd),
                                             "box": overlap,
                                             "offset": (0, 0)})

    # now dealing with first and last patches to see their patchdata overlap
    for ref_pdname, ref_pd in patches[0].patch_datas.items():
        for cmp_pdname, cmp_pd in patches[-1].patch_datas.items():
            if cmp_pdname == ref_pdname:
                gb1 = ref_pd.ghost_box
                gb2 = cmp_pd.ghost_box

                # first check if offseting last to lower left
                # overlaps first

                offset = domain_box.upper + 1
                overlap = gb1 * boxm.shift(gb2, -offset)
                # s = "gb1 {}, gb2 {}, offset {}, shifted gb2 {} and overlap {}"
                # print(s.format(gb1, gb2, offset, shift(gb2,-offset), overlap))
                if overlap is not None:
                    if ref_pd.quantity == 'field':
                        overlap = toFieldBox(overlap, ref_pd)

                    overlaps.append({"pdatas": (ref_pd, cmp_pd),
                                     "box": overlap,
                                     "offset": (0, -offset)})

                    # if it overlaps one way, it should overlap the other
                    other_ovrlp = boxm.shift(gb1, offset) * gb2

                    assert(other_ovrlp is not None)
                    if ref_pd.quantity == 'field':
                        other_ovrlp= toFieldBox(other_ovrlp, ref_pd)


                    overlaps.append({"pdatas":(ref_pd, cmp_pd),
                                     "box" : other_ovrlp,
                                     "offset":offset})


    return overlaps




def hierarchy_overlaps(hierarchy):
    overlaps = {}
    for ilvl, lvl in enumerate(hierarchy.patch_levels):
        overlaps[ilvl] = compute_overlaps(lvl.patches, hierarchy.refined_domain_box(ilvl))
    return overlaps



def get_periodic_list(patches, domain_box):
    """
    given a list of patches and a domain box the function
    returns a list of patches sorted by origin where the
    first (resp. last) patch is the last (resp. first) patch
    of the sorted list, if that patch touches the upper (resp. lower)
    border of the domain
    """
    from copy import copy

    sorted_patches = sorted(patches, key=lambda p: p.origin)

    if touch_domain_border(sorted_patches[-1].box, domain_box, "upper"):
        last_patch = copy(sorted_patches[-1])
        shift_patch(last_patch, -domain_box.upper)
        sorted_patches.insert(0, last_patch)

    if touch_domain_border(sorted_patches[0].box, domain_box, "lower"):
        first_patch = copy(sorted_patches[0])
        shift_patch(first_patch, domain_box.upper)
        sorted_patches.append(first_patch)

    return sorted_patches






def particle_ghost_area_boxes(hierarchy):
    """
    this function returns boxes representing ghost cells for all levels
    return : {level_number : [{"pdata":patch_data1, "boxes":ghost_boxes},
                              {"pdata":patch_data2, "boxes":ghost_boxes}, ...]}
    """
    gaboxes = {}

    for ilvl, lvl in enumerate(hierarchy.patch_levels):
        for patch in lvl.patches:

            patch_data = patch.patch_datas["particles"]
            gbox = patch_data.ghost_box
            box = patch.box

            if ilvl not in gaboxes:
                gaboxes[ilvl] = []

            gaboxes[ilvl] += [{"pdata": patch_data, "boxes": boxm.remove(gbox, box)}]

    return gaboxes


def level_ghost_boxes(hierarchy):
    """
    this function returns boxes representing level ghost cells for all levels
    return : {level_number : [{"pdata":patch_data1, "boxes":lvl_ghost_boxes},
                              {"pdata":patch_data2, "boxes":lvl_ghost_boxes}, ...]}
    """
    gaboxes = particle_ghost_area_boxes(hierarchy)
    lvl_gaboxes = {}

    for ilvl, lvl in enumerate(hierarchy.patch_levels):

        if not is_root_lvl(lvl):  # level ghost do not make sense for periodic root level

            gaboxes_info = gaboxes[ilvl]

            for info in gaboxes_info:

                patch_data = info["pdata"]
                gaboxes = info["boxes"]

                for gabox in gaboxes:

                    # print("gabox: {}".format(gabox))

                    # now loop on all particle patchData
                    # keep only parts of the ghost boxes that do
                    # not intersect other patch data interior

                    if True:  # if periodic, always true for now
                        refined_domain_box = hierarchy.refined_domain_box(ilvl)
                        patches = get_periodic_list(lvl.patches, refined_domain_box)

                    for patch in patches:

                        if patch.patch_datas["particles"] is not patch_data:

                            keep = boxm.remove(gabox, patch.box)

                            # print("{}*{} = {}, keeping {} boxes".format(gabox, patch.box,
                            #                                           gabox*patch.box,
                            #                                          len(keep)))
                            # for k in keep:
                            #    print("keep : {}".format(k))

                            if ilvl not in lvl_gaboxes:
                                lvl_gaboxes[ilvl] = []

                            if len(keep):
                                lvl_gaboxes[ilvl] += [{"pdata": patch_data, "boxes": keep}]

    return lvl_gaboxes

