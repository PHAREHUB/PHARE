from phare.pp.diagnostics import Patch, PatchLevel
from .periodic_overlap import intralevel as periodic_intralevel
from typing import Type
from abc import ABC  # abstract class

from phare.data.wrangler import safe_round


class Overlap(ABC):
    """
    Abstract class for different types of quantity overlaps

    Parameters:
    -------------------

    patch0    : Patch object with ghost box overlap of patch1
    patch1    : Patch object with ghost box overlap of patch0
    data_name : dataset key for patch.data(key) pertaining to this overlap
    nGhosts   : number of ghosts considered in the overlap
    sizes     : size of the overlap box

    """

    def __init__(self, patch0, patch1, data_name, nGhosts, sizes):
        self.patch0, self.patch1 = (patch0, patch1)
        self.data_name = data_name
        self.nGhosts = nGhosts
        self.sizes = sizes


def getOverlaps(OverlapType, diags):
    assert issubclass(OverlapType, Overlap) and isinstance(diags, dict)

    if not OverlapType.patch_data_type.__name__ in diags:
        return []

    diags = diags[OverlapType.patch_data_type.__name__]

    return _intralevel(OverlapType, diags) + _periodic(OverlapType, diags)


def calculate_1d(OverlapType, patch0, patch1, data_name, **kwargs):
    assert issubclass(OverlapType, Overlap)
    assert type(patch0.patch_data) == type(patch1.patch_data)
    assert patch0.patch_level.lvlNbr == patch1.patch_level.lvlNbr

    lower, upper = sorted([patch0, patch1], key=lambda x: x.min_coord("x"))
    nGhosts = patch0.patch_data.nGhosts(data_name)
    level_cell_width = patch0.patch_level.cell_width("x")
    x_gap = upper.min_coord("x") - lower.max_coord("x")

    assert x_gap >= 0  # should never happen

    ghost_area_lgth = level_cell_width * nGhosts

    if x_gap < ghost_area_lgth:

        cell_gap = safe_round(x_gap, level_cell_width)
        return [
            OverlapType(
                lower, upper, data_name, nGhosts, [int(nGhosts - cell_gap)], **kwargs
            )
        ]

    return []


def _intralevel(OverlapType, diags):
    assert issubclass(OverlapType, Overlap)

    overlaps = []
    for diag in diags:
        for lvlNbr, patch_level in diag.levels.items():
            patches = patch_level.patches_list()
            for i, refPatch in enumerate(patches):
                for cmpPatch in patches[i + 1 :]:
                    for data_name in refPatch.patch_data.dataset_names():
                        overlaps += calculate_1d(
                            OverlapType, refPatch, cmpPatch, data_name
                        )

    return overlaps


def _periodic(OverlapType, diags):
    assert issubclass(OverlapType, Overlap)

    overlaps = []
    for diag in diags:
        for lvlNbr, patch_level in diag.levels.items():
            for data_name, border_patches in periodic_intralevel(
                diag, patch_level
            ).items():
                minX, maxX = border_patches
                if len(minX) and len(maxX):
                    overlaps += calculate_1d(OverlapType, minX[0], maxX[0], data_name)

    return overlaps
