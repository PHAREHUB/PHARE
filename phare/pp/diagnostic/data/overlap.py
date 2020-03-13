import math
from phare.pp.diagnostics import Diagnostic, Patch, PatchLevel
from phare.pp.diagnostic.dtype import _DType
from .periodic_overlap import intralevel as periodic_intralevel
from typing import Type


class Overlap:
    def __init__(self, p0: Patch, p1: Patch, key: str, nGhosts: int, sizes):
        self.p0, self.p1 = (p0, p1)
        self.key = key  # data dictionary key for patches
        self.nGhosts = nGhosts
        self.sizes = sizes  # list size per dim


def getOverlaps(clazz, diags):
    assert issubclass(clazz, Overlap)

    if isinstance(diags, dict) and clazz.dType.__name__ in diags:
        diags = diags[clazz.dType.__name__]

    return _intralevel(clazz, diags) + _periodic(clazz, diags)


def _calculate_1d(clazz, patch0, patch1, dataset_key):
    assert patch0.patch_level.idx == patch1.patch_level.idx
    assert issubclass(clazz, Overlap)

    lower, upper = sorted([patch0, patch1], key=lambda x: x.min_coord("x"))
    nGhosts = patch0.dtype.nGhosts(dataset_key)
    level_cell_width = patch0.patch_level.cell_width("x")
    x_gap = upper.min_coord("x") - lower.max_coord("x")

    if x_gap < (level_cell_width * nGhosts):
        cell_gap = round(x_gap / level_cell_width)
        return [clazz(lower, upper, dataset_key, nGhosts, [int(nGhosts - cell_gap)])]
    return []


def _intralevel(clazz, diags):
    assert issubclass(clazz, Overlap)

    overlaps = []
    for diag in diags:
        for patch_level_idx, patch_level in diag.levels.items():
            patches = patch_level.patches
            for i, patch0 in enumerate(patches):
                for j in range(i + 1, len(patches)):
                    for dataset_key in patch0.dtype.keys():
                        overlaps += _calculate_1d(
                            clazz, patch0, patches[j], dataset_key
                        )
    return overlaps


def _periodic(clazz, diags):
    assert issubclass(clazz, Overlap)

    overlaps = []
    for diag in diags:
        for patch_level_idx, patch_level in diag.levels.items():
            for dataset_key, border_patches in periodic_intralevel(
                clazz, diag, patch_level.patches
            ).items():
                minX, maxX = border_patches
                if len(minX) and len(maxX):
                    overlaps += _calculate_1d(clazz, minX[0], maxX[0], dataset_key)
    return overlaps
