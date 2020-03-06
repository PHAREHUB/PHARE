import math
from phare.pp.diagnostics import Diagnostic, Patch, PatchLevel
from phare.pp.diagnostic.dtype import _DType
from .periodic_overlap import PeriodicOverlap
from typing import Type


class Overlap:
    def __init__(self, p0: Patch, p1: Patch, key: str, nGhosts: int, offsets, sizes):
        self.p0, self.p1 = (p0, p1)
        self.key = key  # data dictionary key for patches
        self.nGhosts = nGhosts
        self.offsets = offsets  # list offset per dim
        self.sizes = sizes  # list size per dim

    @staticmethod
    def _float_equal(f0, f1):
        return math.isclose(f0, f1, rel_tol=1e-6)  # acceptable?

    @staticmethod
    def _get(diags, clazz):
        assert issubclass(clazz, Overlap)

        if isinstance(diags, dict) and clazz.dType.__name__ in diags:
            diags = diags[clazz.dType.__name__]

        return Overlap._intralevel(diags, clazz) + Overlap._periodic(diags, clazz)

    @staticmethod
    def _intralevel(diags, clazz):
        assert issubclass(clazz, Overlap)

        overlaps = []
        for diag in diags:
            for patch_level_idx, patch_level in diag.levels.items():
                patches = patch_level.patches
                for i, patch0 in enumerate(patches):
                    for j in range(i + 1, len(patches)):
                        for dataset_key in patch0.dtype.keys():
                            overlaps += Overlap._calculate_1d(
                                patch0, patches[j], dataset_key, clazz
                            )
        return overlaps

    @staticmethod
    def _periodic(diags, clazz):
        assert issubclass(clazz, Overlap)

        overlaps = []
        for diag in diags:
            for patch_level_idx, patch_level in diag.levels.items():
                for dataset_key, border_patches in PeriodicOverlap.intralevel(
                    diag, patch_level.patches
                ).items():
                    minX, maxX = border_patches
                    if len(minX) and len(maxX):
                        overlaps += Overlap._calculate_1d(
                            minX[0], maxX[0], dataset_key, clazz
                        )
        return overlaps

    @staticmethod
    def _calculate_1d(patch0, patch1, dataset_key: str, clazz):
        assert patch0.patch_level.idx == patch1.patch_level.idx
        assert issubclass(clazz, Overlap)

        lower, upper = sorted([patch0, patch1], key=lambda x: x.min_coord("x"))
        nGhosts = patch0.dtype.nGhosts(dataset_key)
        level_cell_width = patch0.patch_level.cell_width("x")
        distance = upper.min_coord("x") - lower.max_coord("x")

        if distance < (level_cell_width * nGhosts):
            cell_gap = round(distance / level_cell_width)

            return [
                clazz(
                    lower,
                    upper,
                    dataset_key,
                    nGhosts,
                    [(lower.cells[0] + cell_gap, 0)],
                    [int(nGhosts - cell_gap)],
                )
            ]
        return []
