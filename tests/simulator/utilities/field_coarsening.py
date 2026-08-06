import numpy as np
from pyphare.core import gridlayout
from pyphare.core.gridlayout import directions
from pyphare.core.phare_utilities import refinement_ratio


def coarsen(qty, coarseField, fineField, coarseBox, fineData, coarseData):
    coarseLayout = coarseField.layout
    fineLayout = fineField.layout
    ndim = coarseLayout.box.ndim

    is_primal = []
    for direction in directions[:ndim]:
        is_primal += [gridlayout.yee_element_is_primal(qty, direction)]

    def coarse_indices(dim):
        return np.arange(coarseBox.lower[dim], coarseBox.upper[dim] + 1)

    def fineLocal(coarse_index, dim):
        return fineLayout.AMRIndexToLocal(dim, coarse_index * refinement_ratio)

    def coarseLocal(index, dim):
        return coarseLayout.AMRIndexToLocal(dim, index)

    if ndim == 1:
        for index in coarse_indices(0):
            fineIndex = fineLocal(index, 0)
            coarseLocalIndex = coarseLocal(index, 0)
            if is_primal[0]:
                coarseData[coarseLocalIndex] = fineData[fineIndex]
            else:
                # ind 1D it is the same formula for By and Bz as for other quantities
                coarseData[coarseLocalIndex] = (
                    fineData[fineIndex] * 0.5 + fineData[fineIndex + 1] * 0.5
                )

    if ndim == 2:
        for indexX in coarse_indices(0):
            fineIndexX = fineLocal(indexX, 0)
            coarseLocalIndexX = coarseLocal(indexX, 0)

            for indexY in coarse_indices(1):
                fineIndexY = fineLocal(indexY, 1)
                coarseLocalIndexY = coarseLocal(indexY, 1)
                left, middle, right = 0, 0, 0
                if all(is_primal):
                    coarseData[coarseLocalIndexX][coarseLocalIndexY] = fineData[
                        fineIndexX
                    ][fineIndexY]

                if is_primal[0] and not is_primal[1]:
                    coarseData[coarseLocalIndexX, coarseLocalIndexY] = 0.5 * (
                        fineData[fineIndexX, fineIndexY]
                        + fineData[fineIndexX, fineIndexY + 1]
                    )

                if not is_primal[0] and is_primal[1]:
                    coarseData[coarseLocalIndexX, coarseLocalIndexY] = 0.5 * (
                        fineData[fineIndexX, fineIndexY]
                        + fineData[fineIndexX + 1, fineIndexY]
                    )

                if not any(is_primal):
                    coarseData[coarseLocalIndexX, coarseLocalIndexY] = 0.25 * (
                        fineData[fineIndexX, fineIndexY]
                        + fineData[fineIndexX, fineIndexY + 1]
                        + fineData[fineIndexX + 1, fineIndexY + 1]
                        + fineData[fineIndexX + 1, fineIndexY]
                    )
