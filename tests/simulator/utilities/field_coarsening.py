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
                if qty == "Bx" or qty == "Ey" or qty == "Ez":
                    coarseData[coarseLocalIndex] = fineData[fineIndex]
                else:
                    coarseData[coarseLocalIndex] = (
                        fineData[fineIndex - 1] * 0.25
                        + fineData[fineIndex] * 0.5
                        + fineData[fineIndex + 1] * 0.25
                    )
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
                    if qty == "Ez":
                        coarseData[coarseLocalIndexX][coarseLocalIndexY] = fineData[
                            fineIndexX
                        ][fineIndexY]
                    else:
                        left += fineData[fineIndexX - 1][fineIndexY - 1] * 0.25
                        left += fineData[fineIndexX - 1][fineIndexY] * 0.5
                        left += fineData[fineIndexX - 1][fineIndexY + 1] * 0.25
                        middle += fineData[fineIndexX][fineIndexY - 1] * 0.25
                        middle += fineData[fineIndexX][fineIndexY] * 0.5
                        middle += fineData[fineIndexX][fineIndexY + 1] * 0.25
                        right += fineData[fineIndexX + 1][fineIndexY - 1] * 0.25
                        right += fineData[fineIndexX + 1][fineIndexY] * 0.5
                        right += fineData[fineIndexX + 1][fineIndexY + 1] * 0.25
                        coarseData[coarseLocalIndexX][coarseLocalIndexY] = (
                            left * 0.25 + middle * 0.5 + right * 0.25
                        )

                if is_primal[0] and not is_primal[1]:
                    if qty == "Bx" or qty == "Ey":
                        coarseData[coarseLocalIndexX, coarseLocalIndexY] = 0.5 * (
                            fineData[fineIndexX, fineIndexY]
                            + fineData[fineIndexX, fineIndexY + 1]
                        )
                    else:
                        left += fineData[fineIndexX - 1][fineIndexY] * 0.5
                        left += fineData[fineIndexX - 1][fineIndexY + 1] * 0.5
                        middle += fineData[fineIndexX][fineIndexY] * 0.5
                        middle += fineData[fineIndexX][fineIndexY + 1] * 0.5
                        right += fineData[fineIndexX + 1][fineIndexY] * 0.5
                        right += fineData[fineIndexX + 1][fineIndexY + 1] * 0.5
                        coarseData[coarseLocalIndexX][coarseLocalIndexY] = (
                            left * 0.25 + middle * 0.5 + right * 0.25
                        )

                if not is_primal[0] and is_primal[1]:
                    if qty == "By" or qty == "Ex":
                        coarseData[coarseLocalIndexX, coarseLocalIndexY] = 0.5 * (
                            fineData[fineIndexX, fineIndexY]
                            + fineData[fineIndexX + 1, fineIndexY]
                        )

                    else:
                        left += fineData[fineIndexX][fineIndexY - 1] * 0.25
                        left += fineData[fineIndexX][fineIndexY] * 0.5
                        left += fineData[fineIndexX][fineIndexY + 1] * 0.25
                        right += fineData[fineIndexX + 1][fineIndexY - 1] * 0.25
                        right += fineData[fineIndexX + 1][fineIndexY] * 0.5
                        right += fineData[fineIndexX + 1][fineIndexY + 1] * 0.25
                        coarseData[coarseLocalIndexX][coarseLocalIndexY] = (
                            left * 0.5 + right * 0.5
                        )

                if not any(is_primal):
                    left += fineData[fineIndexX][fineIndexY] * 0.5
                    left += fineData[fineIndexX][fineIndexY + 1] * 0.5
                    right += fineData[fineIndexX + 1][fineIndexY] * 0.5
                    right += fineData[fineIndexX + 1][fineIndexY + 1] * 0.5
                    coarseData[coarseLocalIndexX][coarseLocalIndexY] = (
                        left * 0.5 + right * 0.5
                    )
