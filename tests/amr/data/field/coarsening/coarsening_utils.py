from pyphare.core import gridlayout
from pyphare.core.gridlayout import directions
from pyphare.core.phare_utilities import refinement_ratio

import numpy as np


def coarsen(qty, coarseField, fineField, coarseBox, fineData, coarseData):
    """
    this function takes data in fineData (the dataset from the fineFied FieldData object)
    and put it in the coarseData dataset

    The coarsening method is the equivalent (and union) of the C++ coarsening methods,
    namely:
    - DefaultFieldCoarsener, used for all fields (moments and E) but the magnetic field B
    - MagneticFieldCoarsener, used for the magnetic field B

    as in the C++ counterpart, we hard-code here that the layout of B follows the yee layout.
    That is Bx is PDD, By is DPD and Bz is DDP
    """
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
                if qty == "Bx":
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
                    if qty == "Bx":
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
                    if qty == "By":
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
                    if qty == "Bz":
                        coarseData[coarseLocalIndexX, coarseLocalIndexY] = 0.25 * (
                            fineData[fineIndexX, fineIndexY]
                            + fineData[fineIndexX, fineIndexY + 1]
                            + fineData[fineIndexX + 1, fineIndexY + 1]
                            + fineData[fineIndexX + 1, fineIndexY]
                        )
                    else:
                        left += fineData[fineIndexX][fineIndexY] * 0.5
                        left += fineData[fineIndexX][fineIndexY + 1] * 0.5
                        right += fineData[fineIndexX + 1][fineIndexY] * 0.5
                        right += fineData[fineIndexX + 1][fineIndexY + 1] * 0.5
                        coarseData[coarseLocalIndexX][coarseLocalIndexY] = (
                            left * 0.5 + right * 0.5
                        )
    if ndim == 3:
        for indexX in coarse_indices(0):
            fineIndexX = fineLocal(indexX, 0)
            coarseLocalIndexX = coarseLocal(indexX, 0)

            for indexY in coarse_indices(1):
                fineIndexY = fineLocal(indexY, 1)
                coarseLocalIndexY = coarseLocal(indexY, 1)

                for indexZ in coarse_indices(2):
                    fineIndexZ = fineLocal(indexZ, 2)
                    coarseLocalIndexZ = coarseLocal(indexZ, 2)

                    # all primal is for moments
                    if all(is_primal):
                        value = 0.0
                        ishift = (-1, 0, 1)
                        iweights = (0.25, 0.5, 0.25)
                        jshift = (-1, 0, 1)
                        jweights = (0.25, 0.5, 0.25)
                        kshift = (-1, 0, 1)
                        kweights = (0.25, 0.5, 0.25)
                        for i, iweight in zip(ishift, iweights):
                            yvalue = 0.0
                            for j, jweight in zip(jshift, jweights):
                                zvalue = 0.0
                                for k, kweight in zip(kshift, kweights):
                                    zvalue += (
                                        fineData[fineIndexX + i][fineIndexY + j][
                                            fineIndexZ + k
                                        ]
                                        * kweight
                                    )
                                yvalue += zvalue * jweight
                            value += yvalue * iweight
                        coarseData[
                            coarseLocalIndexX, coarseLocalIndexY, coarseLocalIndexZ
                        ] = value

                    if is_primal[0] and not is_primal[1] and not is_primal[2]:
                        if qty == "Bx":
                            coarseData[
                                coarseLocalIndexX, coarseLocalIndexY, coarseLocalIndexZ
                            ] = 0.25 * (
                                fineData[fineIndexX, fineIndexY, fineIndexZ]
                                + fineData[fineIndexX, fineIndexY + 1, fineIndexZ]
                                + fineData[fineIndexX, fineIndexY, fineIndexZ + 1]
                                + fineData[fineIndexX, fineIndexY + 1, fineIndexZ + 1]
                            )

                    elif not is_primal[0] and is_primal[1] and not is_primal[2]:
                        if qty == "By":
                            coarseData[
                                coarseLocalIndexX, coarseLocalIndexY, coarseLocalIndexZ
                            ] = 0.25 * (
                                fineData[fineIndexX, fineIndexY, fineIndexZ]
                                + fineData[fineIndexX + 1, fineIndexY, fineIndexZ]
                                + fineData[fineIndexX, fineIndexY, fineIndexZ + 1]
                                + fineData[fineIndexX + 1, fineIndexY, fineIndexZ + 1]
                            )

                    if not is_primal[0] and not is_primal[1] and is_primal[2]:
                        if qty == "Bz":
                            coarseData[
                                coarseLocalIndexX, coarseLocalIndexY, coarseLocalIndexZ
                            ] = 0.25 * (
                                fineData[fineIndexX, fineIndexY, fineIndexZ]
                                + fineData[fineIndexX + 1, fineIndexY, fineIndexZ]
                                + fineData[fineIndexX, fineIndexY + 1, fineIndexZ]
                                + fineData[fineIndexX + 1, fineIndexY + 1, fineIndexZ]
                            )

                    elif not is_primal[0] and is_primal[1] and is_primal[2]:  # Ex
                        value = 0.0
                        # dual in X means taking the two dual fine nodes
                        # around the coarse one we want the value for.
                        ishift = (0, 1)
                        iweights = (0.5, 0.5)

                        jshift = (-1, 0, 1)
                        jweights = (0.25, 0.5, 0.25)

                        kshift = (-1, 0, 1)
                        kweights = (0.25, 0.5, 0.25)
                        for i, iweight in zip(ishift, iweights):
                            yvalue = 0.0
                            for j, jweight in zip(jshift, jweights):
                                zvalue = 0.0
                                for k, kweight in zip(kshift, kweights):
                                    zvalue += (
                                        fineData[fineIndexX + i][fineIndexY + j][
                                            fineIndexZ + k
                                        ]
                                        * kweight
                                    )
                                yvalue += zvalue * jweight
                            value += yvalue * iweight
                        coarseData[
                            coarseLocalIndexX, coarseLocalIndexY, coarseLocalIndexZ
                        ] = value

                    elif (
                        is_primal[0] and not is_primal[1] and is_primal[2]
                    ):  # Ey is  PDP
                        value = 0.0
                        ishift = (-1, 0, 1)
                        iweights = (0.25, 0.5, 0.25)

                        # dual in Y means taking the two dual fine nodes
                        # around the coarse one we want the value for.
                        jshift = (0, 1)
                        jweights = (0.5, 0.5)

                        kshift = (-1, 0, 1)
                        kweights = (0.25, 0.5, 0.25)
                        for i, iweight in zip(ishift, iweights):
                            yvalue = 0.0
                            for j, jweight in zip(jshift, jweights):
                                zvalue = 0.0
                                for k, kweight in zip(kshift, kweights):
                                    zvalue += (
                                        fineData[fineIndexX + i][fineIndexY + j][
                                            fineIndexZ + k
                                        ]
                                        * kweight
                                    )
                                yvalue += zvalue * jweight
                            value += yvalue * iweight
                        coarseData[
                            coarseLocalIndexX, coarseLocalIndexY, coarseLocalIndexZ
                        ] = value

                    elif (
                        is_primal[0] and is_primal[1] and not is_primal[2]
                    ):  # Ez is PPD
                        value = 0.0
                        ishift = (-1, 0, 1)
                        iweights = (0.25, 0.5, 0.25)

                        jshift = (-1, 0, 1)
                        jweights = (0.25, 0.5, 0.25)

                        # dual in Z means taking the two dual fine nodes
                        # around the coarse one we want the value for.
                        kshift = (0, 1)
                        kweights = (0.5, 0.5)
                        for i, iweight in zip(ishift, iweights):
                            yvalue = 0.0
                            for j, jweight in zip(jshift, jweights):
                                zvalue = 0.0
                                for k, kweight in zip(kshift, kweights):
                                    zvalue += (
                                        fineData[fineIndexX + i][fineIndexY + j][
                                            fineIndexZ + k
                                        ]
                                        * kweight
                                    )
                                yvalue += zvalue * jweight
                            value += yvalue * iweight
                        coarseData[
                            coarseLocalIndexX, coarseLocalIndexY, coarseLocalIndexZ
                        ] = value
