#!/usr/bin/env python3
# This scripts aims at generate some signal on two levels (cos, sin, other functions)
# and then to perform a coarsening on a given region

import os
import sys
import numpy as np

from pyphare.core.box import Box
import pyphare.core.box as boxm
from pyphare.core import gridlayout
from pyphare.core.gridlayout import directions
from pyphare.core.phare_utilities import refinement_ratio
from pyphare.pharesee.hierarchy import FieldData

def exec_fn(xyz, fn):
    ndim = len(xyz)

    if ndim == 2:
        xx, yy = np.meshgrid(*xyz, indexing="ij")
        return fn(0, xx) * fn(1, yy)

    if ndim == 3:
        xx, yy, zz = np.meshgrid(*xyz, indexing="ij")
        return fn(0, xx) * fn(1, yy) * fn(2, zz)

    return fn(0, xyz[0])


def sin(xyz, L):
    return exec_fn(xyz, lambda i, v: np.sin((2 * np.pi / L[i]) * v))


def cos(xyz, L):
    return exec_fn(xyz, lambda i, v: np.cos((2 * np.pi / L[i]) * v))


functionList = [cos, sin]


def dump(ndim, path, quantity):

    origin = [0] * ndim

    coarseLayout = gridlayout.GridLayout(
        Box([0] * ndim, [39] * ndim), origin=origin, dl=[0.2] * ndim
    )

    fineLayout = gridlayout.GridLayout(
        Box([18] * ndim, [37] * ndim), origin=origin, dl=[0.1] * ndim
    )

    domainSize = coarseLayout.box.shape * coarseLayout.dl

    coarseCoords, fineCoords, is_primal = [], [], []

    for direction in directions[:ndim]:
        coarseCoords += [coarseLayout.yeeCoordsFor(quantity, direction)]
        fineCoords += [fineLayout.yeeCoordsFor(quantity, direction)]
        is_primal += [gridlayout.yee_element_is_primal(quantity, direction)]

    lower = [10] * ndim
    upper = [15 + primal for primal in is_primal]
    coarseBox = Box(lower, np.asarray(upper) - 1)

    for fn_idx, function in enumerate(functionList):
        fine = function(fineCoords, domainSize)
        coarse = function(coarseCoords, domainSize)
        afterCoarse = np.copy(coarse)

        coarseField = FieldData(coarseLayout, quantity, data=coarse)
        fineField   = FieldData(fineLayout, quantity, data=fine)
        coarsen(quantity, coarseField, fineField, coarseBox, fine, afterCoarse)

        if fn_idx == 0:
            fineData = np.copy(fine)
            coarseData = np.copy(coarse)
            afterCoarseData = np.copy(afterCoarse)
        else:
            fineData = np.vstack((fineData, fine))
            coarseData = np.vstack((coarseData, coarse))
            afterCoarseData = np.vstack((afterCoarseData, afterCoarse))

    file_datas = {
        f"{quantity}_fine_original{ndim}d.txt": fineData,
        f"{quantity}_coarse_original{ndim}d.txt": coarseData,
        f"{quantity}_coarse_linear_coarsed_{ndim}d.txt": afterCoarseData,
    }

    for filename, data in file_datas.items():
        np.savetxt(os.path.join(path, filename), data, delimiter=" ")



def main(path="./"):

    if len(sys.argv) > 1:
        path = sys.argv[1]

    for ndim in [1, 2]:
        for EM in ["E", "B"]:
            for qty in ["x", "y", "z"]:
                dump(ndim, path, EM + qty)



def coarsen(qty, coarseField, fineField, coarseBox, fineData, coarseData):
    coarseLayout = coarseField.layout
    fineLayout = fineField.layout
    ndim = coarseLayout.box.ndim

    nGhosts = coarseLayout.nbrGhostFor(qty)
    coarseStartIndex = coarseLayout.physicalStartIndices(qty)
    fineOffset = fineLayout.box.lower - boxm.refine(coarseLayout.box, 2).lower

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
                coarseData[coarseLocalIndex] = (
                    fineData[fineIndex - 1] * .25 + fineData[fineIndex] * .5 + fineData[fineIndex + 1] * .25
                )
            else:
                coarseData[coarseLocalIndex] = (
                    fineData[fineIndex] * .5 + fineData[fineIndex + 1] * .5
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
                    left += fineData[fineIndexX - 1][fineIndexY - 1]  * .25
                    left += fineData[fineIndexX - 1][fineIndexY]      * .5
                    left += fineData[fineIndexX - 1][fineIndexY + 1]  * .25
                    middle += fineData[fineIndexX][fineIndexY - 1]    * .25
                    middle += fineData[fineIndexX][fineIndexY]        * .5
                    middle += fineData[fineIndexX][fineIndexY + 1]    * .25
                    right += fineData[fineIndexX + 1][fineIndexY - 1] * .25
                    right += fineData[fineIndexX + 1][fineIndexY]     * .5
                    right += fineData[fineIndexX + 1][fineIndexY + 1] * .25
                    coarseData[coarseLocalIndexX][coarseLocalIndexY] = (
                        left * .25 + middle * .5 + right * .25
                    )

                if is_primal[0] and not is_primal[1]:
                    left += fineData[fineIndexX - 1][fineIndexY]      * .5
                    left += fineData[fineIndexX - 1][fineIndexY + 1]  * .5
                    middle += fineData[fineIndexX][fineIndexY]        * .5
                    middle += fineData[fineIndexX][fineIndexY + 1]    * .5
                    right += fineData[fineIndexX + 1][fineIndexY]     * .5
                    right += fineData[fineIndexX + 1][fineIndexY + 1] * .5
                    coarseData[coarseLocalIndexX][coarseLocalIndexY] = (
                        left * .25 + middle * .5 + right * .25
                    )

                if not is_primal[0] and is_primal[1]:
                    left += fineData[fineIndexX][fineIndexY - 1]      * .25
                    left += fineData[fineIndexX][fineIndexY]          * .5
                    left += fineData[fineIndexX][fineIndexY + 1]      * .25
                    right += fineData[fineIndexX + 1][fineIndexY - 1] * .25
                    right += fineData[fineIndexX + 1][fineIndexY]     * .5
                    right += fineData[fineIndexX + 1][fineIndexY + 1] * .25
                    coarseData[coarseLocalIndexX][coarseLocalIndexY] = (left * .5 + right * .5)

                if not any(is_primal):
                    left += fineData[fineIndexX][fineIndexY]          * .5
                    left += fineData[fineIndexX][fineIndexY + 1]      * .5
                    right += fineData[fineIndexX + 1][fineIndexY]     * .5
                    right += fineData[fineIndexX + 1][fineIndexY + 1] * .5
                    coarseData[coarseLocalIndexX][coarseLocalIndexY] = (left * .5 + right * .5)



if __name__ == "__main__":
    main()
