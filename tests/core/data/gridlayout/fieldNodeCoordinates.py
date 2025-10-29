#!/usr/bin/env python
#!coding: utf-8


import os
import sys
import utilities
import numpy as np
import cellCenteredCoordinates

from pyphare.core import gridlayout

# TODO : FieldNode coords is general case of cellCenteredCoord
#        Which means this has to be fully refactor
#        The only difference is names, loop on hq, and which function to call
#        on the gridlayout side


# ---- Get coordinate methods -------------------------
# ds stands for dx or dy or dz
#
# This method returns a point
#
def fieldCoords(primalIndex, quantity, direction, dl, gl):
    halfCell = 0.0

    if gl.qtyCentering(quantity, direction) == "dual":
        halfCell = 0.5

    return (primalIndex + halfCell) * dl


class FieldNodeCoordParams(cellCenteredCoordinates.CenteredCoordParams):
    def __init__(self, dim, interpOrder):
        cellCenteredCoordinates.CenteredCoordParams.__init__(self, dim, interpOrder)


def getQtyCentering(gl, quantity, dimension):
    if dimension == 1:
        return gl.qtyCentering(quantity, "X")
    elif dimension == 2:
        return (gl.qtyCentering(quantity, "X"), gl.qtyCentering(quantity, "Y"))
    elif dimension == 3:
        return (
            gl.qtyCentering(quantity, "X"),
            gl.qtyCentering(quantity, "Y"),
            gl.qtyCentering(quantity, "Z"),
        )


# ---------------------- MAIN CODE -----------------------------------------


def main(path="./"):
    if len(sys.argv) == 2:
        path = sys.argv[1]

    # TODO: refactor this for all python since they are the same

    interpOrders = [1, 2, 3, 4]

    nbDimsList = [1, 2, 3]

    nbrCellXList = [40, 40, 40]
    nbrCellYList = [0, 12, 12]
    nbrCellZList = [0, 0, 12]

    dxList = [0.01, 0.01, 0.01]
    dyList = [0.0, 0.01, 0.01]
    dzList = [0.0, 0.0, 0.01]

    originPosition = [0.0, 0.0, 0.0]

    # TODO : end Todo

    # TODO: FieldCoords and CenteredCoords share a common base, refactor this

    baseNameSummary = "fieldCoords_summary"
    baseNameValues = "fieldCoords_values"

    outFilenameBaseSummary = os.path.join(path, baseNameSummary)
    outFilenameBaseValues = os.path.join(path, baseNameValues)
    outSummaries = []
    outValues = []

    for interpOrder in interpOrders:
        filenamesSum = [
            outFilenameBaseSummary + "_" + str(dim) + "d_O" + str(interpOrder) + ".txt"
            for dim in nbDimsList
        ]
        filenamesVal = [
            outFilenameBaseValues + "_" + str(dim) + "d_O" + str(interpOrder) + ".txt"
            for dim in nbDimsList
        ]
        outSummaries.append([open(f, "w") for f in filenamesSum])
        outValues.append([open(f, "w") for f in filenamesVal])

    for interpOrder, outFilesSumDim, outFilesValDim in zip(
        interpOrders, outSummaries, outValues
    ):
        gl = gridlayout.GridLayout(interp_order=interpOrder)

        for (
            dimension,
            outFileS,
            outFileV,
            nbrCellX,
            nbrCellY,
            nbrCellZ,
            dx,
            dy,
            dz,
        ) in zip(
            nbDimsList,
            outFilesSumDim,
            outFilesValDim,
            nbrCellXList,
            nbrCellYList,
            nbrCellZList,
            dxList,
            dyList,
            dzList,
        ):
            for quantity in gl.hybridQuantities:
                params = FieldNodeCoordParams(dimension, interpOrder)
                params.setNbrCell(nbrCellX, nbrCellY, nbrCellZ)
                params.setDl(dx, dy, dz)

                centering = getQtyCentering(gl, quantity, dimension)

                params.setCoord(gl, originPosition, centering)

                summaryBasePart = "{} {} {} ".format(
                    quantity, params.nbrCell, params.dl
                )

                summaryGridLayoutPart = "{} {} {}\n".format(
                    params.iStart, params.iEnd, params.origin
                )

                outSummaryString = summaryBasePart + summaryGridLayoutPart

                outSummaryString = utilities.removeTupleFormat(outSummaryString)

                outFileS.write(outSummaryString)

                if dimension == 1:
                    for position in np.arange(params.iStart, params.iEnd + 1):
                        outValuesString = "{} {} {}\n".format(
                            quantity,
                            position,
                            fieldCoords(position, quantity, "X", params.dl, gl),
                        )

                        outFileV.write(utilities.removeTupleFormat(outValuesString))

                elif dimension == 2:
                    for positionX in np.arange(params.iStart[0], params.iEnd[0] + 1):
                        for positionY in np.arange(
                            params.iStart[1], params.iEnd[1] + 1
                        ):
                            position = (positionX, positionY)
                            centered = (
                                fieldCoords(positionX, quantity, "X", params.dl[0], gl),
                                fieldCoords(positionY, quantity, "Y", params.dl[1], gl),
                            )

                            outValuesString = "{} {} {}\n".format(
                                quantity, position, centered
                            )

                            outFileV.write(utilities.removeTupleFormat(outValuesString))

                elif dimension == 3:
                    for positionX in np.arange(params.iStart[0], params.iEnd[0] + 1):
                        for positionY in np.arange(
                            params.iStart[1], params.iEnd[1] + 1
                        ):
                            for positionZ in np.arange(
                                params.iStart[2], params.iEnd[2] + 1
                            ):
                                position = (positionX, positionY, positionZ)
                                centered = (
                                    fieldCoords(
                                        positionX, quantity, "X", params.dl[0], gl
                                    ),
                                    fieldCoords(
                                        positionY, quantity, "Y", params.dl[1], gl
                                    ),
                                    fieldCoords(
                                        positionZ, quantity, "Z", params.dl[2], gl
                                    ),
                                )

                                outValuesString = "{} {} {}\n".format(
                                    quantity, position, centered
                                )

                                outFileV.write(
                                    utilities.removeTupleFormat(outValuesString)
                                )

    for outFilesSumDim, outFilesValDim in zip(outSummaries, outValues):
        for f1, f2 in zip(outFilesSumDim, outFilesValDim):
            f1.close()
            f2.close()


if __name__ == "__main__":
    main()
