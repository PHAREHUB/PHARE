#!/usr/bin/env python
#!coding: utf-8

import os
import sys
import math
import numpy as np

from pyphare.core import gridlayout

import utilities
import fieldNodeCoordinates
import cellCenteredCoordinates


class DerivParams(cellCenteredCoordinates.CenteredCoordParams):
    def __init__(self, dim, interpOrder):
        cellCenteredCoordinates.CenteredCoordParams.__init__(self, dim, interpOrder)

        self.iGhostStart = ()
        self.iGhostEnd = ()

    def setCoord(self, gl, origin, centering, derivCentering):
        cellCenteredCoordinates.CenteredCoordParams.setCoord(
            self, gl, origin, derivCentering
        )

        if self.dim == 1:
            self.iGhostStart = gl.ghostStartIndex()
            self.iGhostEnd = gl.ghostEndIndex(self.interpOrder, centering, self.nbrCell)


#     def setFunction(self, functionName):
#         self.functionName = functionName

# def minus_sin(coord):
#     return - math.sin(coord)

# def constant(coord):
#     return 1.0

# def unity(coord):
#     return coord

# derivedTable = { 'cos_sin_cosh' : (minus_sin, math.cos, math.sinh),
#                  'linear3': (constant, constant, constant)}

# functionTable = { 'cos_sin_cosh' : (math.cos, math.sin, math.cosh),
#                   'linear3' : (unity, unity, unity)}


def primitive(coord):
    return math.cos(coord)


# def callFunction(function, coord, dim):
#     if (dim == 1):
#         return function[0](coord)
#     elif (dim == 2):
#         return ( function[0](coord[0]), function[1](coord[1]) )
#     elif (dim == 3):
#         return ( function[0](coord[0]), function[1](coord[1]),
#                  function[2](coord[2]))


# def functions(functionName, coord, dimension):
#     return callFunction(functionTable[functionName], coord, dimension)

# def functionDerived(functionName, coord, dimension):
#     return callFunction(derivedTable[functionName], coord, dimension)


# ---------------------- MAIN CODE -----------------------------------------


def main(path="./"):
    import scipy.misc  # see doc/conventions.md #2.2.1

    if len(sys.argv) == 2:
        path = sys.argv[1]

    # TODO: refactor this for all python since they are the same

    interpOrders = [1, 2, 3]

    nbDimsList = [1]  # , 2, 3]

    nbrCellXList = [40, 40, 40]
    nbrCellYList = [0, 12, 12]
    nbrCellZList = [0, 0, 12]

    dxList = [0.01, 0.01, 0.01]
    dyList = [0.0, 0.01, 0.01]
    dzList = [0.0, 0.0, 0.01]

    originPosition = [0.0, 0.0, 0.0]

    # TODO : end Todo

    # TODO: FieldCoords and CenteredCoords share a common base, refactor this

    outSummary1D = open(os.path.join(path, "deriv_summary_1d.txt"), "w")
    outSummary2D = open(os.path.join(path, "deriv_summary_2d.txt"), "w")
    outSummary3D = open(os.path.join(path, "deriv_summary_3d.txt"), "w")

    outSummaries = [outSummary1D, outSummary2D, outSummary3D]

    outActualValues1D = open(os.path.join(path, "deriv_values_1d.txt"), "w")
    outActualValues2D = open(os.path.join(path, "deriv_values_2d.txt"), "w")
    outActualValues3D = open(os.path.join(path, "deriv_values_3d.txt"), "w")

    outActualValuesFiles = [outActualValues1D, outActualValues2D, outActualValues3D]

    outDerivedValues1D = open(os.path.join(path, "deriv_derived_values_1d.txt"), "w")

    for interpOrder in interpOrders:
        gl = gridlayout.GridLayout(interp_order=interpOrder)

        for (
            dimension,
            outSummary,
            outValues,
            nbrCellX,
            nbrCellY,
            nbrCellZ,
            dx,
            dy,
            dz,
        ) in zip(
            nbDimsList,
            outSummaries,
            outActualValuesFiles,
            nbrCellXList,
            nbrCellYList,
            nbrCellZList,
            dxList,
            dyList,
            dzList,
        ):
            # TODO : the association depend on the direction
            #        Here we derive in the direction X
            quantities = ("Ex", "Ey")
            derivedQuantities = ("rho", "Bz")

            for quantity, derivQuantity in zip(quantities, derivedQuantities):
                params = DerivParams(dimension, interpOrder)
                params.setNbrCell(nbrCellX, nbrCellY, nbrCellZ)
                params.setDl(dx, dy, dz)

                centering = fieldNodeCoordinates.getQtyCentering(
                    gl, quantity, dimension
                )
                derivedCentering = fieldNodeCoordinates.getQtyCentering(
                    gl, derivQuantity, dimension
                )

                params.setCoord(gl, originPosition, centering, derivedCentering)

                summaryBasePart = "{} {} {} {} {} ".format(
                    params.interpOrder,
                    quantity,
                    derivQuantity,
                    params.nbrCell,
                    params.dl,
                )

                summaryGridLayoutPart = " {} {} {} {} {}\n".format(
                    params.iGhostStart,
                    params.iGhostEnd,
                    params.iStart,
                    params.iEnd,
                    params.origin,
                )

                outSummaryString = summaryBasePart + summaryGridLayoutPart

                outSummaryString = utilities.removeTupleFormat(outSummaryString)

                outSummary.write(outSummaryString)

                if dimension == 1:
                    for position in np.arange(params.iGhostStart, params.iGhostEnd):
                        coord = fieldNodeCoordinates.fieldCoords(
                            position, quantity, "X", params.dl, gl
                        )
                        functionValue = primitive(coord)

                        outValuesString = "{} {} {} {}\n".format(
                            params.interpOrder,
                            quantity,
                            position,
                            functionValue,
                        )

                        outValues.write(utilities.removeTupleFormat(outValuesString))

                    for position in np.arange(params.iStart, params.iEnd + 1):
                        coord = fieldNodeCoordinates.fieldCoords(
                            position, derivQuantity, "X", params.dl, gl
                        )

                        derivedValue = scipy.misc.derivative(
                            primitive, coord, dx, order=3
                        )

                        outDerivedValuesString = "{} {} {} {}\n".format(
                            params.interpOrder, derivQuantity, position, derivedValue
                        )

                        outDerivedValues1D.write(
                            utilities.removeTupleFormat(outDerivedValuesString)
                        )

    for outSummary, outValues in zip(outSummaries, outActualValuesFiles):
        outSummary.close()
        outValues.close()

    outDerivedValues1D.close()


if __name__ == "__main__":
    main()
