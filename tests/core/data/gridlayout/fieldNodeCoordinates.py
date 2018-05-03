#!/usr/bin/env python
#!coding: utf-8

import numpy as np
#import math

import sys

import gridlayout
import os
import gridparams
import cellCenteredCoordinates
import utilities

# TODO : FieldNode coords is general case of cellCenteredCoord
#        Which means this has to be fully refactor
#        The only difference is names, loop on hq, and which function to call
#        on the gridlayout side


# ---- Get coordinate methods -------------------------
# ds stands for dx or dy or dz
#
# This method returns a point
#
def fieldCoords(primalIndex, startIndex, quantity, direction, dl, origin):
    halfCell = 0.
    gl  = gridlayout.GridLayout()

    if gl.qtyCentering(quantity, direction) == 'dual':
        halfCell = 0.5

    x = ( (primalIndex - startIndex) + halfCell )*dl + origin

    return x

class FieldNodeCoordParams(cellCenteredCoordinates.CenteredCoordParams):
    def __init__(self, dim, interpOrder):
        cellCenteredCoordinates.CenteredCoordParams.__init__(self, dim, interpOrder)

def getQtyCentering(gl,quantity, dimension):
    if (dimension == 1):
        return gl.qtyCentering(quantity, 'X')
    elif (dimension == 2):
        return (gl.qtyCentering(quantity, 'X'),
                gl.qtyCentering(quantity, 'Y'))
    elif (dimension == 3):
        return (gl.qtyCentering(quantity, 'X'),
                gl.qtyCentering(quantity, 'Y'),
                gl.qtyCentering(quantity, 'Z'))

# ---------------------- MAIN CODE -----------------------------------------

def main(path='./'):

    if len(sys.argv) == 2:
        path = sys.argv[1]

    # TODO: refactor this for all python since they are the same

    interpOrders = [1,2,3,4]

    nbDimsList= [1, 2, 3]


    nbrCellXList=[40, 40, 40]
    nbrCellYList=[ 0, 12, 12]
    nbrCellZList=[ 0,  0, 12]

    dxList=[0.01, 0.01, 0.01]
    dyList=[0. , 0.01, 0.01]
    dzList=[0. , 0. , 0.01]


    originPosition = [0., 0., 0.]

    gl = gridlayout.GridLayout()

    # TODO : end Todo

    # TODO: FieldCoords and CenteredCoords share a common base, refactor this


    outSummary1D = open(os.path.join(path,"fieldCoords_summary_1d.txt"), "w")
    outSummary2D = open(os.path.join(path,"fieldCoords_summary_2d.txt"), "w")
    outSummary3D = open(os.path.join(path,"fieldCoords_summary_3d.txt"), "w")

    outSummaries =[outSummary1D, outSummary2D, outSummary3D]

    outActualValues1D = open(os.path.join(path,"fieldCoords_values_1d.txt"), "w")
    outActualValues2D = open(os.path.join(path,"fieldCoords_values_2d.txt"), "w")
    outActualValues3D = open(os.path.join(path,"fieldCoords_values_3d.txt"), "w")

    outActualValuesFiles = [outActualValues1D, outActualValues2D, outActualValues3D]



    for interpOrder in interpOrders:
        for dimension,outSummary, outValues,nbrCellX,nbrCellY,\
            nbrCellZ,dx,dy,dz in zip(nbDimsList, outSummaries,
                                     outActualValuesFiles,
                                     nbrCellXList,nbrCellYList,
                                     nbrCellZList,dxList,dyList,
                                     dzList):

            for quantity in gl.hybridQuantities:

                params = FieldNodeCoordParams(dimension,interpOrder);
                params.setNbrCell(nbrCellX,nbrCellY,nbrCellZ)
                params.setDl(dx,dy,dz)

                centering = getQtyCentering(gl,quantity,dimension)

                params.setCoord(gl, originPosition, centering)


                summaryBasePart = "{} {} {} {} ".format(
                    params.interpOrder,
                    quantity,
                    params.nbrCell,
                    params.dl)


                summaryGridLayoutPart = "{} {} {}\n".format(
                    params.iStart,
                    params.iEnd,
                    params.origin)

                outSummaryString = summaryBasePart + summaryGridLayoutPart

                outSummaryString = utilities.removeTupleFormat(outSummaryString)


                outSummary.write(outSummaryString)



                if dimension == 1:
                    for position in np.arange(params.iStart, params.iEnd + 1):
                        outValuesString = "{} {} {} {}\n".format(
                            params.interpOrder,
                            quantity,
                            position,
                            fieldCoords(position, params.iStart,quantity,'X',params.dl, params.origin)
                        )

                        outValues.write(utilities.removeTupleFormat(outValuesString))

                elif dimension == 2:
                    for positionX in np.arange(params.iStart[0], params.iEnd[0] + 1):
                        for positionY in np.arange(params.iStart[1], params.iEnd[1] + 1):
                            position = (positionX, positionY)
                            centered = (fieldCoords(positionX, params.iStart[0],quantity,'X', params.dl[0],
                                                       params.origin[0]),
                                        fieldCoords(positionY, params.iStart[1],quantity,'Y', params.dl[1],
                                                       params.origin[1]))

                            outValuesString = "{} {} {} {}\n".format(
                                params.interpOrder,
                                quantity,
                                position,
                                centered
                            )

                            outValues.write(utilities.removeTupleFormat(outValuesString))

                elif dimension == 3:
                    for positionX in np.arange(params.iStart[0], params.iEnd[0] + 1):
                        for positionY in np.arange(params.iStart[1], params.iEnd[1] + 1):
                            for positionZ in np.arange(params.iStart[2], params.iEnd[2] + 1):

                                position = (positionX, positionY, positionZ)
                                centered = (fieldCoords(positionX, params.iStart[0],quantity,'X', params.dl[0],
                                                           params.origin[0]),
                                            fieldCoords(positionY, params.iStart[1],quantity,'Y', params.dl[1],
                                                           params.origin[1]),
                                            fieldCoords(positionZ, params.iStart[2],quantity,'Z', params.dl[2],
                                                           params.origin[2]))

                                outValuesString = "{} {} {} {}\n".format(
                                    params.interpOrder,
                                    quantity,
                                    position,
                                    centered
                                )

                                outValues.write(utilities.removeTupleFormat(outValuesString))




    for outSummary,outValues in zip(outSummaries, outActualValuesFiles):
        outSummary.close()
        outValues.close()



if __name__ == "__main__":
    main()
