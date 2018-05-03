#!/usr/bin/env python
#!coding: utf-8
"""
this script aims at producing expected values for cell center coordinates
it writes the following files:
    - centeredCoords_summary.txt : contains physical Start and End indexes
        as a function of the interpolation order, dimensionality

    - and centeredCoords_ord*.txt which contain the coordinates of all cell
    centers as a function of the interpolation order and dimensionality
"""

import numpy as np

import sys

import gridlayout
import utilities
import os
import gridparams

class CenteredCoordParams(gridparams.GridParams):
    def __init__(self, dim, interpOrder):
        gridparams.GridParams.__init__(self, dim, interpOrder)

        self.origin = ()

        self.iStart = ()
        self.iEnd = ()


    def setCoord(self, gl, originPosition, centering):

        if (self.dim == 1):
            self.origin = originPosition[0]
            self.iStart = gl.physicalStartIndex(self.interpOrder, centering)
            self.iEnd = gl.physicalEndIndex(self.interpOrder, centering, self.nbrCell) -\
                        gl.isDual(centering)


        if (self.dim > 1):

            iStartX = gl.physicalStartIndex(self.interpOrder, centering[0])
            iEndX = gl.physicalEndIndex(self.interpOrder, centering[0], self.nbrCell[0]) -\
                    gl.isDual(centering[0])

            iStartY = gl.physicalStartIndex(self.interpOrder, centering[1])
            iEndY = gl.physicalEndIndex(self.interpOrder, centering[1], self.nbrCell[1]) -\
                    gl.isDual(centering[1])


        if (self.dim > 2):

            iStartZ = gl.physicalStartIndex(self.interpOrder, centering[2])

            iEndZ = gl.physicalEndIndex(self.interpOrder, centering[2], self.nbrCell[2]) -\
                    gl.isDual(centering[2])

        if (self.dim == 2):
            self.origin = (originPosition[0], originPosition[1])
            self.iStart = (iStartX, iStartY)
            self.iEnd = (iEndX, iEndY)

        elif (self.dim == 3):
            self.origin = (originPosition[0], originPosition[1], originPosition[2])
            self.iStart = (iStartX, iStartY, iStartZ)
            self.iEnd = (iEndX, iEndY, iEndZ)


def getCellCentered(dimension):
    if (dimension == 1):
        return 'primal'
    elif (dimension == 2):
        return ('primal','primal')
    elif (dimension == 3):
        return ('primal','primal','primal')




# ---- Get coordinate methods -------------------------
# ds stands for dx or dy or dz
#
# This method returns a point
#
def centeredCoords(primalIndex, startIndex, dl, origin):

    # a cell-centered coordinate is always dual
    halfCell = 0.5

    x = ( (primalIndex - startIndex) + halfCell )*dl + origin

    return x



# ---------------------- MAIN CODE -----------------------------------------




def main(path='./'):

    if len(sys.argv) == 2:
        path = sys.argv[1]

    interpOrders = [1,2,3,4]

    nbDimsList= [1, 2, 3]


    nbrCellXList=[40, 40, 40]
    nbrCellYList=[ 0, 12, 12]
    nbrCellZList=[ 0,  0, 12]

    dxList=[0.1, 0.1, 0.1]
    dyList=[0. , 0.1, 0.1]
    dzList=[0. , 0. , 0.1]


    originPosition = [0., 0., 0.]

    gl = gridlayout.GridLayout()

    # ------- Debug commands -------
    # for icase in icase_l:
    #     idim = dim_l[icase]
    #     order = interpOrder_l[icase]
    #     print( "Interpolation order = %d" % order )
    #     print( "Nbr of cells = %d" %  nbrCells[Direction_l[idim][1]][icase])
    #     print( "Nbr of ghost cells on the primal mesh = %d on each side" %
    #             gl.nbrGhostsPrimal(order) )

    # ------------------------------

    outSummary1D = open(os.path.join(path,"centeredCoords_summary_1d.txt"), "w")
    outSummary2D = open(os.path.join(path,"centeredCoords_summary_2d.txt"), "w")
    outSummary3D = open(os.path.join(path,"centeredCoords_summary_3d.txt"), "w")

    outSummaries =[outSummary1D, outSummary2D, outSummary3D]

    outActualValues1D = open(os.path.join(path,"centeredCoords_values_1d.txt"), "w")
    outActualValues2D = open(os.path.join(path,"centeredCoords_values_2d.txt"), "w")
    outActualValues3D = open(os.path.join(path,"centeredCoords_values_3d.txt"), "w")

    outActualValuesFiles = [outActualValues1D, outActualValues2D, outActualValues3D]

    for interpOrder in interpOrders:
        for dimension,outSummary, outValues,nbrCellX,nbrCellY,\
            nbrCellZ,dx,dy,dz in zip(nbDimsList, outSummaries,
                                     outActualValuesFiles,
                                     nbrCellXList,nbrCellYList,
                                     nbrCellZList,dxList,dyList,
                                     dzList):

            params = CenteredCoordParams(dimension,interpOrder);
            params.setNbrCell(nbrCellX,nbrCellY,nbrCellZ)
            params.setDl(dx,dy,dz)

            centering = getCellCentered(dimension)

            params.setCoord(gl, originPosition, centering)


            summaryBasePart = "{} {} {} ".format(
                               params.interpOrder,
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
                    outValuesString = "{} {} {}\n".format(
                            params.interpOrder,
                            position,
                            centeredCoords(position, params.iStart,params.dl, params.origin)
                        )

                    outValues.write(utilities.removeTupleFormat(outValuesString))

            elif dimension == 2:
                for positionX in np.arange(params.iStart[0], params.iEnd[0] + 1):
                    for positionY in np.arange(params.iStart[1], params.iEnd[1] + 1):
                        position = (positionX, positionY)
                        centered = (centeredCoords(positionX, params.iStart[0], params.dl[0],
                                                   params.origin[0]),
                                    centeredCoords(positionY, params.iStart[1], params.dl[1],
                                                   params.origin[1]))

                        outValuesString = "{} {} {}\n".format(
                            params.interpOrder,
                            position,
                            centered
                        )

                        outValues.write(utilities.removeTupleFormat(outValuesString))

            elif dimension == 3:
                for positionX in np.arange(params.iStart[0], params.iEnd[0] + 1):
                    for positionY in np.arange(params.iStart[1], params.iEnd[1] + 1):
                        for positionZ in np.arange(params.iStart[2], params.iEnd[2] + 1):

                            position = (positionX, positionY, positionZ)
                            centered = (centeredCoords(positionX, params.iStart[0], params.dl[0],
                                                       params.origin[0]),
                                        centeredCoords(positionY, params.iStart[1], params.dl[1],
                                                       params.origin[1]),
                                        centeredCoords(positionZ, params.iStart[2], params.dl[2],
                                                       params.origin[2]))

                            outValuesString = "{} {} {}\n".format(
                                params.interpOrder,
                                position,
                                centered
                            )

                            outValues.write(utilities.removeTupleFormat(outValuesString))




    for outSummary,outValues in zip(outSummaries, outActualValuesFiles):
        outSummary.close()
        outValues.close()


if __name__ == "__main__":
    main()

