#!/usr/bin/env python
#!coding: utf-8
"""
This script prepares 'gridIndexing.txt'. The file contains expected values
for physical/ghostStart/EndIndexes in directions X, Y and Z for all hybrid
quantities and for all interpolation orders (1,2,3,4), in 1D, 2D and 3D.
"""

import numpy as np

import sys

import gridlayout
import os
import utilities
import gridparams

class IndexingParams(gridparams.GridParams):
    def __init__(self,dim,interpOrder):
        gridparams.GridParams.__init__(self,dim,interpOrder)

        self.physicalStart = ()
        self.physicalEnd = ()
        self.ghostStart = ()
        self.ghostEnd = ()

    def setIndexes(self,quantity,gl):

        centeringX = gl.qtyCentering(quantity, 'X')
        centeringY = gl.qtyCentering(quantity, 'Y')
        centeringZ = gl.qtyCentering(quantity, 'Z')

        if self.dim == 1:
            self.physicalStart = (gl.physicalStartIndex(self.interpOrder,centeringX))

            self.physicalEnd = (gl.physicalEndIndex(self.interpOrder,centeringX,self.nbrCell))

            self.ghostStart = (gl.ghostStartIndex())

            self.ghostEnd   = (gl.ghostEndIndex(self.interpOrder,centeringX,self.nbrCell))

        if self.dim == 2:
            self.physicalStart = (gl.physicalStartIndex(self.interpOrder,centeringX),
                             gl.physicalStartIndex(self.interpOrder,centeringY))

            self.physicalEnd = (gl.physicalEndIndex(self.interpOrder,centeringX,self.nbrCell[0]),
                           gl.physicalEndIndex(self.interpOrder,centeringY,self.nbrCell[1]))

            self.ghostStart = (gl.ghostStartIndex(),
                          gl.ghostStartIndex())

            self.ghostEnd   = (gl.ghostEndIndex(self.interpOrder,centeringX,self.nbrCell[0]),
                          gl.ghostEndIndex(self.interpOrder,centeringY,self.nbrCell[1]))

        if self.dim == 3:
            self.physicalStart = (gl.physicalStartIndex(self.interpOrder,centeringX),
                             gl.physicalStartIndex(self.interpOrder,centeringY),
                             gl.physicalStartIndex(self.interpOrder,centeringZ))

            self.physicalEnd = (gl.physicalEndIndex(self.interpOrder,centeringX,self.nbrCell[0]),
                           gl.physicalEndIndex(self.interpOrder,centeringY,self.nbrCell[1]),
                           gl.physicalEndIndex(self.interpOrder,centeringZ,self.nbrCell[2]))

            self.ghostStart = (gl.ghostStartIndex(),
                          gl.ghostStartIndex(),
                          gl.ghostStartIndex())

            self.ghostEnd   = (gl.ghostEndIndex(self.interpOrder,centeringX,self.nbrCell[0]),
                          gl.ghostEndIndex(self.interpOrder,centeringY,self.nbrCell[1]),
                          gl.ghostEndIndex(self.interpOrder,centeringZ,self.nbrCell[2]))

# ---------------------- MAIN CODE -----------------------------------------
def main(path='./'):

    if len(sys.argv) == 2:
        path = sys.argv[1]


    interpOrders=[1, 2, 3, 4]


    gl = gridlayout.GridLayout()

    directions = gl.directions
    quantities = gl.hybridQuantities

    nbrCellXList=[40, 40, 40]
    nbrCellYList=[ 0, 12, 12]
    nbrCellZList=[ 0,  0, 12]

    nbDimsList=[1, 2, 3]

    dxList=[0.1, 0.1, 0.1] # 1D, 2D and 3D cases
    dyList=[0. , 0.1, 0.1]
    dzList=[0. , 0. , 0.1]


    maxNbrDim = 3


    baseName = 'gridIndexing'

    out_1D = open(os.path.join(path, baseName + '_1d.txt'), 'w')
    out_2D = open(os.path.join(path, baseName + '_2d.txt'), 'w')
    out_3D = open(os.path.join(path, baseName + '_3d.txt'), 'w')

    outFiles = [out_1D, out_2D, out_3D]


    for interpOrder in interpOrders:
        for dim,outFile,nbrCellX,nbrCellY,nbrCellZ, dx,dy,dz in zip(nbDimsList, outFiles,nbrCellXList,nbrCellYList,
                                                                    nbrCellZList,dxList,dyList,dzList):

            params = IndexingParams(dim,interpOrder)

            params.setNbrCell(nbrCellX,nbrCellY,nbrCellZ)
            params.setDl(dx,dy,dz)

            for quantity in quantities:

                params.setIndexes(quantity,gl)


                outString =  "{} {} {} {} {} {} {} {}\n".format( interpOrder,
                               quantities.index(quantity),
                               params.nbrCell,
                               params.dl,
                               params.physicalStart,
                               params.physicalEnd,
                               params.ghostStart,
                               params.ghostEnd)

                outFile.write(utilities.removeTupleFormat(outString))


    for file  in outFiles:
        file.close()



if __name__ == "__main__":
    main()
