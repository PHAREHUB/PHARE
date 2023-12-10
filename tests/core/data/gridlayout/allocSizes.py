#!/usr/bin/env python
#!coding: utf-8
"""
this scripts defins 'allocSize.txt', a file containing expected values
for grid sizes, as a function of the interpolation order, the dimensionality (1,2,3)
and the quantity.
"""


import os
import sys

print(sys.path)
import gridparams
import utilities
from pyphare.core import gridlayout


class AllocParams(gridparams.GridParams):
    def __init__(self, dim, interpOrder):
        gridparams.GridParams.__init__(self, dim, interpOrder)
        self.alloc = ()
        self.allocDer = ()

    def setAlloc(self, quantity, gl):
        centeringX = gl.qtyCentering(quantity, "X")
        centeringY = gl.qtyCentering(quantity, "Y")
        centeringZ = gl.qtyCentering(quantity, "Z")

        if self.dim == 1:
            allocX = gl.allocSize(self.interpOrder, centeringX, self.nbrCell)
            allocDerX = gl.allocSizeDerived(self.interpOrder, centeringX, self.nbrCell)
        else:
            allocX = gl.allocSize(self.interpOrder, centeringX, self.nbrCell[0])
            allocDerX = gl.allocSizeDerived(
                self.interpOrder, centeringX, self.nbrCell[0]
            )

        if self.dim > 1:
            allocY = gl.allocSize(self.interpOrder, centeringY, self.nbrCell[1])
            allocDerY = gl.allocSizeDerived(
                self.interpOrder, centeringY, self.nbrCell[1]
            )
        if self.dim > 2:
            allocZ = gl.allocSize(self.interpOrder, centeringZ, self.nbrCell[2])
            allocDerZ = gl.allocSizeDerived(
                self.interpOrder, centeringZ, self.nbrCell[2]
            )

        if self.dim == 1:
            self.alloc = allocX
            self.allocDer = allocDerX

        if self.dim == 2:
            self.alloc = (allocX, allocY)
            self.allocDer = (allocDerX, allocDerY)

        if self.dim == 3:
            self.alloc = (allocX, allocY, allocZ)
            self.allocDer = (allocDerX, allocDerY, allocDerZ)


# ---------------------- MAIN CODE -----------------------------------------
def main(path="./"):
    if len(sys.argv) == 2:
        path = sys.argv[1]

    interpOrders = [1, 2, 3, 4]

    gl = gridlayout.GridLayout()

    directions = gl.directions
    quantities = [
        "Bx",
        "By",
        "Bz",
        "Ex",
        "Ey",
        "Ez",
        "Jx",
        "Jy",
        "Jz",
        "rho",
        "Vx",
        "Vy",
        "Vz",
        "P",
    ]

    nbrCellXList = [40, 40, 40]
    nbrCellYList = [0, 12, 12]
    nbrCellZList = [0, 0, 12]

    nbDimsList = [1, 2, 3]

    dxList = [0.1, 0.1, 0.1]  # 1D, 2D, 3D cases
    dyList = [0.0, 0.1, 0.1]
    dzList = [0.0, 0.0, 0.1]

    maxNbrDim = 3

    baseName = "allocSizes"

    # out_1D = open(os.path.join(path, baseName + '_1d.txt'), 'w')
    # out_2D = open(os.path.join(path, baseName + '_2d.txt'), 'w')
    # out_3D = open(os.path.join(path, baseName + '_3d.txt'), 'w')

    outFilenameBase = os.path.join(path, baseName)
    outFiles = []

    for interpOrder in interpOrders:
        filenames = [
            outFilenameBase + "_" + str(dim) + "d_O" + str(interpOrder) + ".txt"
            for dim in nbDimsList
        ]
        outFiles.append([open(f, "w") for f in filenames])

    for interpOrder, outFilesDim in zip(interpOrders, outFiles):
        for dim, outFile, nbrCellX, nbrCellY, nbrCellZ, dx, dy, dz in zip(
            nbDimsList,
            outFilesDim,
            nbrCellXList,
            nbrCellYList,
            nbrCellZList,
            dxList,
            dyList,
            dzList,
        ):
            params = AllocParams(dim, interpOrder)
            params.setNbrCell(nbrCellX, nbrCellY, nbrCellZ)
            params.setDl(dx, dy, dz)

            for quantity in quantities:
                params.setAlloc(quantity, gl)

                outString = "{} {} {} {} {}\n".format(
                    quantities.index(quantity),
                    params.nbrCell,
                    params.dl,
                    params.alloc,
                    params.allocDer,
                )

                outString = utilities.removeTupleFormat(outString)

                outFile.write(outString)

    for files in outFiles:
        for f in files:
            f.close()


if __name__ == "__main__":
    main()
