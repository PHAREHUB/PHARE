#!/usr/bin/env python
#!coding: utf-8



"""
this script produces expected values for the first order derivatives
implemented in GridLayout with a Yee implementation

We derive a primal quantity and a dual quantity.
"""



import numpy as np
import gridlayout
import os

class TestVariables (object) :

    def __init__(self) :

        """ This class only has a constructor
            in order to set the variables to test the derivatives
            all quantities are 3D

        """

        self.nbrCells = (50, 30, 40)
        self.meshSize = (0.1, 0.2, 0.3)
        self.interpOrders = (1, 2, 3)
        self.ByCentering = ('dual', 'primal', 'dual')
        self.EzCentering = ('primal', 'primal', 'dual')
        self.domainSize = tuple(n * m for n, m in zip(self.nbrCells, self.meshSize))


def test_deriv1D(path):
    # derivative along X:
    # By is dual in X, and derived along X it should lie on Ez
    # Ez is primal in X, and derived along X it should lie on By

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    for interpOrder in tv.interpOrders:

        filename_dxBy = 'dxBy_interpOrder_{}_1d.txt'.format(interpOrder)
        filename_dxEz = 'dxEz_interpOrder_{}_1d.txt'.format(interpOrder)

        By   = np.empty(layout.allocSize(interpOrder, tv.ByCentering[0], tv.nbrCells[0]))
        Ez   = np.empty(layout.allocSize(interpOrder, tv.EzCentering[0], tv.nbrCells[0]))

        dxBy = np.empty(layout.allocSizeDerived(interpOrder, tv.ByCentering[0], tv.nbrCells[0]), dtype=np.float64)
        dxEz = np.empty(layout.allocSizeDerived(interpOrder, tv.EzCentering[0], tv.nbrCells[0]), dtype=np.float64)

        nbrGhost_p = layout.nbrGhosts(interpOrder, 'primal')
        nbrGhost_d = layout.nbrGhosts(interpOrder, 'dual')

        x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
        x_primal = tv.meshSize[0]*np.arange(layout.allocSize(interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p

        By = np.cos(2*np.pi/tv.domainSize[0] * x_dual)
        Ez = np.cos(2*np.pi/tv.domainSize[0] * x_primal)

        psi_d_X = layout.physicalStartIndex(interpOrder, 'dual')
        pei_d_X = layout.physicalEndIndex(  interpOrder, 'dual',   tv.nbrCells[0])
        psi_p_X = layout.physicalStartIndex(interpOrder, 'primal')
        pei_p_X = layout.physicalEndIndex(  interpOrder, 'primal', tv.nbrCells[0])

        dxBy[psi_p_X:pei_p_X+1]  = (By[psi_d_X  :pei_d_X+2] - By[psi_d_X-1:pei_d_X+1])/tv.meshSize[0]
        dxEz[psi_d_X:pei_d_X+1]  = (Ez[psi_p_X+1:pei_p_X+1] - Ez[psi_p_X  :pei_p_X  ])/tv.meshSize[0]

        np.savetxt(os.path.join(path,filename_dxBy), dxBy, delimiter=" ")
        np.savetxt(os.path.join(path,filename_dxEz), dxEz, delimiter=" ")


def test_deriv2D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    for interpOrder in tv.interpOrders:

        filename_dxBy = 'dxBy_interpOrder_{}_2d.txt'.format(interpOrder)
        filename_dyBy = 'dyBy_interpOrder_{}_2d.txt'.format(interpOrder)
        filename_dxEz = 'dxEz_interpOrder_{}_2d.txt'.format(interpOrder)
        filename_dyEz = 'dyEz_interpOrder_{}_2d.txt'.format(interpOrder)

        By   = np.empty([layout.allocSize(interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.ByCentering[1], tv.nbrCells[1])], dtype=np.float64)
        Ez   = np.empty([layout.allocSize(interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)

        dxBy = np.empty([layout.allocSizeDerived(interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.ByCentering[1], tv.nbrCells[1])], dtype=np.float64)
        dyBy = np.empty([layout.allocSize(interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                         layout.allocSizeDerived(interpOrder, tv.ByCentering[1], tv.nbrCells[1])], dtype=np.float64)
        dxEz = np.empty([layout.allocSizeDerived(interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)
        dyEz = np.empty([layout.allocSize(interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSizeDerived(interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)

        nbrGhost_p = layout.nbrGhosts(interpOrder, 'primal')
        nbrGhost_d = layout.nbrGhosts(interpOrder, 'dual')

        x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
        y_dual   = tv.meshSize[1]*np.arange(layout.allocSize(interpOrder, 'dual'  , tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_d + tv.meshSize[1]*0.5
        x_primal = tv.meshSize[0]*np.arange(layout.allocSize(interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
        y_primal = tv.meshSize[1]*np.arange(layout.allocSize(interpOrder, 'primal', tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_p

        By = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_dual),
                          np.sin(2*np.pi/tv.domainSize[1] * y_primal), axes=0)

        Ez = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_primal),
                          np.sin(2*np.pi/tv.domainSize[1] * y_primal), axes=0)

        psi_d_X = layout.physicalStartIndex(interpOrder, 'dual')
        pei_d_X = layout.physicalEndIndex(  interpOrder, 'dual',   tv.nbrCells[0])
        psi_d_Y = layout.physicalStartIndex(interpOrder, 'dual')
        pei_d_Y = layout.physicalEndIndex(  interpOrder, 'dual',   tv.nbrCells[1])
        psi_p_X = layout.physicalStartIndex(interpOrder, 'primal')
        pei_p_X = layout.physicalEndIndex(  interpOrder, 'primal', tv.nbrCells[0])
        psi_p_Y = layout.physicalStartIndex(interpOrder, 'primal')
        pei_p_Y = layout.physicalEndIndex(  interpOrder, 'primal', tv.nbrCells[1])

        dxBy[psi_p_X:pei_p_X+1,:]  = (By[psi_d_X  :pei_d_X+2,:] - By[psi_d_X-1:pei_d_X+1,:])/tv.meshSize[0]
        dxEz[psi_d_X:pei_d_X+1,:]  = (Ez[psi_p_X+1:pei_p_X+1,:] - Ez[psi_p_X  :pei_p_X  ,:])/tv.meshSize[0]
        dyBy[:,psi_d_Y:pei_d_Y+1]  = (By[:,psi_p_Y+1:pei_p_Y+1] - By[:,psi_p_Y  :pei_p_Y  ])/tv.meshSize[1]
        dyEz[:,psi_d_Y:pei_d_Y+1]  = (Ez[:,psi_p_Y+1:pei_p_Y+1] - Ez[:,psi_p_Y  :pei_p_Y  ])/tv.meshSize[1]

        np.savetxt(os.path.join(path,filename_dxBy), dxBy.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dyBy), dyBy.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dxEz), dxEz.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dyEz), dyEz.flatten('C'), delimiter=" ")


def test_deriv3D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    for interpOrder in tv.interpOrders:

        filename_dxBy = 'dxBy_interpOrder_{}_3d.txt'.format(interpOrder)
        filename_dyBy = 'dyBy_interpOrder_{}_3d.txt'.format(interpOrder)
        filename_dzBy = 'dzBy_interpOrder_{}_3d.txt'.format(interpOrder)
        filename_dxEz = 'dxEz_interpOrder_{}_3d.txt'.format(interpOrder)
        filename_dyEz = 'dyEz_interpOrder_{}_3d.txt'.format(interpOrder)
        filename_dzEz = 'dzEz_interpOrder_{}_3d.txt'.format(interpOrder)

        By   = np.empty([layout.allocSize(interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.ByCentering[1], tv.nbrCells[1]),
                         layout.allocSize(interpOrder, tv.ByCentering[2], tv.nbrCells[2])], dtype=np.float64)
        Ez   = np.empty([layout.allocSize(interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                         layout.allocSize(interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)

        dxBy = np.empty([layout.allocSizeDerived(interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.ByCentering[1], tv.nbrCells[1]),
                         layout.allocSize(interpOrder, tv.ByCentering[2], tv.nbrCells[2])], dtype=np.float64)
        dyBy = np.empty([layout.allocSize(interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                         layout.allocSizeDerived(interpOrder, tv.ByCentering[1], tv.nbrCells[1]),
                         layout.allocSize(interpOrder, tv.ByCentering[2], tv.nbrCells[2])], dtype=np.float64)
        dzBy = np.empty([layout.allocSize(interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.ByCentering[1], tv.nbrCells[1]),
                         layout.allocSizeDerived(interpOrder, tv.ByCentering[2], tv.nbrCells[2])], dtype=np.float64)

        dxEz = np.empty([layout.allocSizeDerived(interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                         layout.allocSize(interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)
        dyEz = np.empty([layout.allocSize(interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSizeDerived(interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                         layout.allocSize(interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)
        dzEz = np.empty([layout.allocSize(interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSize(interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                         layout.allocSizeDerived(interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)

        nbrGhost_p = layout.nbrGhosts(interpOrder, 'primal')
        nbrGhost_d = layout.nbrGhosts(interpOrder, 'dual')

        x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
        y_dual   = tv.meshSize[1]*np.arange(layout.allocSize(interpOrder, 'dual'  , tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_d + tv.meshSize[1]*0.5
        z_dual   = tv.meshSize[2]*np.arange(layout.allocSize(interpOrder, 'dual'  , tv.nbrCells[2])) - tv.meshSize[2] * nbrGhost_d + tv.meshSize[2]*0.5
        x_primal = tv.meshSize[0]*np.arange(layout.allocSize(interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
        y_primal = tv.meshSize[1]*np.arange(layout.allocSize(interpOrder, 'primal', tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_p
        z_primal = tv.meshSize[2]*np.arange(layout.allocSize(interpOrder, 'primal', tv.nbrCells[2])) - tv.meshSize[2] * nbrGhost_p

        By = np.tensordot(np.sin(2*np.pi/tv.domainSize[0] * x_dual),
                          np.tensordot(np.cos(2*np.pi/tv.domainSize[1] * y_primal),
                                       np.sin(2*np.pi/tv.domainSize[2] * z_dual), axes=0), axes=0)

        Ez = np.tensordot(np.sin(2*np.pi/tv.domainSize[0] * x_primal),
                          np.tensordot(np.cos(2*np.pi/tv.domainSize[1] * y_primal),
                                       np.sin(2*np.pi/tv.domainSize[2] * z_dual), axes=0), axes=0)

        psi_d_X = layout.physicalStartIndex(interpOrder, 'dual')
        pei_d_X = layout.physicalEndIndex(  interpOrder, 'dual',   tv.nbrCells[0])
        psi_d_Y = layout.physicalStartIndex(interpOrder, 'dual')
        pei_d_Y = layout.physicalEndIndex(  interpOrder, 'dual',   tv.nbrCells[1])
        psi_d_Z = layout.physicalStartIndex(interpOrder, 'dual')
        pei_d_Z = layout.physicalEndIndex(  interpOrder, 'dual',   tv.nbrCells[2])
        psi_p_X = layout.physicalStartIndex(interpOrder, 'primal')
        pei_p_X = layout.physicalEndIndex(  interpOrder, 'primal', tv.nbrCells[0])
        psi_p_Y = layout.physicalStartIndex(interpOrder, 'primal')
        pei_p_Y = layout.physicalEndIndex(  interpOrder, 'primal', tv.nbrCells[1])
        psi_p_Z = layout.physicalStartIndex(interpOrder, 'primal')
        pei_p_Z = layout.physicalEndIndex(  interpOrder, 'primal', tv.nbrCells[2])

        dxBy[psi_p_X:pei_p_X+1,:,:]  = (By[psi_d_X  :pei_d_X+2,:,:] - By[psi_d_X-1:pei_d_X+1,:,:])/tv.meshSize[0]
        dxEz[psi_d_X:pei_d_X+1,:,:]  = (Ez[psi_p_X+1:pei_p_X+1,:,:] - Ez[psi_p_X  :pei_p_X  ,:,:])/tv.meshSize[0]
        dyBy[:,psi_d_Y:pei_d_Y+1,:]  = (By[:,psi_p_Y+1:pei_p_Y+1,:] - By[:,psi_p_Y  :pei_p_Y  ,:])/tv.meshSize[1]
        dyEz[:,psi_d_Y:pei_d_Y+1,:]  = (Ez[:,psi_p_Y+1:pei_p_Y+1,:] - Ez[:,psi_p_Y  :pei_p_Y  ,:])/tv.meshSize[1]
        dzBy[:,:,psi_p_Z:pei_p_Z+1]  = (By[:,:,psi_d_Z  :pei_d_Z+2] - By[:,:,psi_d_Z-1:pei_d_Z+1])/tv.meshSize[2]
        dzEz[:,:,psi_p_Z:pei_p_Z+1]  = (Ez[:,:,psi_d_Z  :pei_d_Z+2] - Ez[:,:,psi_d_Z-1:pei_d_Z+1])/tv.meshSize[2]

        np.savetxt(os.path.join(path,filename_dxBy), dxBy.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dyBy), dyBy.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dzBy), dzBy.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dxEz), dxEz.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dyEz), dyEz.flatten('C'), delimiter=" ")
        np.savetxt(os.path.join(path,filename_dzEz), dzEz.flatten('C'), delimiter=" ")


def main(path='./'):

    test_deriv1D(path)
    test_deriv2D(path)
    test_deriv3D(path)


if __name__ == '__main__':
    main()

