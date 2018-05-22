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

def test_deriv1D(path):
    # derivative along X:
    # By is dual in X, and derived along X it should lie on Ez
    # Ez is primal in X, and derived along X it should lie on By

    layout = gridlayout.GridLayout() # yee layout

    nbrCellsX = 50
    meshSize  = 0.1
    interpOrders = [1,2,3]
    ByCentering = 'dual'
    EzCentering = 'primal'
    x_max = 5.

    for interpOrder in interpOrders:

        filename_dxBy = 'dxBy_interpOrder_{}.txt'.format(interpOrder)
        filename_dxEz = 'dxEz_interpOrder_{}.txt'.format(interpOrder)

        By   = np.arange(layout.allocSize(interpOrder, ByCentering, nbrCellsX))
        Ez   = np.arange(layout.allocSize(interpOrder, EzCentering, nbrCellsX))

        dxBy = np.arange(layout.allocSizeDerived(interpOrder, ByCentering, nbrCellsX), dtype=np.float64)
        dxEz = np.arange(layout.allocSizeDerived(interpOrder, EzCentering, nbrCellsX), dtype=np.float64)

        nbrGhost_p = layout.nbrGhosts(interpOrder, 'primal')
        nbrGhost_d = layout.nbrGhosts(interpOrder, 'dual')

        x_dual   = meshSize*np.arange(layout.allocSize(interpOrder, 'dual', nbrCellsX)) - meshSize * nbrGhost_d + meshSize/2.
        x_primal = meshSize*np.arange(layout.allocSize(interpOrder, 'primal', nbrCellsX)) - meshSize * nbrGhost_p

        psi_d = layout.physicalStartIndex(interpOrder,'dual')
        pei_d = layout.physicalEndIndex(interpOrder, 'dual',nbrCellsX)
        psi_p = layout.physicalStartIndex(interpOrder, 'primal')
        pei_p = layout.physicalEndIndex(interpOrder, 'primal', nbrCellsX)

        By = np.cos(2*np.pi/x_max * x_dual)
        Ez = np.cos(2*np.pi/x_max * x_primal)

        dxBy[psi_p:pei_p+1]  = (By[psi_d:pei_d+2] - By[psi_d-1:pei_d+1])/meshSize
        dxEz[psi_d:pei_d+1]  = (Ez[psi_p+1:pei_p+1] - Ez[psi_p:pei_p])/meshSize

        np.savetxt(os.path.join(path,filename_dxBy), dxBy, delimiter=" ")
        np.savetxt(os.path.join(path,filename_dxEz), dxEz, delimiter=" ")

        return By,Ez,dxBy, dxEz, x_primal, x_dual, psi_d, pei_d, psi_p, pei_p


def main(path='./'):

    test_deriv1D(path)


if __name__ == '__main__':
    main()

