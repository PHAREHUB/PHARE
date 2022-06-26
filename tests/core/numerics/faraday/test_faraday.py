
# the goal of this scrip is to generate a file for the faraday unit test
# the file will contain expected values for the unit testself.

# the code consists in calculating the curl of a function numerically

import os
import sys
import numpy as np

from pyphare.core import gridlayout

class TestVariables (object) :

    def __init__(self) :

        """ This class only has a constructor
            in order to set the variables to test the derivatives
            all quantities are 3D

        """

        self.nbrCells = (50, 30, 40)
        self.meshSize = (0.1, 0.2, 0.3)
        self.interpOrder = 1
        self.BxCentering = ('primal', 'dual', 'dual')
        self.ByCentering = ('dual', 'primal', 'dual')
        self.BzCentering = ('dual', 'dual', 'primal')
        self.ExCentering = ('dual', 'primal', 'primal')
        self.EyCentering = ('primal', 'dual', 'primal')
        self.EzCentering = ('primal', 'primal', 'dual')
        self.domainSize = tuple(n * m for n, m in zip(self.nbrCells, self.meshSize))
        self.dt = 1.0


# ------------------------------------------------------------------------------
# since the derivative is tested for all interporders
# we don't do ampere for all interporders, only for interporder=1
# ------------------------------------------------------------------------------

def test_faraday_yee1D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()


    By = np.zeros(layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]), dtype=np.float64)
    Bz = np.zeros(layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]), dtype=np.float64)

    Ey = np.zeros(layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]), dtype=np.float64)
    Ez = np.zeros(layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]), dtype=np.float64)


    ByNew = np.zeros(layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]), dtype=np.float64)
    BzNew = np.zeros(layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]), dtype=np.float64)

    psi_p_X = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_X = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[0])

    psi_d_X = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_X = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[0])

    nbrGhost_p = layout.nbrGhosts(tv.interpOrder, 'primal')
    nbrGhost_d = layout.nbrGhosts(tv.interpOrder, 'dual')

    x_primal = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
    x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5

    Ey = np.cos(2*np.pi/tv.domainSize[0] * x_primal)
    Ez = np.sin(2*np.pi/tv.domainSize[0] * x_primal)
    By = np.tanh(x_dual - 0.5*tv.domainSize[0])
    Bz = np.tanh(x_dual - 0.5*tv.domainSize[0])

    ByNew[psi_d_X:pei_d_X+1] = By[psi_d_X:pei_d_X+1] + tv.dt*( Ez[psi_p_X+1:pei_p_X+1] - Ez[psi_p_X:pei_p_X])/tv.meshSize[0]
    BzNew[psi_d_X:pei_d_X+1] = Bz[psi_d_X:pei_d_X+1] - tv.dt*( Ey[psi_p_X+1:pei_p_X+1] - Ey[psi_p_X:pei_p_X])/tv.meshSize[0]

    filename_dbydt = "dbydt_yee_1D_order1.txt"
    filename_dbzdt = "dbzdt_yee_1D_order1.txt"

    np.savetxt(os.path.join(path, filename_dbydt), ByNew, delimiter=" ")
    np.savetxt(os.path.join(path, filename_dbzdt), BzNew, delimiter=" ")


def test_faraday_yee2D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    Bx = np.zeros([layout.allocSize(tv.interpOrder, tv.BxCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.BxCentering[1], tv.nbrCells[1])], dtype=np.float64)
    By = np.zeros([layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.ByCentering[1], tv.nbrCells[1])], dtype=np.float64)
    Bz = np.zeros([layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.BzCentering[1], tv.nbrCells[1])], dtype=np.float64)

    Ex = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1])], dtype=np.float64)
    Ey = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1])], dtype=np.float64)
    Ez = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)

    BxNew = np.zeros([layout.allocSize(tv.interpOrder, tv.BxCentering[0], tv.nbrCells[0]),
                      layout.allocSize(tv.interpOrder, tv.BxCentering[1], tv.nbrCells[1])], dtype=np.float64)
    ByNew = np.zeros([layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                      layout.allocSize(tv.interpOrder, tv.ByCentering[1], tv.nbrCells[1])], dtype=np.float64)
    BzNew = np.zeros([layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]),
                      layout.allocSize(tv.interpOrder, tv.BzCentering[1], tv.nbrCells[1])], dtype=np.float64)
    w1 = np.zeros_like(BzNew)
    w2 = np.zeros_like(BzNew)

    psi_p_X = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_X = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[0])
    psi_p_Y = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_Y = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[1])

    psi_d_X = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_X = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[0])
    psi_d_Y = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_Y = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[1])

    nbrGhost_p = layout.nbrGhosts(tv.interpOrder, 'primal')
    nbrGhost_d = layout.nbrGhosts(tv.interpOrder, 'dual')

    x_primal = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
    y_primal = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_p
    x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
    y_dual   = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_d + tv.meshSize[1]*0.5

    Ex = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_dual),
                      np.sin(2*np.pi/tv.domainSize[1] * y_primal), axes=0)
    Ey = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_primal),
                     np.tanh(2*np.pi/tv.domainSize[1] * y_dual), axes=0)
    Ez = np.tensordot(np.sin(2*np.pi/tv.domainSize[0] * x_primal),
                     np.tanh(2*np.pi/tv.domainSize[1] * y_primal), axes=0)

    Bx = np.tensordot(np.tanh(x_primal - 0.5*tv.domainSize[0]),
                      np.tanh(y_dual - 0.5*tv.domainSize[1]), axes=0)
    By = np.tensordot(np.tanh(x_dual - 0.5*tv.domainSize[0]),
                      np.tanh(y_primal - 0.5*tv.domainSize[1]), axes=0)
    Bz = np.tensordot(np.tanh(x_dual - 0.5*tv.domainSize[0]),
                      np.tanh(y_dual - 0.5*tv.domainSize[1]), axes=0)

    BxNew[:,psi_d_Y:pei_d_Y+1] = Bx[:,psi_d_Y:pei_d_Y+1] - tv.dt*( Ez[:,psi_p_Y+1:pei_p_Y+1] - Ez[:,psi_p_Y:pei_p_Y] )/tv.meshSize[1]
    ByNew[psi_d_X:pei_d_X+1,:] = By[psi_d_X:pei_d_X+1,:] + tv.dt*( Ez[psi_p_X+1:pei_p_X+1,:] - Ez[psi_p_X:pei_p_X,:] )/tv.meshSize[0]
    w1[psi_d_X:pei_d_X+1,:] = - tv.dt*( Ey[psi_p_X+1:pei_p_X+1,:] - Ey[psi_p_X:pei_p_X,:])/tv.meshSize[0]
    w2[:,psi_d_Y:pei_d_Y+1] = + tv.dt*( Ex[:,psi_p_Y+1:pei_p_Y+1] - Ex[:,psi_p_Y:pei_p_Y])/tv.meshSize[1]
    BzNew = Bz + w1 + w2

    filename_dbxdt = "dbxdt_yee_2D_order1.txt"
    filename_dbydt = "dbydt_yee_2D_order1.txt"
    filename_dbzdt = "dbzdt_yee_2D_order1.txt"

    np.savetxt(os.path.join(path, filename_dbxdt), BxNew.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_dbydt), ByNew.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_dbzdt), BzNew.flatten('C'), delimiter=" ")


def test_faraday_yee3D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    Bx = np.zeros([layout.allocSize(tv.interpOrder, tv.BxCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.BxCentering[1], tv.nbrCells[1]),
                   layout.allocSize(tv.interpOrder, tv.BxCentering[2], tv.nbrCells[2])], dtype=np.float64)
    By = np.zeros([layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.ByCentering[1], tv.nbrCells[1]),
                   layout.allocSize(tv.interpOrder, tv.ByCentering[2], tv.nbrCells[2])], dtype=np.float64)
    Bz = np.zeros([layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.BzCentering[1], tv.nbrCells[1]),
                   layout.allocSize(tv.interpOrder, tv.BzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    Ex = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1]),
                   layout.allocSize(tv.interpOrder, tv.ExCentering[2], tv.nbrCells[2])], dtype=np.float64)
    Ey = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1]),
                   layout.allocSize(tv.interpOrder, tv.EyCentering[2], tv.nbrCells[2])], dtype=np.float64)
    Ez = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                   layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                   layout.allocSize(tv.interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    BxNew = np.zeros([layout.allocSize(tv.interpOrder, tv.BxCentering[0], tv.nbrCells[0]),
                      layout.allocSize(tv.interpOrder, tv.BxCentering[1], tv.nbrCells[1]),
                      layout.allocSize(tv.interpOrder, tv.BxCentering[2], tv.nbrCells[2])], dtype=np.float64)
    ByNew = np.zeros([layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                      layout.allocSize(tv.interpOrder, tv.ByCentering[1], tv.nbrCells[1]),
                      layout.allocSize(tv.interpOrder, tv.ByCentering[2], tv.nbrCells[2])], dtype=np.float64)
    BzNew = np.zeros([layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]),
                      layout.allocSize(tv.interpOrder, tv.BzCentering[1], tv.nbrCells[1]),
                      layout.allocSize(tv.interpOrder, tv.BzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    u1 = np.zeros_like(BxNew)
    u2 = np.zeros_like(BxNew)
    v1 = np.zeros_like(ByNew)
    v2 = np.zeros_like(ByNew)
    w1 = np.zeros_like(BzNew)
    w2 = np.zeros_like(BzNew)

    psi_p_X = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_X = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[0])
    psi_p_Y = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_Y = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[1])
    psi_p_Z = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_Z = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[2])

    psi_d_X = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_X = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[0])
    psi_d_Y = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_Y = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[1])
    psi_d_Z = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_Z = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[2])

    nbrGhost_p = layout.nbrGhosts(tv.interpOrder, 'primal')
    nbrGhost_d = layout.nbrGhosts(tv.interpOrder, 'dual')

    x_primal = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
    y_primal = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_p
    z_primal = tv.meshSize[2]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[2])) - tv.meshSize[2] * nbrGhost_p
    x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
    y_dual   = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_d + tv.meshSize[1]*0.5
    z_dual   = tv.meshSize[2]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[2])) - tv.meshSize[2] * nbrGhost_d + tv.meshSize[2]*0.5

    Ex = np.tensordot(np.sin(2*np.pi/tv.domainSize[0] * x_dual),
                      np.tensordot(np.cos(2*np.pi/tv.domainSize[1] * y_primal),
                                   np.tanh(2*np.pi/tv.domainSize[2] * z_primal), axes=0), axes=0)
    Ey = np.tensordot(np.tanh(2*np.pi/tv.domainSize[0] * x_primal),
                      np.tensordot(np.sin(2*np.pi/tv.domainSize[1] * y_dual),
                                   np.cos(2*np.pi/tv.domainSize[2] * z_primal), axes=0), axes=0)
    Ez = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_primal),
                      np.tensordot(np.tanh(2*np.pi/tv.domainSize[1] * y_primal),
                                   np.sin(2*np.pi/tv.domainSize[2] * z_dual), axes=0), axes=0)

    Bx = np.tensordot(np.tanh(x_primal - 0.5*tv.domainSize[0]),
                      np.tensordot(np.tanh(y_dual - 0.5*tv.domainSize[1]),
                                   np.tanh(z_dual - 0.5*tv.domainSize[2]), axes=0), axes=0)
    By = np.tensordot(np.tanh(x_dual - 0.5*tv.domainSize[0]),
                      np.tensordot(np.tanh(y_primal - 0.5*tv.domainSize[1]),
                                   np.tanh(z_dual - 0.5*tv.domainSize[2]), axes=0), axes=0)
    Bz = np.tensordot(np.tanh(x_dual - 0.5*tv.domainSize[0]),
                      np.tensordot(np.tanh(y_dual - 0.5*tv.domainSize[1]),
                                   np.tanh(z_primal - 0.5*tv.domainSize[2]), axes=0), axes=0)

    u1[:,psi_d_Y:pei_d_Y+1,:] = - tv.dt*( Ez[:,psi_p_Y+1:pei_p_Y+1,:] - Ez[:,psi_p_Y:pei_p_Y,:] )/tv.meshSize[1]
    u2[:,:,psi_d_Z:pei_d_Z+1] = + tv.dt*( Ey[:,:,psi_p_Z+1:pei_p_Z+1] - Ey[:,:,psi_p_Z:pei_p_Z] )/tv.meshSize[2]
    BxNew = Bx + u1 + u2

    v1[:,:,psi_d_Z:pei_d_Z+1] = - tv.dt*( Ex[:,:,psi_p_Z+1:pei_p_Z+1] - Ex[:,:,psi_p_Z:pei_p_Z] )/tv.meshSize[2]
    v2[psi_d_X:pei_d_X+1,:,:] = + tv.dt*( Ez[psi_p_X+1:pei_p_X+1,:,:] - Ez[psi_p_X:pei_p_X,:,:] )/tv.meshSize[0]
    ByNew = By + v1 + v2

    w1[psi_d_X:pei_d_X+1,:,:] = - tv.dt*( Ey[psi_p_X+1:pei_p_X+1,:,:] - Ey[psi_p_X:pei_p_X,:,:] )/tv.meshSize[0]
    w2[:,psi_d_Y:pei_d_Y+1,:] = + tv.dt*( Ex[:,psi_p_Y+1:pei_p_Y+1,:] - Ex[:,psi_p_Y:pei_p_Y,:] )/tv.meshSize[1]
    BzNew = Bz + w1 + w2

    filename_dbxdt = "dbxdt_yee_3D_order1.txt"
    filename_dbydt = "dbydt_yee_3D_order1.txt"
    filename_dbzdt = "dbzdt_yee_3D_order1.txt"

    np.savetxt(os.path.join(path, filename_dbxdt), BxNew.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_dbydt), ByNew.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_dbzdt), BzNew.flatten('C'), delimiter=" ")


def main(path='./'):

    if len(sys.argv) >1:
        path=sys.argv[1]

    test_faraday_yee1D(path)
    test_faraday_yee2D(path)
    test_faraday_yee3D(path)


if __name__ == '__main__':
    main()
