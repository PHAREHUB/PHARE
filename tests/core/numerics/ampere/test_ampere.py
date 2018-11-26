
# the goal of this script is to generate a file for the ampere unit test
# the file will contain expected values for the unit test.

# the code consists in calculating the curl of a function numerically

import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../data/', 'gridlayout'))


import numpy as np
import gridlayout

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
        self.JxCentering = ('dual', 'primal', 'primal')
        self.JyCentering = ('primal', 'dual', 'primal')
        self.JzCentering = ('primal', 'primal', 'dual')
        self.domainSize = tuple(n * m for n, m in zip(self.nbrCells, self.meshSize))


# ------------------------------------------------------------------------------
# since the derivative is tested for all interporders
# we don't do ampere for all interporders, only for interporder=1
# ------------------------------------------------------------------------------

def test_ampere_yee1D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    By = np.zeros(layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]), dtype=np.float64)
    Bz = np.zeros(layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]), dtype=np.float64)
    Jy = np.zeros(layout.allocSize(tv.interpOrder, tv.JyCentering[0], tv.nbrCells[0]), dtype=np.float64)
    Jz = np.zeros(layout.allocSize(tv.interpOrder, tv.JzCentering[0], tv.nbrCells[0]), dtype=np.float64)


    psi_p_X = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_X = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[0])

    psi_d_X = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_X = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[0])

    nbrGhost_p = layout.nbrGhosts(tv.interpOrder, 'primal')
    nbrGhost_d = layout.nbrGhosts(tv.interpOrder, 'dual')

    x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
    x_primal = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p

    By = np.cos(2*np.pi/tv.domainSize[0] * x_dual)
    Bz = np.sin(2*np.pi/tv.domainSize[0] * x_dual)

    # Jy = -dxBz
    # Jz =  dxBy
    Jy[psi_p_X:pei_p_X+1] = -(Bz[psi_d_X:pei_d_X+2] - Bz[psi_d_X-1:pei_d_X+1])/tv.meshSize[0]
    Jz[psi_p_X:pei_p_X+1] =  (By[psi_d_X:pei_d_X+2] - By[psi_d_X-1:pei_d_X+1])/tv.meshSize[0]

    filename_jy = "jy_yee_1D_order1.txt"
    filename_jz = "jz_yee_1D_order1.txt"

    np.savetxt(os.path.join(path, filename_jy), Jy, delimiter=" ")
    np.savetxt(os.path.join(path, filename_jz), Jz, delimiter=" ")


def test_ampere_yee2D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    Bx   = np.zeros([layout.allocSize(tv.interpOrder, tv.BxCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.BxCentering[1], tv.nbrCells[1])], dtype=np.float64)
    By   = np.zeros([layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.ByCentering[1], tv.nbrCells[1])], dtype=np.float64)
    Bz   = np.zeros([layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.BzCentering[1], tv.nbrCells[1])], dtype=np.float64)

    Jx   = np.zeros([layout.allocSize(tv.interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.JxCentering[1], tv.nbrCells[1])], dtype=np.float64)
    Jy   = np.zeros([layout.allocSize(tv.interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.JyCentering[1], tv.nbrCells[1])], dtype=np.float64)
    Jz   = np.zeros([layout.allocSize(tv.interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.JzCentering[1], tv.nbrCells[1])], dtype=np.float64)
    w1 = np.zeros_like(Jz)
    w2 = np.zeros_like(Jz)

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

    x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
    y_dual   = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_d + tv.meshSize[1]*0.5
    x_primal = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
    y_primal = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_p

    Bx = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_primal),
                      np.sin(2*np.pi/tv.domainSize[1] * y_dual), axes=0)
    By = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_dual),
                     np.tanh(2*np.pi/tv.domainSize[1] * y_primal), axes=0)
    Bz = np.tensordot(np.sin(2*np.pi/tv.domainSize[0] * x_dual),
                     np.tanh(2*np.pi/tv.domainSize[1] * y_dual), axes=0)

    # Jx =  dyBz
    # Jy = -dxBz
    # Jz =  dxBy - dyBx
    Jx[:,psi_p_Y:pei_p_Y+1] =  (Bz[:,psi_d_Y:pei_d_Y+2] - Bz[:,psi_d_Y-1:pei_d_Y+1])/tv.meshSize[1]
    Jy[psi_p_X:pei_p_X+1,:] = -(Bz[psi_d_X:pei_d_X+2,:] - Bz[psi_d_X-1:pei_d_X+1,:])/tv.meshSize[0]
    w1[psi_p_X:pei_p_X+1,:] =  (By[psi_d_X:pei_d_X+2,:] - By[psi_d_X-1:pei_d_X+1,:])/tv.meshSize[0]
    w2[:,psi_p_Y:pei_p_Y+1] = -(Bx[:,psi_d_Y:pei_d_Y+2] - Bx[:,psi_d_Y-1:pei_d_Y+1])/tv.meshSize[1]
    Jz = w1+w2

    filename_jx = "jx_yee_2D_order1.txt"
    filename_jy = "jy_yee_2D_order1.txt"
    filename_jz = "jz_yee_2D_order1.txt"

    np.savetxt(os.path.join(path, filename_jx), Jx.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_jy), Jy.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_jz), Jz.flatten('C'), delimiter=" ")


def test_ampere_yee3D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()

    Bx   = np.zeros([layout.allocSize(tv.interpOrder, tv.BxCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.BxCentering[1], tv.nbrCells[1]),
                     layout.allocSize(tv.interpOrder, tv.BxCentering[2], tv.nbrCells[2])], dtype=np.float64)
    By   = np.zeros([layout.allocSize(tv.interpOrder, tv.ByCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.ByCentering[1], tv.nbrCells[1]),
                     layout.allocSize(tv.interpOrder, tv.ByCentering[2], tv.nbrCells[2])], dtype=np.float64)
    Bz   = np.zeros([layout.allocSize(tv.interpOrder, tv.BzCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.BzCentering[1], tv.nbrCells[1]),
                     layout.allocSize(tv.interpOrder, tv.BzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    Jx   = np.zeros([layout.allocSize(tv.interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.JxCentering[1], tv.nbrCells[1]),
                     layout.allocSize(tv.interpOrder, tv.JxCentering[2], tv.nbrCells[2])], dtype=np.float64)
    Jy   = np.zeros([layout.allocSize(tv.interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.JyCentering[1], tv.nbrCells[1]),
                     layout.allocSize(tv.interpOrder, tv.JyCentering[2], tv.nbrCells[2])], dtype=np.float64)
    Jz   = np.zeros([layout.allocSize(tv.interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
                     layout.allocSize(tv.interpOrder, tv.JzCentering[1], tv.nbrCells[1]),
                     layout.allocSize(tv.interpOrder, tv.JzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    u1 = np.zeros_like(Jx)
    u2 = np.zeros_like(Jx)
    v1 = np.zeros_like(Jy)
    v2 = np.zeros_like(Jy)
    w1 = np.zeros_like(Jz)
    w2 = np.zeros_like(Jz)

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

    x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5
    y_dual   = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_d + tv.meshSize[1]*0.5
    z_dual   = tv.meshSize[2]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[2])) - tv.meshSize[2] * nbrGhost_d + tv.meshSize[2]*0.5
    x_primal = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
    y_primal = tv.meshSize[1]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[1])) - tv.meshSize[1] * nbrGhost_p
    z_primal = tv.meshSize[2]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[2])) - tv.meshSize[2] * nbrGhost_p

    Bx = np.tensordot(np.sin(2*np.pi/tv.domainSize[0] * x_primal),
                      np.tensordot(np.cos(2*np.pi/tv.domainSize[1] * y_dual),
                                   np.tanh(2*np.pi/tv.domainSize[2] * z_dual), axes=0), axes=0)
    By = np.tensordot(np.tanh(2*np.pi/tv.domainSize[0] * x_dual),
                      np.tensordot(np.sin(2*np.pi/tv.domainSize[1] * y_primal),
                                   np.cos(2*np.pi/tv.domainSize[2] * z_dual), axes=0), axes=0)
    Bz = np.tensordot(np.cos(2*np.pi/tv.domainSize[0] * x_dual),
                      np.tensordot(np.tanh(2*np.pi/tv.domainSize[1] * y_dual),
                                   np.sin(2*np.pi/tv.domainSize[2] * z_primal), axes=0), axes=0)

    # Jx = dyBz - dzBy
    # Jy = dzBx - dxBz
    # Jz = dxBy - dyBx
    u1[:,psi_p_Y:pei_p_Y+1,:] =  (Bz[:,psi_d_Y:pei_d_Y+2,:] - Bz[:,psi_d_Y-1:pei_d_Y+1,:])/tv.meshSize[1]
    u2[:,:,psi_p_Z:pei_p_Z+1] = -(By[:,:,psi_d_Z:pei_d_Z+2] - By[:,:,psi_d_Z-1:pei_d_Z+1])/tv.meshSize[2]
    Jx = u1+u2
    v1[:,:,psi_p_Z:pei_p_Z+1] =  (Bx[:,:,psi_d_Z:pei_d_Z+2] - Bx[:,:,psi_d_Z-1:pei_d_Z+1])/tv.meshSize[2]
    v2[psi_p_X:pei_p_X+1,:,:] = -(Bz[psi_d_X:pei_d_X+2,:,:] - Bz[psi_d_X-1:pei_d_X+1,:,:])/tv.meshSize[0]
    Jy = v1+v2
    w1[psi_p_X:pei_p_X+1,:,:] =  (By[psi_d_X:pei_d_X+2,:,:] - By[psi_d_X-1:pei_d_X+1,:,:])/tv.meshSize[0]
    w2[:,psi_p_Y:pei_p_Y+1,:] = -(Bx[:,psi_d_Y:pei_d_Y+2,:] - Bx[:,psi_d_Y-1:pei_d_Y+1,:])/tv.meshSize[1]
    Jz = w1+w2

    filename_jx = "jx_yee_3D_order1.txt"
    filename_jy = "jy_yee_3D_order1.txt"
    filename_jz = "jz_yee_3D_order1.txt"

    np.savetxt(os.path.join(path, filename_jx), Jx.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_jy), Jy.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_jz), Jz.flatten('C'), delimiter=" ")


def main(path='./'):

    if len(sys.argv) >1:
        path=sys.argv[1]

    test_ampere_yee1D(path)
    test_ampere_yee2D(path)
    test_ampere_yee3D(path)


if __name__ == '__main__':
    main()
