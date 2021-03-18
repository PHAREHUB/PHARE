
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

        self.nbrCells = (50, 40, 30)
        self.meshSize = (0.1, 0.2, 0.3)
        self.interpOrder = 1
        self.BxCentering = ('primal', 'dual', 'dual')
        self.ByCentering = ('dual', 'primal', 'dual')
        self.BzCentering = ('dual', 'dual', 'primal')
        self.JxCentering = ('dual', 'primal', 'primal')
        self.JyCentering = ('primal', 'dual', 'primal')
        self.JzCentering = ('primal', 'primal', 'dual')
        self.ExCentering = ('dual', 'primal', 'primal')
        self.EyCentering = ('primal', 'dual', 'primal')
        self.EzCentering = ('primal', 'primal', 'dual')
        self.MomentsCentering = ('primal', 'primal', 'primal')
        self.domainSize = tuple(n * m for n, m in zip(self.nbrCells, self.meshSize))


# ------------------------------------------------------------------------------
# since the derivative is tested for all interporders
# we don't do ohm for all interporders, only for interporder=1
# ------------------------------------------------------------------------------

def test_ohm_yee1D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()
    eta = 1.0
    nu = 0.01

    idealx = np.zeros(layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]), dtype=np.float64)
    idealy = np.zeros(layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]), dtype=np.float64)
    idealz = np.zeros(layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]), dtype=np.float64)

    pressx = np.zeros(layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]), dtype=np.float64)
    pressy = np.zeros(layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]), dtype=np.float64)
    pressz = np.zeros(layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]), dtype=np.float64)

    resistx = np.zeros(layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]), dtype=np.float64)
    resisty = np.zeros(layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]), dtype=np.float64)
    resistz = np.zeros(layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]), dtype=np.float64)

    viscousx = np.zeros(layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]), dtype=np.float64)
    viscousy = np.zeros(layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]), dtype=np.float64)
    viscousz = np.zeros(layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]), dtype=np.float64)

    psi_p_X = layout.physicalStartIndex(tv.interpOrder, 'primal')
    pei_p_X = layout.physicalEndIndex(  tv.interpOrder, 'primal', tv.nbrCells[0])

    psi_d_X = layout.physicalStartIndex(tv.interpOrder, 'dual')
    pei_d_X = layout.physicalEndIndex(  tv.interpOrder, 'dual',   tv.nbrCells[0])

    nbrGhost_p = layout.nbrGhosts(tv.interpOrder, 'primal')
    nbrGhost_d = layout.nbrGhosts(tv.interpOrder, 'dual')

    x_primal = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'primal', tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_p
    x_dual   = tv.meshSize[0]*np.arange(layout.allocSize(tv.interpOrder, 'dual'  , tv.nbrCells[0])) - tv.meshSize[0] * nbrGhost_d + tv.meshSize[0]*0.5

    # analytical profiles of density, velocity...
    n  = np.cosh(0.5*x_primal)
    Vx = np.sinh(0.2*x_primal)
    Vy = np.sinh(0.3*x_primal)
    Vz = np.sinh(0.4*x_primal)
    P  = np.cosh(0.5*x_primal)
    Bx = np.cosh(0.2*x_primal)
    By = np.cosh(0.3*x_dual)
    Bz = np.cosh(0.4*x_dual)
    Jx = np.tanh(0.2*x_dual)
    Jy = np.tanh(0.3*x_primal)
    Jz = np.tanh(0.4*x_primal)

    # ideal term
    idealx[psi_d_X:pei_d_X+1] = 0.5*(Vz[psi_p_X:pei_p_X]+Vz[psi_p_X+1:pei_p_X+1])*By[psi_d_X:pei_d_X+1]\
                              - 0.5*(Vy[psi_p_X:pei_p_X]+Vy[psi_p_X+1:pei_p_X+1])*Bz[psi_d_X:pei_d_X+1]

    idealy[psi_p_X:pei_p_X+1] = 0.5*Vx[psi_p_X:pei_p_X+1]*(Bz[psi_d_X-1:pei_d_X+1]+Bz[psi_d_X:pei_d_X+2])\
                              - Vz[psi_p_X:pei_p_X+1]*Bx[psi_p_X:pei_p_X+1]

    idealz[psi_p_X:pei_p_X+1] = Vy[psi_p_X:pei_p_X+1]*Bx[psi_p_X:pei_p_X+1]\
                              - 0.5*Vx[psi_p_X:pei_p_X+1]*(By[psi_d_X-1:pei_d_X+1]+By[psi_d_X:pei_d_X+2])

    # pressure term
    pressx[psi_d_X:pei_d_X+1] = ((P[psi_p_X+1:pei_p_X+1]-P[psi_p_X:pei_p_X])/tv.meshSize[0])\
                                /(0.5*(n[psi_p_X:pei_p_X]+n[psi_p_X+1:pei_p_X+1]))
    pressy[psi_p_X:pei_p_X+1] = 0
    pressz[psi_p_X:pei_p_X+1] = 0

    # resistive term
    resistx[psi_d_X:pei_d_X+1] = eta*Jx[psi_d_X:pei_d_X+1]
    resisty[psi_p_X:pei_p_X+1] = eta*Jy[psi_p_X:pei_p_X+1]
    resistz[psi_p_X:pei_p_X+1] = eta*Jz[psi_p_X:pei_p_X+1]

    # viscous term
    viscousx[psi_d_X:pei_d_X+1] = -nu*(Jx[psi_d_X-1:pei_d_X]-2.0*Jx[psi_d_X:pei_d_X+1]+Jx[psi_d_X+1:pei_d_X+2])/(tv.meshSize[0]*tv.meshSize[0])
    viscousy[psi_p_X:pei_p_X+1] = -nu*(Jy[psi_p_X-1:pei_p_X]-2.0*Jy[psi_p_X:pei_p_X+1]+Jy[psi_p_X+1:pei_p_X+2])/(tv.meshSize[0]*tv.meshSize[0])
    viscousz[psi_p_X:pei_p_X+1] = -nu*(Jz[psi_p_X-1:pei_p_X]-2.0*Jz[psi_p_X:pei_p_X+1]+Jz[psi_p_X+1:pei_p_X+2])/(tv.meshSize[0]*tv.meshSize[0])

    ExNew = idealx+pressx+resistx+viscousx
    EyNew = idealy+pressy+resisty+viscousy
    EzNew = idealz+pressz+resistz+viscousz

    filename_ohmx = "ohmx_yee_1D_order1.txt"
    filename_ohmy = "ohmy_yee_1D_order1.txt"
    filename_ohmz = "ohmz_yee_1D_order1.txt"

    np.savetxt(os.path.join(path, filename_ohmx), ExNew, delimiter=" ")
    np.savetxt(os.path.join(path, filename_ohmy), EyNew, delimiter=" ")
    np.savetxt(os.path.join(path, filename_ohmz), EzNew, delimiter=" ")



def test_ohm_yee2D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()
    eta = 1.0
    nu = 0.01

    idealx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1])], dtype=np.float64)
    idealy = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1])], dtype=np.float64)
    idealz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)

    pressx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1])], dtype=np.float64)
    pressy = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1])], dtype=np.float64)
    pressz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)

    resistx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                        layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1])], dtype=np.float64)
    resisty = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                        layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1])], dtype=np.float64)
    resistz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                        layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)

    viscousx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                         layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1])], dtype=np.float64)
    viscousy = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                         layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1])], dtype=np.float64)
    viscousz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1])], dtype=np.float64)

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

    # analytical profiles of density, velocity...
    n  = np.tensordot(np.cosh(0.5*x_primal),
                      np.cosh(0.5*y_primal), axes=0)
    Vx = np.tensordot(np.sinh(0.2*x_primal),
                      np.sinh(0.2*y_primal), axes=0)
    Vy = np.tensordot(np.sinh(0.3*x_primal),
                      np.sinh(0.3*y_primal), axes=0)
    Vz = np.tensordot(np.sinh(0.4*x_primal),
                      np.sinh(0.4*y_primal), axes=0)
    P  = np.tensordot(np.cosh(0.5*x_primal),
                      np.cosh(0.5*y_primal), axes=0)
    Bx = np.tensordot(np.cosh(0.2*x_primal),
                      np.cosh(0.2*y_dual), axes=0)
    By = np.tensordot(np.cosh(0.3*x_dual),
                      np.cosh(0.3*y_primal), axes=0)
    Bz = np.tensordot(np.cosh(0.4*x_dual),
                      np.cosh(0.4*y_dual), axes=0)
    Jx = np.tensordot(np.tanh(0.2*x_dual),
                      np.tanh(0.2*y_primal), axes=0)
    Jy = np.tensordot(np.tanh(0.3*x_primal),
                      np.tanh(0.3*y_dual), axes=0)
    Jz = np.tensordot(np.tanh(0.4*x_primal),
                      np.tanh(0.4*y_primal), axes=0)

    # ideal term
    idealx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1] = \
 -0.25*(Vy[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1]\
       +Vy[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1])\
      *(Bz[psi_d_X  :pei_d_X+1, psi_d_Y-1:pei_d_Y+1]\
       +Bz[psi_d_X  :pei_d_X+1, psi_d_Y  :pei_d_Y+2])\
 +0.5 *(Vz[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1]\
       +Vz[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1])\
      * By[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1]

    idealy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1] = \
 -0.5 *(Vz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  ]\
       +Vz[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1])\
      *(Bx[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1])\
 +0.25*(Vx[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  ]\
       +Vx[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1])\
      *(Bz[psi_d_X-1:pei_d_X+1, psi_d_Y  :pei_d_Y+1]\
       +Bz[psi_d_X  :pei_d_X+2, psi_d_Y  :pei_d_Y+1])

    idealz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1] = \
 -0.5  *Vx[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1]\
      *(By[psi_d_X-1:pei_d_X+1, psi_p_Y  :pei_p_Y+1]\
       +By[psi_d_X  :pei_d_X+2, psi_p_Y  :pei_p_Y+1])\
 +0.5  *Vy[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1]\
      *(Bx[psi_p_X  :pei_p_X+1, psi_d_Y-1:pei_d_Y+1]\
       +Bx[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+2])

    # pressure term
    pressx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1] = \
       ((P[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1]\
        -P[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1])/tv.meshSize[0])\
  /(0.5*(n[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1]\
        +n[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1]))

    pressy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1] = \
       ((P[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1]\
        -P[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  ])/tv.meshSize[1])\
  /(0.5*(n[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  ]\
        +n[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1]))

    pressz[psi_p_X:pei_p_X+1, psi_p_Y:pei_p_Y+1] = 0

    # resistive term
    resistx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1] = \
     eta*Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1]

    resisty[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1] = \
     eta*Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1]

    resistz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1] = \
     eta*Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1]

    # viscous term
    viscousx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1] = \
     -nu*(Jx[psi_d_X-1:pei_d_X  , psi_p_Y  :pei_p_Y+1]\
     -2.0*Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1]\
         +Jx[psi_d_X+1:pei_d_X+2, psi_p_Y  :pei_p_Y+1])\
        /(tv.meshSize[0]*tv.meshSize[0])\
     -nu*(Jx[psi_d_X  :pei_d_X+1, psi_p_Y-1:pei_p_Y  ]\
     -2.0*Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1]\
         +Jx[psi_d_X  :pei_d_X+1, psi_p_Y+1:pei_p_Y+2])\
        /(tv.meshSize[1]*tv.meshSize[1])\

    viscousy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1] = \
     -nu*(Jy[psi_p_X-1:pei_p_X  , psi_d_Y  :pei_d_Y+1]\
     -2.0*Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1]\
         +Jy[psi_p_X+1:pei_p_X+2, psi_d_Y  :pei_d_Y+1])\
        /(tv.meshSize[0]*tv.meshSize[0])\
     -nu*(Jy[psi_p_X  :pei_p_X+1, psi_d_Y-1:pei_d_Y  ]\
     -2.0*Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1]\
         +Jy[psi_p_X  :pei_p_X+1, psi_d_Y+1:pei_d_Y+2])\
        /(tv.meshSize[1]*tv.meshSize[1])\

    viscousz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1] = \
     -nu*(Jz[psi_p_X-1:pei_p_X  , psi_p_Y  :pei_p_Y+1]\
     -2.0*Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1]\
         +Jz[psi_p_X+1:pei_p_X+2, psi_p_Y  :pei_p_Y+1])\
        /(tv.meshSize[0]*tv.meshSize[0])\
     -nu*(Jz[psi_p_X  :pei_p_X+1, psi_p_Y-1:pei_p_Y  ]\
     -2.0*Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1]\
         +Jz[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+2])\
        /(tv.meshSize[1]*tv.meshSize[1])\

    ExNew = idealx+pressx+resistx+viscousx
    EyNew = idealy+pressy+resisty+viscousy
    EzNew = idealz+pressz+resistz+viscousz

    filename_ohmx = "ohmx_yee_2D_order1.txt"
    filename_ohmy = "ohmy_yee_2D_order1.txt"
    filename_ohmz = "ohmz_yee_2D_order1.txt"

    np.savetxt(os.path.join(path, filename_ohmx), ExNew, delimiter=" ")
    np.savetxt(os.path.join(path, filename_ohmy), EyNew, delimiter=" ")
    np.savetxt(os.path.join(path, filename_ohmz), EzNew, delimiter=" ")



def test_ohm_yee3D(path):

    layout = gridlayout.GridLayout() # yee layout

    tv = TestVariables()
    eta = 1.0
    nu = 0.01


    idealx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1]),
                       layout.allocSize(tv.interpOrder, tv.ExCentering[2], tv.nbrCells[2])], dtype=np.float64)
    idealy = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1]),
                       layout.allocSize(tv.interpOrder, tv.EyCentering[2], tv.nbrCells[2])], dtype=np.float64)
    idealz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                       layout.allocSize(tv.interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    pressx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1]),
                       layout.allocSize(tv.interpOrder, tv.ExCentering[2], tv.nbrCells[2])], dtype=np.float64)
    pressy = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1]),
                       layout.allocSize(tv.interpOrder, tv.EyCentering[2], tv.nbrCells[2])], dtype=np.float64)
    pressz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                       layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                       layout.allocSize(tv.interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    resistx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                        layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1]),
                        layout.allocSize(tv.interpOrder, tv.ExCentering[2], tv.nbrCells[2])], dtype=np.float64)
    resisty = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                        layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1]),
                        layout.allocSize(tv.interpOrder, tv.EyCentering[2], tv.nbrCells[2])], dtype=np.float64)
    resistz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                        layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                        layout.allocSize(tv.interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)

    viscousx = np.zeros([layout.allocSize(tv.interpOrder, tv.ExCentering[0], tv.nbrCells[0]),
                         layout.allocSize(tv.interpOrder, tv.ExCentering[1], tv.nbrCells[1]),
                         layout.allocSize(tv.interpOrder, tv.ExCentering[2], tv.nbrCells[2])], dtype=np.float64)
    viscousy = np.zeros([layout.allocSize(tv.interpOrder, tv.EyCentering[0], tv.nbrCells[0]),
                         layout.allocSize(tv.interpOrder, tv.EyCentering[1], tv.nbrCells[1]),
                         layout.allocSize(tv.interpOrder, tv.EyCentering[2], tv.nbrCells[2])], dtype=np.float64)
    viscousz = np.zeros([layout.allocSize(tv.interpOrder, tv.EzCentering[0], tv.nbrCells[0]),
                         layout.allocSize(tv.interpOrder, tv.EzCentering[1], tv.nbrCells[1]),
                         layout.allocSize(tv.interpOrder, tv.EzCentering[2], tv.nbrCells[2])], dtype=np.float64)


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

    # analytical profiles of density, velocity...
    n  = np.tensordot(np.cosh(0.5*x_primal),
                      np.tensordot(np.cosh(0.5*y_primal),
                                   np.cosh(0.5*z_primal), axes=0), axes=0)
    Vx = np.tensordot(np.sinh(0.2*x_primal),
                      np.tensordot(np.sinh(0.2*y_primal),
                                   np.sinh(0.2*z_primal), axes=0), axes=0)
    Vy = np.tensordot(np.sinh(0.3*x_primal),
                      np.tensordot(np.sinh(0.3*y_primal),
                                   np.sinh(0.3*z_primal), axes=0), axes=0)
    Vz = np.tensordot(np.sinh(0.4*x_primal),
                      np.tensordot(np.sinh(0.4*y_primal),
                                   np.sinh(0.4*z_primal), axes=0), axes=0)
    P  = np.tensordot(np.cosh(0.5*x_primal),
                      np.tensordot(np.cosh(0.5*y_primal),
                                   np.cosh(0.5*z_primal), axes=0), axes=0)
    Bx = np.tensordot(np.cosh(0.2*x_primal),
                      np.tensordot(np.cosh(0.2*y_dual),
                                   np.cosh(0.2*z_dual), axes=0), axes=0)
    By = np.tensordot(np.cosh(0.3*x_dual),
                      np.tensordot(np.cosh(0.3*y_primal),
                                   np.cosh(0.3*z_dual), axes=0), axes=0)
    Bz = np.tensordot(np.cosh(0.4*x_dual),
                      np.tensordot(np.cosh(0.4*y_dual),
                                   np.cosh(0.4*z_primal), axes=0), axes=0)
    Jx = np.tensordot(np.tanh(0.2*x_dual),
                      np.tensordot(np.tanh(0.2*y_primal),
                                   np.tanh(0.2*z_primal), axes=0), axes=0)
    Jy = np.tensordot(np.tanh(0.3*x_primal),
                      np.tensordot(np.tanh(0.3*y_dual),
                                   np.tanh(0.3*z_primal), axes=0), axes=0)
    Jz = np.tensordot(np.tanh(0.4*x_primal),
                      np.tensordot(np.tanh(0.4*y_primal),
                                   np.tanh(0.4*z_dual), axes=0), axes=0)

    # ideal term
    idealx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1] = \
 -0.25*(Vy[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
       +Vy[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1])\
      *(Bz[psi_d_X  :pei_d_X+1, psi_d_Y-1:pei_d_Y+1, psi_p_Z  :pei_p_Z+1]\
       +Bz[psi_d_X  :pei_d_X+1, psi_d_Y  :pei_d_Y+2, psi_p_Z  :pei_p_Z+1])\
 +0.25*(Vz[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
       +Vz[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1])\
      *(By[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z-1:pei_d_Z+1]\
       +By[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+2])

    idealy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1] = \
 -0.25*(Vz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  , psi_p_Z  :pei_p_Z+1]\
       +Vz[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1, psi_p_Z  :pei_p_Z+1])\
      *(Bx[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_d_Z-1:pei_d_Z+1]\
       +Bx[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_d_Z  :pei_d_Z+2])\
 +0.25*(Vx[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  , psi_p_Z  :pei_p_Z+1]\
       +Vx[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1, psi_p_Z  :pei_p_Z+1])\
      *(Bz[psi_d_X-1:pei_d_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1]\
       +Bz[psi_d_X  :pei_d_X+2, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1])

    idealz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1] = \
 -0.25*(Vx[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z  ]\
       +Vx[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z+1:pei_p_Z+1])\
      *(By[psi_d_X-1:pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1]\
       +By[psi_d_X  :pei_d_X+2, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1])\
 +0.25*(Vy[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z  ]\
       +Vy[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z+1:pei_p_Z+1])\
      *(Bx[psi_p_X  :pei_p_X+1, psi_d_Y-1:pei_d_Y+1, psi_d_Z  :pei_d_Z+1]\
       +Bx[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+2, psi_d_Z  :pei_d_Z+1])

    # pressure term
    pressx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1] = \
       ((P[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
        -P[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1])/tv.meshSize[0])\
  /(0.5*(n[psi_p_X  :pei_p_X  , psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
        +n[psi_p_X+1:pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]))

    pressy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1] = \
       ((P[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
        -P[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  , psi_p_Z  :pei_p_Z+1])/tv.meshSize[1])\
  /(0.5*(n[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y  , psi_p_Z  :pei_p_Z+1]\
        +n[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+1, psi_p_Z  :pei_p_Z+1]))

    pressz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1] = \
       ((P[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z+1:pei_p_Z+1]\
        -P[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z  ])/tv.meshSize[2])\
  /(0.5*(n[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z  ]\
        +n[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z+1:pei_p_Z+1]))

    # resistive term
    resistx[psi_d_X:pei_d_X+1, psi_p_Y:pei_p_Y+1, psi_p_Z  :pei_p_Z+1] = \
     eta*Jx[psi_d_X:pei_d_X+1, psi_p_Y:pei_p_Y+1, psi_p_Z  :pei_p_Z+1]

    resisty[psi_p_X:pei_p_X+1, psi_d_Y:pei_d_Y+1, psi_p_Z  :pei_p_Z+1] = \
     eta*Jy[psi_p_X:pei_p_X+1, psi_d_Y:pei_d_Y+1, psi_p_Z  :pei_p_Z+1]

    resistz[psi_p_X:pei_p_X+1, psi_p_Y:pei_p_Y+1, psi_d_Z  :pei_d_Z+1] = \
     eta*Jz[psi_p_X:pei_p_X+1, psi_p_Y:pei_p_Y+1, psi_d_Z  :pei_d_Z+1]

    # viscous term
    viscousx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1] = \
     -nu*(Jx[psi_d_X-1:pei_d_X  , psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
     -2.0*Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
         +Jx[psi_d_X+1:pei_d_X+2, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1])\
        /(tv.meshSize[0]*tv.meshSize[0])\
     -nu*(Jx[psi_d_X  :pei_d_X+1, psi_p_Y-1:pei_p_Y  , psi_p_Z  :pei_p_Z+1]\
     -2.0*Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
         +Jx[psi_d_X  :pei_d_X+1, psi_p_Y+1:pei_p_Y+2, psi_p_Z  :pei_p_Z+1])\
        /(tv.meshSize[1]*tv.meshSize[1])\
     -nu*(Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z-1:pei_p_Z  ]\
     -2.0*Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z  :pei_p_Z+1]\
         +Jx[psi_d_X  :pei_d_X+1, psi_p_Y  :pei_p_Y+1, psi_p_Z+1:pei_p_Z+2])\
        /(tv.meshSize[2]*tv.meshSize[2])

    viscousy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1] = \
     -nu*(Jy[psi_p_X-1:pei_p_X  , psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1]\
     -2.0*Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1]\
         +Jy[psi_p_X+1:pei_p_X+2, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1])\
        /(tv.meshSize[0]*tv.meshSize[0])\
     -nu*(Jy[psi_p_X  :pei_p_X+1, psi_d_Y-1:pei_d_Y  , psi_p_Z  :pei_p_Z+1]\
     -2.0*Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1]\
         +Jy[psi_p_X  :pei_p_X+1, psi_d_Y+1:pei_d_Y+2, psi_p_Z  :pei_p_Z+1])\
        /(tv.meshSize[1]*tv.meshSize[1])\
     -nu*(Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z-1:pei_p_Z  ]\
     -2.0*Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z  :pei_p_Z+1]\
         +Jy[psi_p_X  :pei_p_X+1, psi_d_Y  :pei_d_Y+1, psi_p_Z+1:pei_p_Z+2])\
        /(tv.meshSize[2]*tv.meshSize[2])

    viscousz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1] = \
     -nu*(Jz[psi_p_X-1:pei_p_X  , psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1]\
     -2.0*Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1]\
         +Jz[psi_p_X+1:pei_p_X+2, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1])\
        /(tv.meshSize[0]*tv.meshSize[0])\
     -nu*(Jz[psi_p_X  :pei_p_X+1, psi_p_Y-1:pei_p_Y  , psi_d_Z  :pei_d_Z+1]\
     -2.0*Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1]\
         +Jz[psi_p_X  :pei_p_X+1, psi_p_Y+1:pei_p_Y+2, psi_d_Z  :pei_d_Z+1])\
        /(tv.meshSize[1]*tv.meshSize[1])\
     -nu*(Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z-1:pei_d_Z  ]\
     -2.0*Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z  :pei_d_Z+1]\
         +Jz[psi_p_X  :pei_p_X+1, psi_p_Y  :pei_p_Y+1, psi_d_Z+1:pei_d_Z+2])\
        /(tv.meshSize[2]*tv.meshSize[2])

    ExNew = idealx+pressx+resistx+viscousx
    EyNew = idealy+pressy+resisty+viscousy
    EzNew = idealz+pressz+resistz+viscousz

    filename_ohmx = "ohmx_yee_3D_order1.txt"
    filename_ohmy = "ohmy_yee_3D_order1.txt"
    filename_ohmz = "ohmz_yee_3D_order1.txt"

    np.savetxt(os.path.join(path, filename_ohmx), ExNew.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_ohmy), EyNew.flatten('C'), delimiter=" ")
    np.savetxt(os.path.join(path, filename_ohmz), EzNew.flatten('C'), delimiter=" ")



def main(path='./'):

    if len(sys.argv) >1:
        path=sys.argv[1]

    test_ohm_yee1D(path)
    test_ohm_yee2D(path)
    test_ohm_yee3D(path)


if __name__ == '__main__':
    main()

