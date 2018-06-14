
# the goal of this scrip is to generate a file for the faraday unit test
# the file will contain expected values for the unit testself.

# the code consists in calculating the curl of a function numerically 


import numpy as np
import sys
import os

def test_faraday_yee1D(path):

    # since the derivative is tested for all interporders
    # we don't do faraday for all interporders, only for interporder=1

    ghost = 5

    nbrCellsX = 50
    dualAllocSize = nbrCellsX + 2*ghost
    primalAllocSize = nbrCellsX + 1 + 2*ghost
    startIndexX = 0
    endIndexX = startIndexX + 50
    meshSize=0.1
    xmax = 5.

    dBydt = np.zeros(dualAllocSize, dtype=np.float64)
    dBzdt = np.zeros(dualAllocSize, dtype=np.float64)

    psi_p = ghost
    pei_p = nbrCellsX + ghost
    psi_d = ghost
    pei_d = nbrCellsX + ghost -1

    x_dual = meshSize*np.arange(startIndexX - ghost, endIndexX + ghost) + meshSize/2.
    x_primal = meshSize*np.arange(startIndexX - ghost, endIndexX + ghost + 1)

    Ey = np.cos(2*np.pi/xmax * x_primal)
    Ez = np.sin(2*np.pi/xmax * x_primal)
    By = np.tanh(x_dual-  xmax/2.)
    Bz = np.tanh(x_dual-  xmax/2.)

    # dBxdt = -(dyEz-dzEy) = 0
    # dBydt = -(dzEx - dxEz) = dxEz
    # dBzdt = -(dxEy - dyEx) = -dxEy
    dBydt[psi_d:pei_d+1] = By[psi_d:pei_d+1] +  ( Ez[psi_p+1:pei_p+1]  - Ez[psi_p:pei_p])/meshSize
    dBzdt[psi_d:pei_d+1] = Bz[psi_d:pei_d+1] -( Ey[psi_p+1:pei_p+1]  - Ey[psi_p:pei_p])/meshSize


    filename_dbydt = "dbydt_yee_1D_order1.txt"
    filename_dbzdt = "dbzdt_yee_1D_order1.txt"

    np.savetxt(os.path.join(path, filename_dbydt), dBydt, delimiter=" ")
    np.savetxt(os.path.join(path, filename_dbzdt), dBzdt, delimiter=" ")


def main(path='./'):

    if len(sys.argv) >1:
        path=sys.argv[1]

    test_faraday_yee1D(path)




if __name__ == '__main__':
    main()
