
# the goal of this script is to generate a file for the ampere unit test
# the file will contain expected values for the unit test.

# the code consists in calculating the curl of a function numerically 


import numpy as np
import sys
import os

def test_ampere_yee1D(path):

    # since the derivative is tested for all interporders
    # we don't do ampere for all interporders, only for interporder=1

    ghost = 5

    nbrCellsX = 50
    dualAllocSize = nbrCellsX + 2* ghost
    primalAllocSize = nbrCellsX + 1 + 2* ghost
    startIndexX = 0
    endIndexX = startIndexX + nbrCellsX
    meshSize=0.1
    xmax = 5.

    Jx = np.zeros(dualAllocSize, dtype=np.float64)
    Jy = np.zeros(primalAllocSize, dtype=np.float64)
    Jz = np.zeros(primalAllocSize, dtype=np.float64)

    psi_p = ghost
    pei_p = nbrCellsX + ghost
    psi_d = ghost
    pei_d = nbrCellsX + ghost -1

    x_dual = meshSize*np.arange(startIndexX - ghost, endIndexX + ghost) + meshSize/2.
    x_primal = meshSize*np.arange(startIndexX -ghost, endIndexX + ghost + 1)

    By = np.cos(2*np.pi/xmax * x_dual)
    Bz = np.sin(2*np.pi/xmax * x_dual)

    # Jy = dzBx - dxBz  -dxBz
    # Jz = dxBy - dyBx   dxBy
    Jy[psi_p:pei_p+1] = -( Bz[psi_d:pei_d+2]  - Bz[psi_d-1:pei_d+1])/meshSize
    Jz[psi_p:pei_p+1] =  ( By[psi_d:pei_d+2]  - By[psi_d-1:pei_d+1])/meshSize


    filename_jy = "jy_yee_1D_order1.txt"
    filename_jz = "jz_yee_1D_order1.txt"

    np.savetxt(os.path.join(path, filename_jy), Jy, delimiter=" ")
    np.savetxt(os.path.join(path, filename_jz), Jz, delimiter=" ")


def main(path='./'):

    if len(sys.argv) >1:
        path=sys.argv[1]

    test_ampere_yee1D(path)




if __name__ == '__main__':
    main()
