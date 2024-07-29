# the goal of this script is to generate a file for the laplacian unit test
# the file will contain expected values for the unit testself.


import os
import sys
import numpy as np
from pyphare.core import gridlayout


class TestVariables(object):
    def __init__(self):
        """This class only has a constructor
        in order to set the variables to test the derivatives
        all quantities are 3D

        """

        self.nbrCells = (50, 30, 40)
        self.meshSize = (0.1, 0.2, 0.3)
        self.interpOrders = {1, 2, 3}
        self.JxCentering = ("dual", "primal", "primal")
        self.JyCentering = ("primal", "dual", "primal")
        self.JzCentering = ("primal", "primal", "dual")
        self.domainSize = tuple(n * m for n, m in zip(self.nbrCells, self.meshSize))


# ------------------------------------------------------------------------------
# since the derivative is tested for all interporders
# we don't do ohm for all interporders, only for interporder=1
# ------------------------------------------------------------------------------


def test_laplacian_yee1D(path):
    layout = gridlayout.GridLayout()  # yee layout

    tv = TestVariables()

    for interpOrder in tv.interpOrders:
        Jx = np.zeros(
            layout.allocSize(interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
            dtype=np.float64,
        )
        Jy = np.zeros(
            layout.allocSize(interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
            dtype=np.float64,
        )
        Jz = np.zeros(
            layout.allocSize(interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
            dtype=np.float64,
        )

        lapJx = np.zeros(
            layout.allocSize(interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
            dtype=np.float64,
        )
        lapJy = np.zeros(
            layout.allocSize(interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
            dtype=np.float64,
        )
        lapJz = np.zeros(
            layout.allocSize(interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
            dtype=np.float64,
        )

        psi_p_X = layout.physicalStartIndex(interpOrder, "primal")
        pei_p_X = layout.physicalEndIndex(interpOrder, "primal", tv.nbrCells[0])

        psi_d_X = layout.physicalStartIndex(interpOrder, "dual")
        pei_d_X = layout.physicalEndIndex(interpOrder, "dual", tv.nbrCells[0])

        nbrGhost_p = layout.nbrGhosts(interpOrder, "primal")
        nbrGhost_d = layout.nbrGhosts(interpOrder, "dual")

        x_primal = (
            tv.meshSize[0]
            * np.arange(layout.allocSize(interpOrder, "primal", tv.nbrCells[0]))
            - tv.meshSize[0] * nbrGhost_p
        )
        x_dual = (
            tv.meshSize[0]
            * np.arange(layout.allocSize(interpOrder, "dual", tv.nbrCells[0]))
            - tv.meshSize[0] * nbrGhost_d
            + tv.meshSize[0] * 0.5
        )

        Jx = np.sinh(0.1 * x_dual)
        Jy = np.sinh(0.3 * x_primal)
        Jz = np.sinh(0.2 * x_primal)

        lapJx[psi_d_X : pei_d_X + 1] = (
            Jx[psi_d_X + 1 : pei_d_X + 2]
            - 2.0 * Jx[psi_d_X : pei_d_X + 1]
            + Jx[psi_d_X - 1 : pei_d_X]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        lapJy[psi_p_X : pei_p_X + 1] = (
            Jy[psi_p_X + 1 : pei_p_X + 2]
            - 2.0 * Jy[psi_p_X : pei_p_X + 1]
            + Jy[psi_p_X - 1 : pei_p_X]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        lapJz[psi_p_X : pei_p_X + 1] = (
            Jz[psi_p_X + 1 : pei_p_X + 2]
            - 2.0 * Jz[psi_p_X : pei_p_X + 1]
            + Jz[psi_p_X - 1 : pei_p_X]
        ) / (tv.meshSize[0] * tv.meshSize[0])

        filename_lapJx = "lapJx_interpOrder_{}_1d.txt".format(interpOrder)
        filename_lapJy = "lapJy_interpOrder_{}_1d.txt".format(interpOrder)
        filename_lapJz = "lapJz_interpOrder_{}_1d.txt".format(interpOrder)

        np.savetxt(os.path.join(path, filename_lapJx), lapJx, delimiter=" ")
        np.savetxt(os.path.join(path, filename_lapJy), lapJy, delimiter=" ")
        np.savetxt(os.path.join(path, filename_lapJz), lapJz, delimiter=" ")


def test_laplacian_yee2D(path):
    layout = gridlayout.GridLayout()  # yee layout

    tv = TestVariables()

    for interpOrder in tv.interpOrders:
        Jx = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JxCentering[1], tv.nbrCells[1]),
            ],
            dtype=np.float64,
        )
        Jy = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JyCentering[1], tv.nbrCells[1]),
            ],
            dtype=np.float64,
        )
        Jz = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JzCentering[1], tv.nbrCells[1]),
            ],
            dtype=np.float64,
        )
        Jx_x = np.zeros_like(Jx)
        Jx_y = np.zeros_like(Jx)
        Jy_x = np.zeros_like(Jy)
        Jy_y = np.zeros_like(Jy)
        Jz_x = np.zeros_like(Jz)
        Jz_y = np.zeros_like(Jz)

        lapJx = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JxCentering[1], tv.nbrCells[1]),
            ],
            dtype=np.float64,
        )
        lapJy = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JyCentering[1], tv.nbrCells[1]),
            ],
            dtype=np.float64,
        )
        lapJz = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JzCentering[1], tv.nbrCells[1]),
            ],
            dtype=np.float64,
        )

        psi_p_X = layout.physicalStartIndex(interpOrder, "primal")
        pei_p_X = layout.physicalEndIndex(interpOrder, "primal", tv.nbrCells[0])
        psi_p_Y = layout.physicalStartIndex(interpOrder, "primal")
        pei_p_Y = layout.physicalEndIndex(interpOrder, "primal", tv.nbrCells[1])

        psi_d_X = layout.physicalStartIndex(interpOrder, "dual")
        pei_d_X = layout.physicalEndIndex(interpOrder, "dual", tv.nbrCells[0])
        psi_d_Y = layout.physicalStartIndex(interpOrder, "dual")
        pei_d_Y = layout.physicalEndIndex(interpOrder, "dual", tv.nbrCells[1])

        nbrGhost_p = layout.nbrGhosts(interpOrder, "primal")
        nbrGhost_d = layout.nbrGhosts(interpOrder, "dual")

        x_primal = (
            tv.meshSize[0]
            * np.arange(layout.allocSize(interpOrder, "primal", tv.nbrCells[0]))
            - tv.meshSize[0] * nbrGhost_p
        )
        x_dual = (
            tv.meshSize[0]
            * np.arange(layout.allocSize(interpOrder, "dual", tv.nbrCells[0]))
            - tv.meshSize[0] * nbrGhost_d
            + tv.meshSize[0] * 0.5
        )
        y_primal = (
            tv.meshSize[1]
            * np.arange(layout.allocSize(interpOrder, "primal", tv.nbrCells[1]))
            - tv.meshSize[1] * nbrGhost_p
        )
        y_dual = (
            tv.meshSize[1]
            * np.arange(layout.allocSize(interpOrder, "dual", tv.nbrCells[1]))
            - tv.meshSize[1] * nbrGhost_d
            + tv.meshSize[1] * 0.5
        )

        Jx = np.tensordot(np.sinh(0.1 * x_dual), np.cosh(0.1 * y_primal), axes=0)
        Jy = np.tensordot(np.sinh(0.3 * x_primal), np.cosh(0.3 * y_dual), axes=0)
        Jz = np.tensordot(np.sinh(0.2 * x_primal), np.cosh(0.2 * y_primal), axes=0)

        Jx_x[psi_d_X : pei_d_X + 1, :] = (
            Jx[psi_d_X + 1 : pei_d_X + 2, :]
            - 2.0 * Jx[psi_d_X : pei_d_X + 1, :]
            + Jx[psi_d_X - 1 : pei_d_X, :]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        Jx_y[:, psi_p_Y : pei_p_Y + 1] = (
            Jx[:, psi_p_Y + 1 : pei_p_Y + 2]
            - 2.0 * Jx[:, psi_p_Y : pei_p_Y + 1]
            + Jx[:, psi_p_Y - 1 : pei_p_Y]
        ) / (tv.meshSize[1] * tv.meshSize[1])

        Jy_x[psi_p_X : pei_p_X + 1, :] = (
            Jy[psi_p_X + 1 : pei_p_X + 2, :]
            - 2.0 * Jy[psi_p_X : pei_p_X + 1, :]
            + Jy[psi_p_X - 1 : pei_p_X, :]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        Jy_y[:, psi_d_Y : pei_d_Y + 1] = (
            Jy[:, psi_d_Y + 1 : pei_d_Y + 2]
            - 2.0 * Jy[:, psi_d_Y : pei_d_Y + 1]
            + Jy[:, psi_d_Y - 1 : pei_d_Y]
        ) / (tv.meshSize[1] * tv.meshSize[1])

        Jz_x[psi_p_X : pei_p_X + 1, :] = (
            Jz[psi_p_X + 1 : pei_p_X + 2, :]
            - 2.0 * Jz[psi_p_X : pei_p_X + 1, :]
            + Jz[psi_p_X - 1 : pei_p_X, :]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        Jz_y[:, psi_p_Y : pei_p_Y + 1] = (
            Jz[:, psi_p_Y + 1 : pei_p_Y + 2]
            - 2.0 * Jz[:, psi_p_Y : pei_p_Y + 1]
            + Jz[:, psi_p_Y - 1 : pei_p_Y]
        ) / (tv.meshSize[1] * tv.meshSize[1])

        lapJx = Jx_x + Jx_y
        lapJy = Jy_x + Jy_y
        lapJz = Jz_x + Jz_y

        filename_lapJx = "lapJx_interpOrder_{}_2d.txt".format(interpOrder)
        filename_lapJy = "lapJy_interpOrder_{}_2d.txt".format(interpOrder)
        filename_lapJz = "lapJz_interpOrder_{}_2d.txt".format(interpOrder)

        np.savetxt(os.path.join(path, filename_lapJx), lapJx, delimiter=" ")
        np.savetxt(os.path.join(path, filename_lapJy), lapJy, delimiter=" ")
        np.savetxt(os.path.join(path, filename_lapJz), lapJz, delimiter=" ")


def test_laplacian_yee3D(path):
    layout = gridlayout.GridLayout()  # yee layout

    tv = TestVariables()

    for interpOrder in tv.interpOrders:
        Jx = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JxCentering[1], tv.nbrCells[1]),
                layout.allocSize(interpOrder, tv.JxCentering[2], tv.nbrCells[2]),
            ],
            dtype=np.float64,
        )
        Jy = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JyCentering[1], tv.nbrCells[1]),
                layout.allocSize(interpOrder, tv.JyCentering[2], tv.nbrCells[2]),
            ],
            dtype=np.float64,
        )
        Jz = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JzCentering[1], tv.nbrCells[1]),
                layout.allocSize(interpOrder, tv.JzCentering[2], tv.nbrCells[2]),
            ],
            dtype=np.float64,
        )
        Jx_x = np.zeros_like(Jx)
        Jx_y = np.zeros_like(Jx)
        Jx_z = np.zeros_like(Jx)
        Jy_x = np.zeros_like(Jy)
        Jy_y = np.zeros_like(Jy)
        Jy_z = np.zeros_like(Jy)
        Jz_x = np.zeros_like(Jz)
        Jz_y = np.zeros_like(Jz)
        Jz_z = np.zeros_like(Jz)

        lapJx = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JxCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JxCentering[1], tv.nbrCells[1]),
                layout.allocSize(interpOrder, tv.JxCentering[2], tv.nbrCells[2]),
            ],
            dtype=np.float64,
        )
        lapJy = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JyCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JyCentering[1], tv.nbrCells[1]),
                layout.allocSize(interpOrder, tv.JyCentering[2], tv.nbrCells[2]),
            ],
            dtype=np.float64,
        )
        lapJz = np.zeros(
            [
                layout.allocSize(interpOrder, tv.JzCentering[0], tv.nbrCells[0]),
                layout.allocSize(interpOrder, tv.JzCentering[1], tv.nbrCells[1]),
                layout.allocSize(interpOrder, tv.JzCentering[2], tv.nbrCells[2]),
            ],
            dtype=np.float64,
        )

        psi_p_X = layout.physicalStartIndex(interpOrder, "primal")
        pei_p_X = layout.physicalEndIndex(interpOrder, "primal", tv.nbrCells[0])
        psi_p_Y = layout.physicalStartIndex(interpOrder, "primal")
        pei_p_Y = layout.physicalEndIndex(interpOrder, "primal", tv.nbrCells[1])
        psi_p_Z = layout.physicalStartIndex(interpOrder, "primal")
        pei_p_Z = layout.physicalEndIndex(interpOrder, "primal", tv.nbrCells[2])

        psi_d_X = layout.physicalStartIndex(interpOrder, "dual")
        pei_d_X = layout.physicalEndIndex(interpOrder, "dual", tv.nbrCells[0])
        psi_d_Y = layout.physicalStartIndex(interpOrder, "dual")
        pei_d_Y = layout.physicalEndIndex(interpOrder, "dual", tv.nbrCells[1])
        psi_d_Z = layout.physicalStartIndex(interpOrder, "dual")
        pei_d_Z = layout.physicalEndIndex(interpOrder, "dual", tv.nbrCells[2])

        nbrGhost_p = layout.nbrGhosts(interpOrder, "primal")
        nbrGhost_d = layout.nbrGhosts(interpOrder, "dual")

        x_primal = (
            tv.meshSize[0]
            * np.arange(layout.allocSize(interpOrder, "primal", tv.nbrCells[0]))
            - tv.meshSize[0] * nbrGhost_p
        )
        x_dual = (
            tv.meshSize[0]
            * np.arange(layout.allocSize(interpOrder, "dual", tv.nbrCells[0]))
            - tv.meshSize[0] * nbrGhost_d
            + tv.meshSize[0] * 0.5
        )
        y_primal = (
            tv.meshSize[1]
            * np.arange(layout.allocSize(interpOrder, "primal", tv.nbrCells[1]))
            - tv.meshSize[1] * nbrGhost_p
        )
        y_dual = (
            tv.meshSize[1]
            * np.arange(layout.allocSize(interpOrder, "dual", tv.nbrCells[1]))
            - tv.meshSize[1] * nbrGhost_d
            + tv.meshSize[1] * 0.5
        )
        z_primal = (
            tv.meshSize[2]
            * np.arange(layout.allocSize(interpOrder, "primal", tv.nbrCells[2]))
            - tv.meshSize[2] * nbrGhost_p
        )
        z_dual = (
            tv.meshSize[2]
            * np.arange(layout.allocSize(interpOrder, "dual", tv.nbrCells[2]))
            - tv.meshSize[2] * nbrGhost_d
            + tv.meshSize[2] * 0.5
        )

        Jx = np.tensordot(
            np.sinh(0.1 * x_dual),
            np.tensordot(np.cosh(0.1 * y_primal), np.tanh(0.1 * z_primal), axes=0),
            axes=0,
        )
        Jy = np.tensordot(
            np.sinh(0.3 * x_primal),
            np.tensordot(np.cosh(0.3 * y_dual), np.tanh(0.3 * z_primal), axes=0),
            axes=0,
        )
        Jz = np.tensordot(
            np.sinh(0.2 * x_primal),
            np.tensordot(np.cosh(0.2 * y_primal), np.tanh(0.2 * z_dual), axes=0),
            axes=0,
        )

        Jx_x[psi_d_X : pei_d_X + 1, :, :] = (
            Jx[psi_d_X + 1 : pei_d_X + 2, :, :]
            - 2.0 * Jx[psi_d_X : pei_d_X + 1, :, :]
            + Jx[psi_d_X - 1 : pei_d_X, :, :]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        Jx_y[:, psi_p_Y : pei_p_Y + 1, :] = (
            Jx[:, psi_p_Y + 1 : pei_p_Y + 2, :]
            - 2.0 * Jx[:, psi_p_Y : pei_p_Y + 1, :]
            + Jx[:, psi_p_Y - 1 : pei_p_Y, :]
        ) / (tv.meshSize[1] * tv.meshSize[1])
        Jx_z[:, :, psi_p_Z : pei_p_Z + 1] = (
            Jx[:, :, psi_p_Z + 1 : pei_p_Z + 2]
            - 2.0 * Jx[:, :, psi_p_Z : pei_p_Z + 1]
            + Jx[:, :, psi_p_Z - 1 : pei_p_Z]
        ) / (tv.meshSize[2] * tv.meshSize[2])

        Jy_x[psi_p_X : pei_p_X + 1, :, :] = (
            Jy[psi_p_X + 1 : pei_p_X + 2, :, :]
            - 2.0 * Jy[psi_p_X : pei_p_X + 1, :, :]
            + Jy[psi_p_X - 1 : pei_p_X, :, :]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        Jy_y[:, psi_d_Y : pei_d_Y + 1, :] = (
            Jy[:, psi_d_Y + 1 : pei_d_Y + 2, :]
            - 2.0 * Jy[:, psi_d_Y : pei_d_Y + 1, :]
            + Jy[:, psi_d_Y - 1 : pei_d_Y, :]
        ) / (tv.meshSize[1] * tv.meshSize[1])
        Jy_z[:, :, psi_p_Z : pei_p_Z + 1] = (
            Jy[:, :, psi_p_Z + 1 : pei_p_Z + 2]
            - 2.0 * Jy[:, :, psi_p_Z : pei_p_Z + 1]
            + Jy[:, :, psi_p_Z - 1 : pei_p_Z]
        ) / (tv.meshSize[2] * tv.meshSize[2])

        Jz_x[psi_p_X : pei_p_X + 1, :, :] = (
            Jz[psi_p_X + 1 : pei_p_X + 2, :, :]
            - 2.0 * Jz[psi_p_X : pei_p_X + 1, :, :]
            + Jz[psi_p_X - 1 : pei_p_X, :, :]
        ) / (tv.meshSize[0] * tv.meshSize[0])
        Jz_y[:, psi_p_Y : pei_p_Y + 1, :] = (
            Jz[:, psi_p_Y + 1 : pei_p_Y + 2, :]
            - 2.0 * Jz[:, psi_p_Y : pei_p_Y + 1, :]
            + Jz[:, psi_p_Y - 1 : pei_p_Y, :]
        ) / (tv.meshSize[1] * tv.meshSize[1])
        Jz_z[:, :, psi_d_Z : pei_d_Z + 1] = (
            Jz[:, :, psi_d_Z + 1 : pei_d_Z + 2]
            - 2.0 * Jz[:, :, psi_d_Z : pei_d_Z + 1]
            + Jz[:, :, psi_d_Z - 1 : pei_d_Z]
        ) / (tv.meshSize[2] * tv.meshSize[2])

        lapJx = Jx_x + Jx_y + Jx_z
        lapJy = Jy_x + Jy_y + Jy_z
        lapJz = Jz_x + Jz_y + Jz_z

        filename_lapJx = "lapJx_interpOrder_{}_3d.txt".format(interpOrder)
        filename_lapJy = "lapJy_interpOrder_{}_3d.txt".format(interpOrder)
        filename_lapJz = "lapJz_interpOrder_{}_3d.txt".format(interpOrder)

        np.savetxt(
            os.path.join(path, filename_lapJx), lapJx.flatten("C"), delimiter=" "
        )
        np.savetxt(
            os.path.join(path, filename_lapJy), lapJy.flatten("C"), delimiter=" "
        )
        np.savetxt(
            os.path.join(path, filename_lapJz), lapJz.flatten("C"), delimiter=" "
        )


def main(path="./"):
    if len(sys.argv) > 1:
        path = sys.argv[1]

    test_laplacian_yee1D(path)
    test_laplacian_yee2D(path)
    test_laplacian_yee3D(path)


if __name__ == "__main__":
    main()
