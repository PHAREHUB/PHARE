#!/usr/bin/env python3


import numpy as np

from .box import Box
from .phare_utilities import is_scalar, listify

directions = ["x", "y", "z"]
direction_to_dim = {direction: idx for idx, direction in enumerate(directions)}

yee_centering = {
    "x": {
        "Bx": "primal",
        "By": "dual",
        "Bz": "dual",
        "Ex": "dual",
        "Ey": "primal",
        "Ez": "primal",
        "Jx": "dual",
        "Jy": "primal",
        "Jz": "primal",
        "Vx": "primal",
        "Vy": "primal",
        "Vz": "primal",
        "Fx": "primal",
        "Fy": "primal",
        "Fz": "primal",
        "rho": "primal",
        "P": "primal",
        "Vthx": "primal",
        "Vthy": "primal",
        "Vthz": "primal",
        "Mxx": "primal",
        "Mxy": "primal",
        "Mxz": "primal",
        "Myy": "primal",
        "Myz": "primal",
        "Mzz": "primal",
        "Pxx": "primal",
        "Pxy": "primal",
        "Pxz": "primal",
        "Pyy": "primal",
        "Pyz": "primal",
        "Pzz": "primal",
        "tags": "dual",
    },
    "y": {
        "Bx": "dual",
        "By": "primal",
        "Bz": "dual",
        "Ex": "primal",
        "Ey": "dual",
        "Ez": "primal",
        "Jx": "primal",
        "Jy": "dual",
        "Jz": "primal",
        "Vx": "primal",
        "Vy": "primal",
        "Vz": "primal",
        "Fx": "primal",
        "Fy": "primal",
        "Fz": "primal",
        "rho": "primal",
        "P": "primal",
        "Vthx": "primal",
        "Vthy": "primal",
        "Vthz": "primal",
        "Mxx": "primal",
        "Mxy": "primal",
        "Mxz": "primal",
        "Myy": "primal",
        "Myz": "primal",
        "Mzz": "primal",
        "Pxx": "primal",
        "Pxy": "primal",
        "Pxz": "primal",
        "Pyy": "primal",
        "Pyz": "primal",
        "Pzz": "primal",
        "tags": "dual",
    },
    "z": {
        "Bx": "dual",
        "By": "dual",
        "Bz": "primal",
        "Ex": "primal",
        "Ey": "primal",
        "Ez": "dual",
        "Jx": "primal",
        "Jy": "primal",
        "Jz": "dual",
        "Vx": "primal",
        "Vy": "primal",
        "Vz": "primal",
        "Fx": "primal",
        "Fy": "primal",
        "Fz": "primal",
        "rho": "primal",
        "P": "primal",
        "Vthx": "primal",
        "Vthy": "primal",
        "Vthz": "primal",
        "Mxx": "primal",
        "Mxy": "primal",
        "Mxz": "primal",
        "Myy": "primal",
        "Myz": "primal",
        "Mzz": "primal",
        "Pxx": "primal",
        "Pxy": "primal",
        "Pxz": "primal",
        "Pyy": "primal",
        "Pyz": "primal",
        "Pzz": "primal",
        "tags": "dual",
    },
}
yee_centering_lower = {
    direction: {key.lower(): centering for key, centering in key_dict.items()}
    for direction, key_dict in yee_centering.items()
}


def yee_element_is_primal(key, direction="x"):
    if key in yee_centering_lower[direction]:
        return yee_centering_lower[direction][key] == "primal"
    return yee_centering[direction][key] == "primal"


class YeeCentering(object):
    def __init__(self):
        self.centerX = yee_centering["x"]
        self.centerY = yee_centering["y"]
        self.centerZ = yee_centering["z"]


def yeeCoordsFor(
    origin, nbrGhosts, dl, nbrCells, qty, direction, withGhosts=False, **kwargs
):
    assert direction in direction_to_dim, f"direction ({direction} not supported)"
    if qty in yee_centering_lower[direction] and qty not in yee_centering[direction]:
        qty = qty[0].upper() + qty[1:]

    if "centering" in kwargs:
        centering = kwargs["centering"]
    else:
        centering = yee_centering[direction][qty]

    offset = 0
    dim = direction_to_dim[direction]
    if withGhosts:
        size = nbrCells[dim] + (nbrGhosts * 2)
    else:
        size = nbrCells[dim]

    if centering == "dual":
        offset = 0.5 * dl[dim]
    else:
        size += 1

    if withGhosts:
        return origin[dim] - nbrGhosts * dl[dim] + np.arange(size) * dl[dim] + offset
    else:
        return origin[dim] + np.arange(size) * dl[dim] + offset


class GridLayout(object):
    """
    field_ghosts_nbr is a parameter to support pyphare geometry tests having hard coded 5 ghosts
    initialized default to -1 as an invalid value allowing the override mechanism. Using None
    results in a pylint error elsewhere
    """

    def __init__(
        self, box=Box(0, 0), origin=0, dl=0.1, interp_order=1, field_ghosts_nbr=-1
    ):
        self.box = box

        self.dl = listify(dl)
        assert len(self.dl) == box.ndim

        self.origin = listify(origin)
        assert len(self.origin) == box.ndim

        self.interp_order = interp_order
        self.impl = "yee"
        self.directions = ["X", "Y", "Z"]

        self.hybridQuantities = [
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
            "Fx",
            "Fy",
            "Fz",
        ]

        self.yeeCentering = YeeCentering()

        self.centering = {
            "X": self.yeeCentering.centerX,
            "Y": self.yeeCentering.centerY,
            "Z": self.yeeCentering.centerZ,
        }

        self.field_ghosts_nbr = field_ghosts_nbr  # allows override

    @property
    def ndim(self):
        return self.box.ndim

    def qtyCentering(self, quantity, direction):
        return self.centering[direction][quantity]

    def particleGhostNbr(self, interp_order):
        return 1 if interp_order == 1 else 2

    def nbrGhosts(self, interpOrder, centering):
        if self.field_ghosts_nbr == -1:
            return int((interpOrder + 1) / 2) + self.particleGhostNbr(interpOrder)
        return self.field_ghosts_nbr

    def nbrGhostsPrimal(self, interpOrder):
        return self.nbrGhosts(interpOrder, "primal")

    def qtyIsDual(self, qty, direction):
        return self.isDual(self.qtyCentering(qty, direction))

    def isDual(self, centering):
        if centering == "dual":
            return 1
        else:
            return 0

    def ghostStartIndex(self):
        return 0

    def ghostEndIndex(self, interpOrder, centering, nbrCells):
        index = self.physicalEndIndex(
            interpOrder, centering, nbrCells
        ) + self.nbrGhosts(interpOrder, centering)
        return index

    def physicalStartIndex(self, interpOrder, centering):
        index = self.ghostStartIndex() + self.nbrGhosts(interpOrder, centering)
        return index

    def physicalEndIndex(self, interpOrder, centering, nbrCells):
        index = (
            self.physicalStartIndex(interpOrder, centering)
            + nbrCells
            - self.isDual(centering)
        )
        return index

    def physicalStartIndices(self, qty):
        assert qty in self.hybridQuantities
        indices = np.zeros(self.box.ndim)
        for i, direction in enumerate(directions[: self.box.ndim]):
            centering = yee_centering[direction][qty]
            indices[i] = self.physicalStartIndex(self.interp_order, centering)
        return indices

    def physicalEndIndices(self, qty):
        assert qty in self.hybridQuantities
        indices = np.zeros(self.box.ndim)
        for i, direction in enumerate(directions[: self.box.ndim]):
            centering = yee_centering[direction][qty]
            indices[i] = self.physicalEndIndex(
                self.interp_order, centering, self.box.shape[i]
            )
        return indices

    def nbrGhostFor(self, qty):
        assert qty in self.hybridQuantities
        nGhosts = np.zeros(self.box.ndim, dtype=np.int32)
        for i, direction in enumerate(directions[: self.box.ndim]):
            centering = yee_centering[direction][qty]
            nGhosts[i] = self.nbrGhosts(self.interp_order, centering)
        return nGhosts

    # ---- Start / End   primal methods ------
    def physicalStartPrimal(self, interpOrder):
        index = self.ghostStartIndex() + self.nbrGhostsPrimal(interpOrder)
        return index

    def physicalEndPrimal(self, interpOrder, nbrCells):
        index = self.physicalStartPrimal(interpOrder) + nbrCells
        return index

    # ---- Alloc methods -------------------------

    def allocSize(self, interpOrder, centering, nbrCells):
        size = (
            nbrCells
            + 1
            + 2 * self.nbrGhosts(interpOrder, centering)
            - self.isDual(centering)
        )
        return size

    # 1st derivative
    def allocSizeDerived(self, interpOrder, centering, nbrCells):
        newCentering = self.changeCentering(centering, 1)

        size = (
            nbrCells
            + 1
            + 2 * self.nbrGhosts(interpOrder, newCentering)
            - self.isDual(newCentering)
        )
        return size

    def AMRPointToLocal(self, point):
        return (
            point - self.box.lower + self.physicalStartIndex(self.interp_order, "dual")
        )

    def AMRBoxToLocal(self, box):
        local = box.copy()
        local.lower = self.AMRPointToLocal(local.lower)
        local.upper = self.AMRPointToLocal(local.upper)
        return local

    def AMRToLocal(self, toLocal):
        assert not is_scalar(toLocal)
        if isinstance(toLocal, Box):
            return self.AMRBoxToLocal(toLocal)
        return self.AMRPointToLocal(toLocal)

    def AMRIndexToLocal(self, dim, index):
        assert is_scalar(index)
        dualStart = self.physicalStartIndex(self.interp_order, "dual")
        return index - self.box.lower[dim] + dualStart

    # ---- Yee coordinate methods -------------------------
    # knode : a primal or dual node index
    #
    # The centering deduced from qty and direction tells
    # whether knode is primal or dual
    #
    # ds stands for dx or dy or dz
    # This method returns a point
    #
    def yeeCoords(self, knode, iStart, centering, direction, ds, origin, derivOrder):
        halfCell = 0.0

        newCentering = self.changeCentering(centering, derivOrder)

        if newCentering == "dual":
            halfCell = 0.5

        x = ((knode - iStart) + halfCell) * ds + origin

        return x

    def yeeCoordsFor(self, qty, direction, withGhosts=True, **kwargs):
        """
        from a qty and a direction, returns a 1d array containing
        the coordinates where the qty is defined, including the ghost nodes

        :param qty: the quantity (can be primal or dual)
        :param direction: can only be a single one
        """

        if (
            qty in yee_centering_lower[direction]
            and qty not in yee_centering[direction]
        ):
            qty = qty[0].upper() + qty[1:]

        if "centering" in kwargs:
            centering = kwargs["centering"]
        else:
            centering = yee_centering[direction][qty]

        return yeeCoordsFor(
            self.origin,
            self.nbrGhosts(self.interp_order, centering),
            self.dl,
            self.box.shape,
            qty,
            direction,
            withGhosts=withGhosts,
            centering=centering,
        )

    # ---- Get coordinate methods -------------------------
    # knode : a primal or dual node index
    #
    # The centering deduced from qty and direction tells
    # whether knode is primal or dual
    #
    # ds stands for dx or dy or dz
    # This method returns a point
    #
    def fieldCoords(self, knode, iStart, qty, direction, ds, origin, derivOrder):
        halfCell = 0.0

        newCentering = self.changeCentering(
            self.qtyCentering(qty, direction), derivOrder
        )

        if newCentering == "dual":
            halfCell = 0.5

        x = ((knode - iStart) + halfCell) * ds + origin

        return x

    # ---- Change centering method -------------------------
    #
    # Use case:
    #   changeCentering( qtyCentering(qty, direct), 1 )
    #
    def changeCentering(self, centering, derivOrder=0):
        newCentering = centering

        # if derivOrder is odd the centering is changed
        if derivOrder % 2 != 0:
            newCentering = self.swapCentering(centering)

        return newCentering

    # -------------------------------------------------------
    def swapCentering(self, centering):
        newCentering = "primal"

        if centering == "primal":
            newCentering = "dual"

        return newCentering
