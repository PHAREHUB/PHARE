#!/usr/bin/env python3

import math
import numpy as np

from .phare_utilities import listify
from .box import Box

directions = ["x", "y", "z"]
direction_to_dim = {direction: idx for idx, direction in enumerate(directions)}

yee_centering = {
    'x': {
        'Bx':'primal', 'By':'dual', 'Bz':'dual',
        'Ex':'dual', 'Ey':'primal', 'Ez':'primal',
        'Jx':'dual', 'Jy':'primal', 'Jz':'primal',
        'Vx':'primal','Vy':'primal','Vz':'primal',
        'Fx':'primal', 'Fy':'primal', 'Fz':'primal',
        'rho':'primal', 'P':'primal',
        'Vthx':'primal', 'Vthy':'primal', 'Vthz':'primal',
    },
    'y' : {
        'Bx':'dual', 'By':'primal', 'Bz':'dual',
        'Ex':'primal', 'Ey':'dual', 'Ez':'primal',
        'Jx':'primal', 'Jy':'dual', 'Jz':'primal',
        'Vx':'primal','Vy':'primal', 'Vz':'primal',
        'Fx':'primal', 'Fy':'primal', 'Fz':'primal',
        'rho':'primal', 'P':'primal',
        'Vthx':'primal', 'Vthy':'primal', 'Vthz':'primal',
    },
    'z' : {
        'Bx':'dual', 'By':'dual', 'Bz':'primal',
        'Ex':'primal', 'Ey':'primal', 'Ez':'dual',
        'Jx':'primal', 'Jy':'primal', 'Jz':'dual',
        'Vx':'primal','Vy':'primal', 'Vz':'primal',
        'Fx':'primal', 'Fy':'primal', 'Fz':'primal',
        'rho':'primal', 'P':'primal',
        'Vthx':'primal', 'Vthy':'primal', 'Vthz':'primal',
    }
}
yee_centering_lower = {
  direction : {
    key.lower() : centering
    for key, centering in key_dict.items()
  }
  for direction, key_dict in yee_centering.items()
}

def yee_element_is_primal(key, direction = 'x'):
    if key in yee_centering_lower[direction]:
        return yee_centering_lower[direction][key] == 'primal'
    return yee_centering[direction][key] == 'primal'


class YeeCentering(object):
    def __init__(self):
        self.centerX = yee_centering["x"]
        self.centerY = yee_centering["y"]
        self.centerZ = yee_centering["z"]

class GridLayout(object):

    def __init__(self, box=Box(0,0), origin=0, dl=0.1, interp_order=1):
        self.box = box

        self.dl = listify(dl)
        assert len(self.dl) == box.ndim

        self.origin = listify(origin)
        assert len(self.origin) == box.ndim

        self.interp_order = interp_order
        self.impl = "yee"
        self.directions = ['X','Y','Z']

        self.hybridQuantities = ['Bx','By','Bz',
                                 'Ex','Ey','Ez',
                                 'Jx','Jy','Jz',
                                 'rho','Vx','Vy',
                                 'Vz','P', 'Fx', 'Fy', 'Fz']


        self.yeeCentering = YeeCentering()

        self.centering = {'X' : self.yeeCentering.centerX,
                          'Y' : self.yeeCentering.centerY,
                          'Z' : self.yeeCentering.centerZ
                         }

    @property
    def ndim(self):
        return self.box.ndim

    def qtyCentering(self, quantity, direction):
        return self.centering[direction][quantity]


    def particleGhostNbr(self, interp_order):
        return 1 if interp_order == 1 else 2

    def nbrGhosts(self, interpOrder, centering):
        minNbrGhost = 5
        if centering == 'primal':
            if interpOrder == 1:
                return max(math.floor((interpOrder+1)/2), minNbrGhost)
            else:
                return max(math.floor( interpOrder/2 ), minNbrGhost)
        else:
            return max(math.floor( (interpOrder +1)/2 ), minNbrGhost)


    def nbrGhostsPrimal(self, interpOrder):
        minNbrGhost = 5
        if interpOrder == 1:
            return max(math.floor( (interpOrder+1)/2 ), minNbrGhost)
        else:
            return max(math.floor( interpOrder/2 ), minNbrGhost)



    def isDual(self, centering):
        if centering == 'dual':
            return 1
        else:
            return 0



    def ghostStartIndex(self):
        return 0;

    def ghostEndIndex(self, interpOrder, centering, nbrCells):
        index = self.physicalEndIndex(interpOrder, centering, nbrCells) \
              + self.nbrGhosts(interpOrder, centering)
        return index


    def physicalStartIndex(self, interpOrder, centering):
        index = self.ghostStartIndex() + self.nbrGhosts(interpOrder, centering)
        return index



    def physicalEndIndex(self, interpOrder, centering, nbrCells):
        index = self.physicalStartIndex(interpOrder, centering) \
                + nbrCells - self.isDual(centering)
        return index



    def physicalStartIndices(self, qty):
        assert qty in self.hybridQuantities
        indices = np.zeros(self.box.ndim)
        for i, direction in enumerate(directions[:self.box.ndim]):
            centering = yee_centering[direction][qty]
            indices[i] = self.physicalStartIndex(self.interp_order, centering)
        return indices


    def physicalEndIndices(self, qty):
        assert qty in self.hybridQuantities
        indices = np.zeros(self.box.ndim)
        for i, direction in enumerate(directions[:self.box.ndim]):
            centering = yee_centering[direction][qty]
            indices[i] = self.physicalEndIndex(self.interp_order, centering, self.box.shape[i])
        return indices


    def nbrGhostFor(self, qty):
        assert qty in self.hybridQuantities
        nGhosts = np.zeros(self.box.ndim)
        for i, direction in enumerate(directions[:self.box.ndim]):
            centering = yee_centering[direction][qty]
            nGhosts[i] = self.nbrGhosts(self.interp_order, centering)
        return nGhosts


    # ---- Start / End   primal methods ------
    def physicalStartPrimal(self, interpOrder):
        index = self.ghostStartIndex() + self.nbrGhostsPrimal(interpOrder)
        return index



    def physicalEndPrimal(self,interpOrder, nbrCells):
        index = self.physicalStartPrimal(interpOrder) + nbrCells
        return index

    # ---- Alloc methods -------------------------

    def allocSize(self, interpOrder, centering, nbrCells):
        size = nbrCells + 1 + 2*self.nbrGhosts(interpOrder, centering) \
               - self.isDual(centering)
        return size



    # 1st derivative
    def allocSizeDerived(self, interpOrder, centering, nbrCells):
        newCentering = self.changeCentering( centering, 1 )

        size = nbrCells + 1 + 2*self.nbrGhosts(interpOrder, newCentering) \
             - self.isDual(newCentering)
        return size


    def AMRIndexToLocal(self, dim, index):
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
        halfCell = 0.

        newCentering = self.changeCentering( centering, derivOrder )

        if newCentering == 'dual':
            halfCell = 0.5

        x = ( (knode - iStart) + halfCell )*ds + origin

        return x


    def yeeCoordsFor(self, qty, direction):
        """
        from a qty and a direction, returns a 1d array containing
        the coordinates where the qty is defined, including the ghost nodes

        :param qty: the quantity (can be primal or dual)
        :param direction: can only be a single one
        """

        assert direction in direction_to_dim, f"direction ({direction} not supported)"
        assert qty in yee_centering[direction] or qty in yee_centering_lower[direction], f"qty ({qty} not supported)"
        if qty in yee_centering_lower[direction] and qty not in yee_centering[direction]:
            qty = qty[0].upper() + qty[1:]

        centering = yee_centering[direction][qty]
        nbrGhosts = self.nbrGhosts(self.interp_order, centering)

        offset = 0
        dim = direction_to_dim[direction]
        size = self.box.shape[dim] + (nbrGhosts * 2)

        if centering == 'dual':
            offset = 0.5*self.dl[dim]
        else:
            size += 1

        return self.origin[dim] - nbrGhosts * self.dl[dim] + np.arange(size) * self.dl[dim] + offset





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
        halfCell = 0.

        newCentering = self.changeCentering( self.qtyCentering(qty, direction), derivOrder )

        if newCentering == 'dual':
            halfCell = 0.5

        x = ( (knode - iStart) + halfCell )*ds + origin

        return x



    # ---- Change centering method -------------------------
    #
    # Use case:
    #   changeCentering( qtyCentering(qty, direct), 1 )
    #
    def changeCentering(self, centering, derivOrder = 0):

        newCentering = centering

        # if derivOrder is odd the centering is changed
        if derivOrder % 2 != 0:
            newCentering = self.swapCentering( centering )

        return newCentering



    # -------------------------------------------------------
    def swapCentering(self, centering ):

        newCentering = 'primal'

        if centering == 'primal':
            newCentering = 'dual'

        return newCentering

