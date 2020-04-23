#!/usr/bin/env python3

import math
from .phare_utilities import listify
from ..pharesee.box import Box

yee_centering = {
    'x': {
        'Bx':'primal', 'By':'dual', 'Bz':'dual',
        'Ex':'dual', 'Ey':'primal', 'Ez':'primal',
        'Jx':'dual', 'Jy':'primal', 'Jz':'primal',
        'rho':'primal', 'Vx':'primal','Vy':'primal',
        'Vz':'primal','P':'primal'
    },
    'y' : {
        'Bx':'dual', 'By':'primal', 'Bz':'dual',
        'Ex':'primal', 'Ey':'dual', 'Ez':'primal',
        'Jx':'primal', 'Jy':'dual', 'Jz':'primal',
        'rho':'primal', 'Vx':'primal','Vy':'primal',
        'Vz':'primal','P':'primal'
    },
    'z' : {
        'Bx':'dual', 'By':'dual', 'Bz':'primal',
        'Ex':'primal', 'Ey':'primal', 'Ez':'dual',
        'Jx':'primal', 'Jy':'primal', 'Jz':'dual',
        'rho':'primal', 'Vx':'primal','Vy':'primal',
        'Vz':'primal', 'P':'primal'
    }
}
yee_centering_lower = {
  direction : {
    key.lower() : prumal
    for key, prumal in key_dict.items()
  }
  for direction, key_dict in yee_centering.items()
}

def yee_element_is_primal(key, direction = 'x'):
    if key in yee_centering_lower[direction]:
        return yee_centering_lower[direction][key] is 'primal'
    return yee_centering[direction][key] is 'primal'


class YeeCentering(object):
    def __init__(self):
        self.centerX = yee_centering["x"]
        self.centerY = yee_centering["y"]
        self.centerZ = yee_centering["z"]

class GridLayout(object):

    def __init__(self, box=Box(0,0), origin=0, dl=0.1, interp_order=1):
        self.box = box
        self.origin = listify(origin)
        self.dl = listify(dl)
        self.interp_order = interp_order

        self.directions = ['X','Y','Z']

        self.hybridQuantities = ['Bx','By','Bz',
                                 'Ex','Ey','Ez',
                                 'Jx','Jy','Jz',
                                 'rho','Vx','Vy',
                                 'Vz','P']


        self.yeeCentering = YeeCentering()

        self.centering = {'X' : self.yeeCentering.centerX,
                          'Y' : self.yeeCentering.centerY,
                          'Z' : self.yeeCentering.centerZ
                         }

    def qtyCentering(self, quantity, direction):
        return self.centering[direction][quantity]


    def particleGhostNbr(self, interp_order):
        return 1 if interp_order == 1 else 2

    def nbrGhosts(self,interpOrder, centering):
        minNbrGhost = 5
        if centering == 'primal':
            if interpOrder == 1:
                return max(math.floor((interpOrder+1)/2), minNbrGhost)
            else:
                return max(math.floor( interpOrder/2 ), minNbrGhost)
        else:
            return max(math.floor( (interpOrder +1)/2 ), minNbrGhost)


    def nbrGhostsPrimal(self,interpOrder):
        minNbrGhost = 5
        if interpOrder == 1:
            return max(math.floor( (interpOrder+1)/2 ), minNbrGhost)
        else:
            return max(math.floor( interpOrder/2 ), minNbrGhost)



    def isDual(self,centering):
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


    def physicalStartIndex(self,interpOrder, centering):
        index = self.ghostStartIndex() + self.nbrGhosts(interpOrder, centering)
        return index



    def physicalEndIndex(self, interpOrder, centering, nbrCells):
        index = self.physicalStartIndex(interpOrder, centering) \
        + nbrCells - self.isDual(centering)
        return index




    # ---- Start / End   primal methods ------
    def physicalStartPrimal(self,interpOrder):
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
    def changeCentering(self, centering, derivOrder ):

    #    print( "Inputs\n   centering : %s\n   derivOrder : %d" % (centering, derivOrder) )

        newCentering = centering

        # if derivOrder is odd the centering
        # is changed
        if derivOrder%2 != 0:
            newCentering = self.swapCentering( centering )

    #    print( "   output : %s" % newCentering )

        return newCentering



    # -------------------------------------------------------
    def swapCentering(self, centering ):

        newCentering = 'primal'

        if centering == 'primal':
            newCentering = 'dual'

        return newCentering


