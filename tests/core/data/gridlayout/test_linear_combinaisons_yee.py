#!/usr/bin/env pyhton
#!coding: utf-8

import numpy as np
import utilities
import os
import sys

# this script writes the following file
# in 1D, in 2D and in 3D :
# {dim} {interpOrder_i} ExToMoment
# {dim} {interpOrder_i} EyToMoment
# {dim} {interpOrder_i} EzToMoment
# {dim} {interpOrder_i} MomentToEx
# {dim} {interpOrder_i} MomentToEy
# {dim} {interpOrder_i} MomentToEz
# {dim} {interpOrder_i} ByToEx
# {dim} {interpOrder_i} BzToEx
# {dim} {interpOrder_i} BxToEy
# {dim} {interpOrder_i} BzToEy
# {dim} {interpOrder_i} BxToEz
# {dim} {interpOrder_i} ByToEz

# the script works only for the Yee Layout


def dualToPrimal(interpOrder):
    if interpOrder ==1 or interpOrder==2:
        return -1
    else:
        return 1


def primalToDual(interpOrder):
    if interpOrder == 1 or interpOrder == 2:
        return 1
    else:
        return -1


# each function returns the coordinates of points for linear interpolation
# there are 1, 2 or 3 integer coordinate for each point
# and the last number is a float coef for the linear combinaison


# in 1D there are either 2 int + 1 float (P1, P2, coef) or 1 int and a float (P1, coef=1)
# in 2D there are 4 int+ 1 float (P1, P2, coef) or 2 int + 1 float (P1 coef)
# in 3D there are 6 int + 1 float (P1, P2, coef) or 3 int and 1 float

# the file will be formated in the following way :
# example : 3 1 2 0 0 0 0 -1 0 0.5
# means 3D, interpOrder=1, 2 points, P1(0,0,0), P2(0,-1,0), coef=0.5


def tupleToString(t):
    return ' '.join(map(str,t))


def ExToMoment(dim, interpOrder):
    """
    Dpp to Ppp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,2, 0, dualToPrimal(interpOrder), 0.5))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0, 0, dualToPrimal(interpOrder),0, 0.5))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0,0,0, dualToPrimal(interpOrder), 0, 0, 0.5))


def EyToMoment(dim, interpOrder):
    """
    pDp to pPp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1,0, 1.))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0,0, 0, dualToPrimal(interpOrder), 0.5))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0,0,0, 0,dualToPrimal(interpOrder), 0, 0.5))


def EzToMoment(dim, interpOrder):
    """
    ppD to ppP
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1.))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0,0, 1))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0,0,0, 0,0,dualToPrimal(interpOrder), 0.5))


def JxToMoment(dim, interpOrder):
    """
    Dpp to Ppp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,2, 0, dualToPrimal(interpOrder), 0.5))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0, 0, dualToPrimal(interpOrder), 0, 0.5))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0, 0, 0, dualToPrimal(interpOrder), 0, 0, 0.5))


def JyToMoment(dim, interpOrder):
    """
    pDp to pPp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1.))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0, 0, 0, dualToPrimal(interpOrder), 0.5))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0, 0, 0, 0, dualToPrimal(interpOrder), 0, 0.5))


def JzToMoment(dim, interpOrder):
    """
    ppD to ppP
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1.))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0, 0, 1))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0, 0, 0, 0, 0, dualToPrimal(interpOrder), 0.5))


def MomentToEx(dim, interpOrder):
    """
    Ppp to Dpp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,2, 0, primalToDual(interpOrder), 0.5))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0,0, primalToDual(interpOrder), 0, 0.5))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0,0,0, primalToDual(interpOrder),0,0, 0.5))


def MomentToEy(dim, interpOrder):
    """
    pPp to pDp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0,0, 0, primalToDual(interpOrder), 0.5))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0,0,0, 0,primalToDual(interpOrder),0, 0.5))


def MomentToEz(dim, interpOrder):
    """
    ppP to ppD
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0,0, 1))
    elif dim == 3:
        return tupleToString((dim,interpOrder,2, 0,0,0,  0, 0,primalToDual(interpOrder), 0.5))


def ByToEx(dim, interpOrder):
    """
    dpD to dpP
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0,0, 1))
    elif dim==3:
        return tupleToString((dim,interpOrder,2, 0,0,0, 0,0,dualToPrimal(interpOrder), 0.5))


def BzToEx(dim, interpOrder):
    """
    dDp to dPp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0,0, 0, dualToPrimal(interpOrder), 0.5))
    elif dim==3:
        return tupleToString((dim,interpOrder,2, 0,0,0, 0,dualToPrimal(interpOrder),0, 0.5))


def BxToEz(dim, interpOrder):
    """
    pDd to pPd
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0,0, 0, dualToPrimal(interpOrder), 0.5))
    elif dim==3:
        return tupleToString((dim,interpOrder,2, 0,0,0, 0,dualToPrimal(interpOrder), 0, 0.5))


def ByToEz(dim, interpOrder):
    """
    Dpd to Ppd
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,2, 0, dualToPrimal(interpOrder), 0.5))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0,0, dualToPrimal(interpOrder),0, 0.5))
    elif dim==3:
        return tupleToString((dim,interpOrder,2, 0,0,0, dualToPrimal(interpOrder), 0, 0, 0.5))


def BxToEy(dim, interpOrder):
    """
    pdD to pdP
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0,0, 1))
    elif dim==3:
        return tupleToString((dim,interpOrder,2, 0,0,0, 0,0,dualToPrimal(interpOrder), 0.5))


def BzToEy(dim, interpOrder):
    """
    Ddp to PdP
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,2, 0, dualToPrimal(interpOrder), 0.5))
    elif dim == 2:
        return tupleToString((dim,interpOrder,2, 0, 0, dualToPrimal(interpOrder),0, 0.5))
    elif dim==3:
        return tupleToString((dim,interpOrder,2, 0,0,0, dualToPrimal(interpOrder), 0, 0, 0.5))


def JxToEx(dim, interpOrder):
    """
    dpp to dpp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1.0))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0,0, 1.0))
    elif dim==3:
        return tupleToString((dim,interpOrder,1, 0,0,0, 1.0))


def JyToEy(dim, interpOrder):
    """
    pdp to pdp
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1.0))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0,0, 1.0))
    elif dim==3:
        return tupleToString((dim,interpOrder,1, 0,0,0, 1.0))


def JzToEz(dim, interpOrder):
    """
    ppd to ppd
    """
    if dim == 1:
        return tupleToString((dim,interpOrder,1, 0, 1.0))
    elif dim == 2:
        return tupleToString((dim,interpOrder,1, 0,0, 1.0))
    elif dim==3:
        return tupleToString((dim,interpOrder,1, 0,0,0, 1.0))



def main(path='./'):

    if len(sys.argv) == 2:
        path = sys.argv[1]

    baseName = "linear_coefs_yee_"
    cases    = {"momentToEx":MomentToEx,
                "momentToEy":MomentToEy,
                "momentToEz":MomentToEz,
                "ExToMoment":ExToMoment,
                "EyToMoment":EyToMoment,
                "EzToMoment":EzToMoment,
                "JxToMoment":JxToMoment,
                "JyToMoment":JyToMoment,
                "JzToMoment":JzToMoment,
                "ByToEx":ByToEx,
                "BzToEx":BzToEx,
                "BxToEy":BxToEy,
                "BzToEy":BzToEy,
                "BxToEz":BxToEz,
                "ByToEz":ByToEz,
                "JxToEx":JxToEx,
                "JyToEy":JyToEy,
                "JzToEz":JzToEz}


    for case in cases:
        filename = baseName + case+".txt"
        outFile = open(os.path.join(path,filename), "w")

        for dim in [1,2, 3]:
            for interpOrder in [1,2,3]:
                outFile.write(cases[case](dim, interpOrder)+"\n")

        outFile.close()


if __name__ == '__main__':
    main()
