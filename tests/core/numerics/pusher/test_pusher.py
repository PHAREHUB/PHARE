#!/usr/bin/env python3
#!coding : utf-8

import os
import sys

from scipy.integrate import odeint
import numpy as np


def Bx(x,y,z):
    return 1.

def By(x,y,z):
    return 1.

def Bz(x,y,z):
    return 1.

def Ex(x,y,z):
    return 0.01

def Ey(x,y,z):
    return -0.05

def Ez(x,y,z):
    return 0.05

def deriv(vec, t0, q, m):

    x  = vec[0]
    y  = vec[1]
    z  = vec[2]
    vx = vec[3]
    vy = vec[4]
    vz = vec[5]

    return [vx,
            vy,
            vz,
            q/m * (vy*Bz(x,y,z) - vz*By(x,y,z) + Ex(x,y,z)),
            q/m * (vz*Bx(x,y,z) - vx*Bz(x,y,z) + Ey(x,y,z)),
            q/m * (vx*By(x,y,z) - vy*Bx(x,y,z) + Ez(x,y,z))
           ]


def main():
    path = sys.argv[1]

    rv0 = [0.25,0.25,0.25,0.,10.,0.]
    q   = 1
    m   = 1
    tstart = 0.
    tend   = 10.
    dt = 0.0001
    nt = int((tend-tstart)/dt)+1
    t      = np.arange(0,nt*dt, dt)

    sol = odeint(deriv, rv0, t, args=(q,m))
    print(path + os.path.sep + "pusher_test_in.txt")
    np.savetxt(path + os.path.sep + "pusher_test_in.txt", sol,delimiter=" ")


if __name__ == '__main__':
    main()
