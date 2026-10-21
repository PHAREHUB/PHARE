#!/usr/bin/env python3

import numpy as np

import pyphare.pharein as ph
from pyphare.core.box import Box2D

cells = (40, 40)

ph.Simulation(
    time_step=0.001,
    time_step_nbr=10,
    cells=cells,
    dl=(1.0 / cells[0], 1.0 / cells[1]),
    # level 2 two level-1 cells inside level 1, so its ghosts reach level 1 only
    refinement_boxes={"L0": {"B0": Box2D(10, 29)}, "L1": {"B0": Box2D(22, 57)}},
    smallest_patch_size=6,
    largest_patch_size=10,
    nesting_buffer=1,
    hyper_mode="spatial",
    eta=0.0,
    nu=0.0,
    gamma=5.0 / 3.0,
    reconstruction="WENOZ",
    limiter="None",
    riemann="Rusanov",
    mhd_timestepper="TVDRK3",
    hall=False,
    model_options=["MHDModel"],
    max_mhd_level=3,
)

# Orszag-Tang
B0 = 1.0 / np.sqrt(4.0 * np.pi)


def density(x, y):
    return 25.0 / (36.0 * np.pi) + 0.0 * x


def vx(x, y):
    return -np.sin(2.0 * np.pi * y)


def vy(x, y):
    return np.sin(2.0 * np.pi * x)


def bx(x, y):
    return -B0 * np.sin(2.0 * np.pi * y)


def by(x, y):
    return B0 * np.sin(4.0 * np.pi * x)


def zero(x, y):
    return 0.0 * x


def p(x, y):
    return 5.0 / (12.0 * np.pi) + 0.0 * x


ph.MHDModel(density=density, vx=vx, vy=vy, vz=zero, bx=bx, by=by, bz=zero, p=p)
