#!/usr/bin/env python3

import pyphare.pharein
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics
from pyphare.pharein import ElectronModel

# configure the simulation

Simulation(
    smallest_patch_size=10,
    largest_patch_size=64,
    time_step_nbr=1000,        # number of time steps (not specified if time_step and final_time provided)
    final_time=1.,             # simulation final time (not specified if time_step and time_step_nbr provided)
    boundary_types="periodic", # boundary condition, string or tuple, length == len(cell) == len(dl)
    cells=65,                  # integer or tuple length == dimension
    dl=1./65,                  # mesh size of the root level, float or tuple
    refinement_boxes = {"L0":{"B0":[(10,),(50,)]}},
    diag_options = {"format":"phareh5", "options": {"dir": "phare_outputs","mode":"overwrite"}}
)

density = lambda x: 2.

bx, by, bz = (lambda x: x  for i in range(3))
ex, ey, ez = (lambda x: x  for i in range(3))
vx, vy, vz = (lambda x: 1. for i in range(3))

vthx, vthy, vthz = (lambda x: 1. for i in range(3))

vvv = {
    "vbulkx":vx, "vbulky":vy, "vbulkz":vz,
    "vthx":vthx, "vthy":vthy, "vthz":vthz
}

MaxwellianFluidModel(
    bx=bx, by=by, bz=bz,
    protons={"charge":1, "density":density, **vvv, "init":{"seed":1337}},
    alpha={"charge":1, "density":density, **vvv, "init":{"seed":2}},
)

ElectronModel(closure="isothermal",Te = 0.12)


