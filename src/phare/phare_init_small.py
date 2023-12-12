#!/usr/bin/env python3

import numpy as np

import pyphare.pharein as ph

# configure the simulation

ph.Simulation(
    smallest_patch_size=20,
    largest_patch_size=20,
    time_step_nbr=300,  # number of time steps (not specified if time_step and final_time provided)
    time_step=0.001,  # simulation final time (not specified if time_step and time_step_nbr provided)
    boundary_types="periodic",  # boundary condition, string or tuple, length == len(cell) == len(dl)
    cells=40,  # integer or tuple length == dimension
    dl=0.3,  # mesh size of the root level, float or tuple
    # max_nbr_levels=2,          # (default=1) max nbr of levels in the AMR hierarchy
    # refinement = "tagging",
    refinement_boxes={"L0": {"B0": [(10,), (20,)]}},
    diag_options={
        "format": "phareh5",
        "options": {"dir": "phare_outputs", "mode": "overwrite"},
    },
)


def density(x):
    return 1.0


def by(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()
    return 0.1 * np.cos(2 * np.pi * x / L[0])


def bz(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()
    return 0.1 * np.sin(2 * np.pi * x / L[0])


def bx(x):
    return 1.0


def vx(x):
    return 0.0


def vy(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()
    return 0.1 * np.cos(2 * np.pi * x / L[0])


def vz(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()
    return 0.1 * np.sin(2 * np.pi * x / L[0])


def vthx(x):
    return 0.01


def vthy(x):
    return 0.01


def vthz(x):
    return 0.01


vvv = {
    "vbulkx": vx,
    "vbulky": vy,
    "vbulkz": vz,
    "vthx": vthx,
    "vthy": vthy,
    "vthz": vthz,
}

ph.MaxwellianFluidModel(
    bx=bx,
    by=by,
    bz=bz,
    protons={"charge": 1, "density": density, **vvv, "init": {"seed": 1337}},
)

ph.ElectronModel(closure="isothermal", Te=0.12)


sim = ph.global_vars.sim

timestamps = np.arange(0, sim.final_time + sim.time_step, 10 * sim.time_step)


for quantity in ["E", "B"]:
    ph.ElectromagDiagnostics(
        quantity=quantity,
        write_timestamps=timestamps,
        compute_timestamps=timestamps,
    )


for quantity in ["density", "bulkVelocity"]:
    ph.FluidDiagnostics(
        quantity=quantity,
        write_timestamps=timestamps,
        compute_timestamps=timestamps,
    )

pops = [
    "protons",
]
