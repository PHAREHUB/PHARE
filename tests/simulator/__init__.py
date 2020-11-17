

import pyphare.pharein as ph, numpy as np
from pyphare.pharein import ElectronModel

# Block accidental dictionary key rewrites
class NoOverwriteDict(dict):
    def __init__(self, dict):
        for k, v in dict.items():
            self[k] = v

    def __setitem__(self, k, v):
        if k in self.keys():
            raise ValueError("Key is already present")
        else:
            return super(NoOverwriteDict, self).__setitem__(k, v)


def basicSimulatorArgs(dim: int, interp: int, **kwargs):
    from pyphare.pharein.simulation import valid_refined_particle_nbr

    cells = [20 for i in range(dim)]
    if "cells" in kwargs:
        cells = kwargs["cells"]
    if not isinstance(cells, list):
        cells = [cells]
    dl = [1.0 / v for v in cells]
    b0 = [[3 for i in range(dim)], [8 for i in range(dim)]]
    boundary = ["periodic" for i in range(dim)]
    args = {
        "interp_order": interp,
        "smallest_patch_size": 5,
        "largest_patch_size": 20,
        "time_step_nbr": 1000,
        "final_time": 1.0,
        "boundary_types": boundary,
        "cells": cells,
        "dl": dl,
        "refinement_boxes": {"L0": {"B0": b0}},
        "refined_particle_nbr": valid_refined_particle_nbr[dim][interp][0],
        "diag_options": {},
        "nesting_buffer": 0,
    }
    for k, v in kwargs.items():
        if k in args:
            args[k] = v

    return args

def meshify(*xyz):
    if all([isinstance(v, np.ndarray) for v in xyz]):
        return xyz
    return np.meshgrid(*xyz,indexing="ij")

def pi_over_max_domain():
    return [np.pi / max_domain for max_domain in ph.global_vars.sim.simulation_domain()]

def fn_1d_periodic(sim, x):
    pi_over_xmax = pi_over_max_domain()[0]
    return np.sin(1 * pi_over_xmax * x)

def fn_2d_periodic(sim, x, y):
    xmax, ymax = sim.simulation_domain()
    xx, yy = meshify(x, y)
    kx, ky = 3, 6
    zx = np.cos(kx * 2 * np.pi / xmax * xx)
    zy = np.sin(ky * 2 * np.pi / ymax * yy)
    return zx * zy

# def fn_3d_periodic(sim, x, y, z):
#     xmax, ymax, zmax = sim.simulation_domain()
#     xx, yy, zz = meshify(x, y, z)
#     kx, ky, kz = 3, 6, 3
#     zx = np.cos(kx * 2 * np.pi / xmax * xx)
#     zy = np.sin(ky * 2 * np.pi / ymax * yy)
#     zz = np.sin(kz * 2 * np.pi / zmax * zz)
#     r = zx * zy * zz
#     return r


def density_1d_periodic(sim, x):
    xmax = sim.simulation_domain()[0]
    background_particles = 0.3  # avoids 0 density
    return 1.0 / np.cosh((x - xmax * 0.5) ** 2 + background_particles)

def density_2d_periodic(sim, x, y):
    xmax, ymax = sim.simulation_domain()
    background_particles = 0.3  # avoids 0 density
    xx, yy = meshify(x, y)
    return np.exp(-(xx-0.5*xmax)**2)*np.exp(-(yy-ymax/2.)**2) + background_particles

# def density_3d_periodic(sim, x, y, z):
#     xmax, ymax, zmax = sim.simulation_domain()
#     background_particles = 0.3  # avoids 0 density
#     xx, yy, zz = meshify(x, y, z)
#     r = np.exp(-(xx-0.5*xmax)**2)*np.exp(-(yy-ymax/2.)**2)*np.exp(-(zz-zmax/2.)**2) + background_particles
#     return r


def defaultPopulationSettings(sim, density_fn, vbulk_fn):
    return {
        "charge": 1,
        "density": lambda *xyz: density_fn(sim, *xyz),
        "vbulkx": lambda *xyz: vbulk_fn(sim, *xyz) + 0.01,
        "vbulky": lambda *xyz: vbulk_fn(sim, *xyz) + 0.02,
        "vbulkz": lambda *xyz: vbulk_fn(sim, *xyz) + 0.03,
        "vthx": lambda *xyz: 1,
        "vthy": lambda *xyz: 1,
        "vthz": lambda *xyz: 1,
    }


def makeBasicModel(extra_pops={}):
    sim = ph.global_vars.sim
    _density_fn_periodic = globals()["density_"+str(sim.dims)+"d_periodic"]
    _fn_periodic = globals()["fn_"+str(sim.dims)+"d_periodic"]

    pops = {
        "protons": {
            **defaultPopulationSettings(sim, _density_fn_periodic, _fn_periodic),
            "nbr_part_per_cell": 100,
            "init": {"seed": 1337},
        },
        "alpha": {
            **defaultPopulationSettings(sim, _density_fn_periodic, _fn_periodic),
            "nbr_part_per_cell": 100,
            "init": {"seed": 13337},
        },
    }
    pops.update(extra_pops)
    return ph.MaxwellianFluidModel(
        bx= lambda *xyz: _fn_periodic(sim, *xyz) + 0.04,
        by= lambda *xyz: _fn_periodic(sim, *xyz) + 0.05,
        bz= lambda *xyz: _fn_periodic(sim, *xyz) + 0.06,
        **pops
    )


def populate_simulation(dim, interp, **input):
    ph.global_vars.sim = None
    simulation = ph.Simulation(**basicSimulatorArgs(dim, interp, **input))
    extra_pops = {}
    if "populations" in input:
        for pop, vals in input["populations"].items():
            extra_pops[pop] = defaultPopulationSettings()
            extra_pops[pop].update(vals)

    model = makeBasicModel(extra_pops)
    if "diags_fn" in input:
        input["diags_fn"](model)

    ElectronModel(closure="isothermal", Te=0.12)

    return simulation
