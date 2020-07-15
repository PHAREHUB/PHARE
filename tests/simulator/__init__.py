

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

    cells = [65 for i in range(dim)]
    if "cells" in kwargs:
        cells = kwargs["cells"]
    if not isinstance(cells, list):
        cells = [cells]
    dl = [1.0 / v for v in cells]
    b0 = [[20 for i in range(dim)], [25 for i in range(dim)]]
    boundary = ["periodic" for i in range(dim)]
    args = {
        "interp_order": interp,
        "smallest_patch_size": 5,
        "largest_patch_size": 64,
        "time_step_nbr": 1000,
        "final_time": 1.0,
        "boundary_types": boundary,
        "cells": cells,
        "dl": dl,
        "max_nbr_levels": 2,
        "refinement_boxes": {"L0": {"B0": b0}},
        "refined_particle_nbr": valid_refined_particle_nbr[dim][interp][0],
        "diag_options": {},
    }
    for k, v in kwargs.items():
        if k in args:
            args[k] = v
    return args


def pi_over_max_domain():
    return [np.pi / max_domain for max_domain in ph.globals.sim.simulation_domain()]


def fn_2d_periodic(sim, x, y):

    nx, ny = sim.cells
    xmax, ymax = sim.simulation_domain()
    dx, dy = (xmax / nx, ymax / ny)
    orix, oriy = sim.origin
    xx, yy = np.meshgrid(x,y,indexing="ij")
    kx, ky = 3, 6
    zx = np.cos(kx * 2 * np.pi / xmax * xx)
    zy = np.sin(ky * 2 * np.pi / ymax * yy)
    return zx * zy


def density_2d_periodic(sim, x, y):

    background_particles = 0.3  # avoids 0 density

    nx, ny = sim.cells
    xmax, ymax = sim.simulation_domain()
    dx, dy = (xmax / nx, ymax / ny)
    orix, oriy = sim.origin
    xx, yy = np.meshgrid(x,y,indexing="ij")

    return np.exp(-(xx-0.5*xmax)**2)*np.exp(-(yy-ymax/2.)**2) + background_particles


def defaultPopulationSettings():
    background_particles = 0.1  # avoids 0 density

    sim = ph.globals.sim
    dim = sim.dims

    xmax = sim.simulation_domain()[0]
    pi_over_xmax = pi_over_max_domain()[0]
    func_per_dim = {
        1: {
            "density": lambda x: 1.0 / np.cosh((x - xmax * 0.5)) ** 2
            + background_particles,
            "vbulkx": lambda x: np.sin(1 * pi_over_xmax * x),
            "vbulky": lambda x: np.sin(1 * pi_over_xmax * x),
            "vbulkz": lambda x: np.sin(1 * pi_over_xmax * x),
        },
        2: {
            "density": lambda x, y: density_2d_periodic(sim, x, y),
            "vbulkx": lambda x, y: fn_2d_periodic(sim, x, y) + 0.01,
            "vbulky": lambda x, y: fn_2d_periodic(sim, x, y) + 0.02,
            "vbulkz": lambda x, y: fn_2d_periodic(sim, x, y) + 0.03,
        },
    }
    assert dim in func_per_dim
    return {
        "charge": 1,
        "density": func_per_dim[dim]["density"],
        "vbulkx": func_per_dim[dim]["vbulkx"],
        "vbulky": func_per_dim[dim]["vbulky"],
        "vbulkz": func_per_dim[dim]["vbulkz"],
        "vthx": lambda *xyz: 1,
        "vthy": lambda *xyz: 1,
        "vthz": lambda *xyz: 1,
    }


def makeBasicModel(extra_pops={}):
    sim = ph.globals.sim
    dim = sim.dims
    pi_over_xmax = pi_over_max_domain()[0]
    func_per_dim = {
        1: {
            "bx": lambda x: np.cos(2 * pi_over_xmax * x),
            "by": lambda x: np.sin(1 * pi_over_xmax * x),
            "bz": lambda x: np.cos(2 * pi_over_xmax * x),
        },
        2: {
            "bx": lambda x, y: fn_2d_periodic(sim, x, y) + 0.04,
            "by": lambda x, y: fn_2d_periodic(sim, x, y) + 0.05,
            "bz": lambda x, y: fn_2d_periodic(sim, x, y) + 0.06,
        },
    }
    assert dim in func_per_dim

    pops = {
        "protons": {
            **defaultPopulationSettings(),
            "nbr_part_per_cell": 100,
            "init": {"seed": 1337},
        },
        "alpha": {
            **defaultPopulationSettings(),
            "nbr_part_per_cell": 100,
            "init": {"seed": 13337},
        },
    }
    pops.update(extra_pops)
    return ph.MaxwellianFluidModel(
        bx=func_per_dim[dim]["bx"],
        by=func_per_dim[dim]["by"],
        bz=func_per_dim[dim]["bz"],
        **pops
    )


def populate_simulation(dim, interp, **input):
    ph.globals.sim = None
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
