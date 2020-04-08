from tests.simulator import test_simulator as tst


from phare import cpp
import unittest
import phare.pharein as ph, numpy as np, math
from phare.pharein import ElectronModel

from phare.pp.diagnostics import extract_diagnostics


class InitValueValidation(unittest.TestCase):
    diag_options = lambda diag_out_dir: {
        "diag_options": {"format": "phareh5", "options": {"dir": diag_out_dir},}
    }

    def getSimulation(self):
        return ph.globals.sim

    def runAndDump(self, dim, interp, input):
        self.dman, self.sim, self.hier = create_simulator(dim, interp, **input)
        timestamp = 0
        timestep = 1
        self.dman.dump(timestamp, timestep)
        del self.dman, self.sim, self.hier  # force hdf5 flush
        return extract_diagnostics(ph.globals.sim.diag_options["options"]["dir"])

    def tearDown(self):
        for k in ["dman", "sim", "hier"]:
            if hasattr(self, k):
                delattr(self, k)
        cpp.reset()


def basicSimulatorArgs(dim: int, interp: int, **kwargs):
    def _checkSetList(name, val):
        if name in kwargs:
            val = kwargs[name]
        if not isinstance(val, list):
            val = [val]
        kwargs[name] = val
        return val

    cells = _checkSetList("cells", [65 for i in range(dim)])
    origin = _checkSetList("origin", [0 for i in range(dim)])

    dl = [1.0 / v for v in cells]
    args = {
        "interp_order": interp,
        "smallest_patch_size": 5,
        "largest_patch_size": 64,
        "time_step_nbr": 1000,
        "final_time": 1.0,
        "boundary_types": "periodic",
        "cells": cells,
        "origin": origin,
        "dl": dl,
        "max_nbr_levels": 2,
        "refinement_boxes": {"L0": {"B0": [(10,), (50,)]}},
        "diag_options": {},
    }
    for k, v in kwargs.items():
        if k in args:
            args[k] = v

    return args


def defaultPopulationSettings():
    background_particles = 0.1  # avoids 0 density
    xmax = ph.globals.sim.simulation_domain()[0]
    pi_over_xmax = np.pi / xmax
    return {
        "charge": 1,
        "density": lambda x: 1.0 / np.cosh((x - xmax * 0.5)) ** 2
        + background_particles,
        "vbulkx": lambda x: np.sin(1 * pi_over_xmax * x),
        "vbulky": lambda x: np.sin(1 * pi_over_xmax * x),
        "vbulkz": lambda x: np.sin(1 * pi_over_xmax * x),
        "vthx": lambda x: 1,
        "vthy": lambda x: 1,
        "vthz": lambda x: 1,
    }


def makeBasicModel(extra_pops={}):
    xmax = ph.globals.sim.simulation_domain()[0]
    pi_over_xmax = np.pi / xmax
    pops = {
        "protons": {
            **defaultPopulationSettings(),
            "nbr_part_per_cell": 100,
            "init": {"seed": 1337},
        },
        "alpha": {
            **defaultPopulationSettings(),
            "nbr_part_per_cell": 1000,
            "init": {"seed": 13337},
        },
    }
    pops.update(extra_pops)
    return ph.MaxwellianFluidModel(
        bx=lambda x: np.cos(2 * pi_over_xmax * x),
        by=lambda x: np.sin(1 * pi_over_xmax * x),
        bz=lambda x: np.cos(2 * pi_over_xmax * x),
        ex=lambda x: np.sin(1 * pi_over_xmax * x),
        ey=lambda x: np.cos(2 * pi_over_xmax * x),
        ez=lambda x: np.sin(1 * pi_over_xmax * x),
        **pops
    )


def create_simulator(dim, interp, **input):

    cpp.reset()
    ph.globals.sim = None
    ph.Simulation(**basicSimulatorArgs(dim, interp, **input))
    extra_pops = {}
    if "populations" in input:
        for pop, vals in input["populations"].items():
            extra_pops[pop] = defaultPopulationSettings()
            extra_pops[pop].update(vals)

    model = makeBasicModel(extra_pops)
    if "diags_fn" in input:
        input["diags_fn"](model)

    ElectronModel(closure="isothermal", Te=0.12)

    ph.populateDict()
    hier = cpp.make_hierarchy()
    sim = cpp.make_simulator(hier)
    sim.initialize()
    return [cpp.make_diagnostic_manager(sim, hier), sim, hier]
