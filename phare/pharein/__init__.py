
from .uniform_model import UniformModel
from .maxwellian_fluid_model import MaxwellianFluidModel
from .diagnostics import DiagnosticInfo, FluidDiagnostics, ElectromagDiagnostics, ParticleDiagnostics
from .simulation import Simulation
from . import globals

import os


def prepare_job():
    """
    prepare a simulation for a run:
        - creates the simulation directory [simulation.path] if it does not exist yet
        - writes an INI file in the directory [simulation.path]
        - for each diagnostic registered in the simulation:
            - creates the diagnostic directory 'diag['path']' under [simulation.path]/
              if it does not exist yet
    """

    if globals.sim is None:
        raise RuntimeError("No Simulation created")

    sim = globals.sim

    print("preparing job...")
    if not os.path.exists(sim.path):
        print("mkdir "+sim.path)
        os.makedirs(sim.path)
    print("writing ini file in " + sim.path + os.path.sep)
    sim.write_ini_file()

    for diag in sim.diagnostics:
        path = diag.path
        full_path = os.path.join(sim.path, path)
        if not os.path.exists(full_path):
            print("mkdir " + full_path)
            os.makedirs(full_path)
