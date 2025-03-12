import os
import subprocess
import sys

# public export
from pyphare.pharein.load_balancer import LoadBalancer

from .diagnostics import (
    ElectromagDiagnostics,
    FluidDiagnostics,
    InfoDiagnostics,
    MetaDiagnostics,
    ParticleDiagnostics,
)
from .electron_model import ElectronModel
from .maxwellian_fluid_model import MaxwellianFluidModel
from .simulation import Simulation
from .uniform_model import UniformModel

__all__ = [
    "UniformModel",
    "MaxwellianFluidModel",
    "ElectronModel",
    "FluidDiagnostics",
    "ElectromagDiagnostics",
    "ParticleDiagnostics",
    "MetaDiagnostics",
    "InfoDiagnostics",
    "Simulation",
    "LoadBalancer",
]

# This exists to allow a condition variable for when we are running PHARE from C++ via phare-exe
#  It is configured to "True" in pyphare/pyphare/pharein/init.py::get_user_inputs(jobname)
#    which is called from Pybind when phare-exe (or similar) is in use
PHARE_EXE = False

venv_path = os.environ.get("VIRTUAL_ENV")

if venv_path is not None:
    pythonexe = os.path.join(venv_path, "bin/python3")
    arg = "import sys; print(sys.path)"
    p = subprocess.run([pythonexe, "-c", arg], stdout=subprocess.PIPE)
    s = p.stdout.decode()
    s = (
        s.replace("[", "")
        .replace("]", "")
        .replace('"', "")
        .replace("\n", "")
        .replace("'", "")
    )
    s = s.split(",")[1:]
    pythonpath = [ss.strip() for ss in s]
    sys.path = sys.path + pythonpath


def NO_GUI():
    """prevents issues when command line only and no desktop etc"""
    import matplotlib as mpl

    mpl.use("Agg")


def clearDict():
    """
    dict may contain dangling references from a previous simulation unless cleared
    """
    import pybindlibs.dictator as pp

    pp.stop()


def populateDict():
    from . import initialize
    from .global_vars import sim

    initialize.general.populateDict(sim)

    if not sim.init_options:
        sim.init_options = ["hybrid"]

    if "hybrid" in sim.init_options:
        initialize.hybrid.populateDict(sim)
    if "mhd" in sim.init_options:
        initialize.mhd.populateDict(sim)

    if not ("hybrid" in sim.init_options or "mhd" in sim.init_options):
        raise ValueError("Unknown simulation type")
