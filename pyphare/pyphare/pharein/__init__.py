import os
import sys
import subprocess


from .diagnostics import (
    ElectromagDiagnostics,
    FluidDiagnostics,
    InfoDiagnostics,
    MetaDiagnostics,
    MHDDiagnostics,
    ParticleDiagnostics,
    FieldDiagnosticSlice,
)
from .load_balancer import LoadBalancer
from .electron_model import ElectronModel
from .maxwellian_fluid_model import MaxwellianFluidModel
from .mhd_model import MHDModel
from .simulation import Simulation
from .uniform_model import UniformModel

__all__ = [
    "UniformModel",
    "MaxwellianFluidModel",
    "ElectronModel",
    "FluidDiagnostics",
    "MHDModel",
    "MHDDiagnostics",
    "ElectromagDiagnostics",
    "ParticleDiagnostics",
    "MetaDiagnostics",
    "InfoDiagnostics",
    "Simulation",
    "LoadBalancer",
    "FieldDiagnosticSlice",
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


def add_vector_string(path, val):
    import pybindlibs.dictator as pp

    pp.add_vector_string(path, list(val))


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


def dict_populator():
    import pybindlibs.dictator as pp

    def add_size_t(path, val):
        casted = int(val)
        if casted < 0:
            raise RuntimeError("pyphare.__init__::add_size_t received negative value")
        pp.add_size_t(path, casted)

    def add_vector_int(path, val):
        pp.add_vector_int(path, list(val))

    class DictPopulator:
        def __init__(self):
            self.add_int = add_int
            self.add_bool = add_bool
            self.add_double = add_double
            self.add_size_t = add_size_t
            self.add_vector_int = add_vector_int
            self.add_string = pp.add_string

    return DictPopulator()


def populateDict(simulation=None):
    from . import initialize

    if simulation is None:
        from .global_vars import sim

        simulation = sim

    # dp = dict_populator()

    initialize.general.populateDict(simulation)
    add_vector_string("simulation/models", simulation.model_options)

    if "HybridModel" in simulation.model_options:
        initialize.hybrid.populateDict(simulation)
    if "MHDModel" in simulation.model_options:
        initialize.mhd.populateDict(simulation)

    # add_int = dp.add_int
    # add_bool = dp.add_bool
    # add_double = dp.add_double
    # add_size_t = dp.add_size_t
    # add_vector_int = dp.add_vector_int
    # add_string = dp.add_string
