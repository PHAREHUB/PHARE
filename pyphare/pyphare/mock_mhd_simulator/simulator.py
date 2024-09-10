from . import clearDict
from . import populateDict

def pyMHD():
    import importlib

    return importlib.import_module("pybindlibs.pyMHD")

def make_cpp_simulator(dim, interp):
    import pybindlibs.pyMHD

    make_sim = f"make_mhd_mock_simulator_{dim}_{interp}"
    return getattr(pyMHD(), make_sim)()

class MHDMockSimulator:
    def __init__(self, simulation):
        self.cpp_sim = None
        self.simulation = simulation

    def __del__(self):
        self.reset()

    def reset(self):
        if self.cpp_sim is not None:
            clearDict()
        self.cpp_sim = None
        return self

    def initialize(self):
        if self.cpp_sim is not None:
            raise ValueError(
                "Simulator already initialized: requires reset to re-initialize"
            )

        populateDict()
        self.cpp_sim = make_cpp_simulator(
            self.simulation.ndim,
            self.simulation.order,
        )

    def _check_init(self):
        if self.cpp_sim is None:
            self.initialize()

    def run(self, filename, dumpfrequency=1):
        self._check_init()
        self.cpp_sim.advance(filename, dumpfrequency)
        return self
