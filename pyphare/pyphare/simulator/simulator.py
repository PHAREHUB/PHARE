

from pybindlibs import cpp
import pyphare.pharein as ph

life_cycles = {}

class Simulator:
    def __init__(self, simulation):
        assert isinstance(simulation, ph.Simulation)
        self.simulation = simulation
        self.cpp_dman = None
        self.cpp_sim  = None
        self.cpp_hier = None

    def __del__(self):
        self.reset()

    def initialize(self):
        if self.cpp_sim is not None:
            raise ValueError("Simulator already initialized: requires reset to re-initialize")
        try:
            if "samrai" not in life_cycles:
                life_cycles["samrai"] = cpp.SamraiLifeCycle()
            self.cpp_hier = cpp.make_hierarchy()
            self.cpp_sim = cpp.make_simulator(self.cpp_hier)
            self.cpp_sim.initialize()
        except:
            import sys
            print('Exception caught in "Simulator.initialize()": {}'.format(sys.exc_info()[0]))
            raise ValueError("no")

    def advance(self):
        self._check_init()
        self.cpp_sim.advance()

    def diagnostics(self):
        self._check_init()
        if self.cpp_dman is None:
            self.cpp_dman = cpp.make_diagnostic_manager(self.cpp_sim, self.cpp_hier)
        return self.cpp_dman

    def reset(self):
        self.cpp_dman = None
        self.cpp_sim = None
        self.cpp_hier = None
        if "samrai" in life_cycles:
            type(life_cycles["samrai"]).reset()

    def timeStep(self):
        self._check_init()
        return self.cpp_sim.timeStep()

    def currentTime(self):
        self._check_init()
        return self.cpp_sim.currentTime()

    def _check_init(self):
        if self.cpp_sim is None:
            raise ValueError("Simulator not initialized")