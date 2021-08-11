
import atexit
import time as timem
import numpy as np


life_cycles = {}

@atexit.register
def simulator_shutdown():
    from ._simulator import obj
    if obj is not None: # needs to be killed before MPI
        obj.reset()
    life_cycles.clear()


def make_cpp_simulator(dim, interp, nbrRefinedPart, hier, offloading):
    from pyphare.cpp import cpp_lib
    print(f"sim offload {offloading}")
    offload = 1 if offloading == "cuda" else 0
    make_sim = f"make_simulator_{dim}_{interp}_{nbrRefinedPart}"
    if offloading != None and offloading != 0:
        make_sim += f"_{offload}"
    print(f"make_sim {make_sim}")
    return getattr(cpp_lib(), make_sim)(hier)


def startMPI():
    if "samrai" not in life_cycles:
        from pyphare.cpp import cpp_lib
        life_cycles["samrai"] = cpp_lib().SamraiLifeCycle()


class Simulator:
    def __init__(self, simulation, auto_dump=True, **kwargs):
        import pyphare.pharein as ph #lgtm [py/import-and-import-from]
        assert isinstance(simulation, ph.Simulation)
        self.simulation = simulation
        self.cpp_hier = None   # HERE
        self.cpp_sim  = None   # BE
        self.cpp_dw   = None   # DRAGONS, i.e. use weakrefs if you have to ref these.
        self.post_advance = kwargs.get("post_advance", None)
        self.auto_dump = auto_dump
        import pyphare.simulator._simulator as _simulator
        _simulator.obj = self


    def __del__(self):
        self.reset()


    def initialize(self):
        if self.cpp_sim is not None:
            raise ValueError("Simulator already initialized: requires reset to re-initialize")
        try:
            from pyphare.cpp import cpp_lib
            from pyphare.pharein import populateDict
            startMPI()
            populateDict()
            self.cpp_hier = cpp_lib().make_hierarchy()
            sim = self.simulation
            self.cpp_sim = make_cpp_simulator(
              sim.ndim, sim.interp_order, sim.refined_particle_nbr, self.cpp_hier, sim.offload
            )

            self.cpp_sim.initialize()
            self._auto_dump() # first dump might be before first advance
            return self
        except:
            import sys
            print('Exception caught in "Simulator.initialize()": {}'.format(sys.exc_info()[0]))
            raise ValueError("Error in Simulator.initialize(), see previous error")

    def _throw(self, e):
        import sys
        from pyphare.cpp import cpp_lib
        if cpp_lib().mpi_rank() == 0:
            print(e)
        sys.exit(1)


    def advance(self, dt = None):
        self._check_init()
        if dt is None:
            dt = self.timeStep()
        
        try:
            self.cpp_sim.advance(dt)
        except (RuntimeError, TypeError, NameError, ValueError) as e:
            self._throw(f"Exception caught in simulator.py::advance: \n{e}")
        except KeyboardInterrupt as e:
            self._throw(f"KeyboardInterrupt in simulator.py::advance: \n{e}")

        if self._auto_dump() and self.post_advance != None:
            self.post_advance(self.cpp_sim.currentTime())
        return self

    def times(self):
        return np.arange(self.cpp_sim.startTime(),
                         self.cpp_sim.endTime() + self.timeStep(),
                         self.timeStep())

    def run(self):
        from pyphare.cpp import cpp_lib
        self._check_init()
        perf = []
        end_time = self.cpp_sim.endTime()
        t = 0.
        while t < end_time:
            tick  = timem.time()
            self.advance()
            tock = timem.time()
            ticktock = tock-tick
            perf.append(ticktock)
            t = self.cpp_sim.currentTime()
            if cpp_lib().mpi_rank() == 0:
                print("t = {:8.5f}  -  {:6.5f}sec  - total {:7.4}sec".format(t, ticktock, np.sum(perf)))

        print("mean advance time = {}".format(np.mean(perf)))
        print("total advance time = {}".format(np.sum(perf)))

        return self.reset()


    def _auto_dump(self):
        return self.auto_dump and len(self.simulation.diagnostics) > 0 and self.dump()


    def dump(self, *args):
        assert len(args) == 0 or len(args) == 2
        if len(args) == 0:
            return self.cpp_sim.dump(timestamp=self.currentTime(), timestep=self.timeStep())
        return self.cpp_sim.dump(timestamp=args[0], timestep=args[1])


    def data_wrangler(self):
        self._check_init()
        if self.cpp_dw is None:
            from pyphare.data.wrangler import DataWrangler
            self.cpp_dw = DataWrangler(self)
        return self.cpp_dw

    def reset(self):
        if self.cpp_dw is not None:
            self.cpp_dw.kill()
        self.cpp_dw = None
        self.cpp_sim = None
        self.cpp_hier = None
        if "samrai" in life_cycles:
            type(life_cycles["samrai"]).reset()
        return self

    def timeStep(self):
        self._check_init()
        return self.cpp_sim.timeStep()

    def currentTime(self):
        self._check_init()
        return self.cpp_sim.currentTime()

    def domain_box(self):
        self._check_init()
        return self.cpp_sim.domain_box()

    def cell_width(self):
        self._check_init()
        return self.cpp_sim.cell_width()

    def interp_order(self):
        self._check_init()
        return self.cpp_sim.interp_order # constexpr static value

    def _check_init(self):
        if self.cpp_sim is None:
            self.initialize()
