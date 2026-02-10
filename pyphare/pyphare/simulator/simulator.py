#
#
#

import os
import sys
import datetime
import atexit
import time as timem
import numpy as np
import pyphare.pharein as ph
from pathlib import Path
from . import monitoring as mon

from pyphare import cpp
import pyphare.pharein.restarts as restarts


exit_on_exception = True
life_cycles = {}
SIM_MONITOR = os.getenv("PHARE_SIM_MON", "False").lower() in ("true", "1", "t")
SCOPE_TIMING = os.getenv("PHARE_SCOPE_TIMING", "False").lower() in ("true", "1", "t")


@atexit.register
def simulator_shutdown():
    from ._simulator import obj

    if obj is not None:  # needs to be killed before MPI
        obj.reset()
    life_cycles.clear()


def make_cpp_simulator(cpp_lib, hier):
    if SCOPE_TIMING:
        mon.timing_setup()

    make_sim = "make_simulator"
    assert hasattr(cpp_lib, make_sim)
    return getattr(cpp_lib, make_sim)(hier)


def startMPI():
    if "samrai" not in life_cycles:
        life_cycles["samrai"] = cpp.cpp_etc_lib().SamraiLifeCycle()


def print_rank0(*args, **kwargs):
    if cpp.mpi_rank() == 0:
        print(*args, **kwargs)


def plot_timestep_time(timestep_times):
    if cpp.mpi_rank() == 0:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(timestep_times)
        plt.ylabel("timestep time")
        plt.xlabel("timestep")
        fig.savefig("timestep_times.png")

    cpp.mpi_barrier()


class Simulator:
    """

    **Mandatory arguments**

        *  **simulation**: a `Simulation` object


    **Optional expert arguments**

        These arguments have good default, change them at your own risk.

        *  **print_one_line**: (``bool``), default False, will print simulator info per advance on one line (erasing the previous)
        *  **auto_dump**: (``bool``), if True (default), will dump diagnostics automatically at requested timestamps
        *  **post_advance**: (``Function``),  default None. A python function to execute after each advance()
        *  **log_to_file**: if True (default), will log prints made from C++ code per MPI rank to the .log directory

    """

    def __init__(self, simulation, auto_dump=True, **kwargs):
        assert isinstance(simulation, ph.Simulation)  # pylint: disable=no-member
        self.simulation = simulation
        self.cpp_hier = None  # HERE
        self.cpp_sim = None  # BE
        self.cpp_dw = None  # DRAGONS, i.e. use weakrefs if you have to ref these.
        self.post_advance = kwargs.get("post_advance", None)
        self.initialized = False
        self.print_eol = "\n"
        if kwargs.get("print_one_line", False):
            self.print_eol = "\r"
        self.print_eol = kwargs.get("print_eol", self.print_eol)
        self.log_to_file = kwargs.get("log_to_file", True)
        self.report = ""
        self.auto_dump = auto_dump
        import pyphare.simulator._simulator as _simulator

        _simulator.obj = self

    def __del__(self):
        self.reset()

    def setup(self):
        # mostly to detach C++ class construction/dict parsing from C++ Simulator::init
        try:
            startMPI()

            if all([not self.simulation.dry_run, self.simulation.write_reports]):
                # not necessary during testing
                cpp.validate.log_runtime_config()
            cpp.validate.check_build_config_is_runtime_compatible()

            if self.log_to_file:
                self._log_to_file()
            ph.populateDict(self.simulation)

            self.cpp_lib = cpp.cpp_lib(self.simulation)
            self.cpp_hier = cpp.cpp_etc_lib().make_hierarchy()
            self.cpp_sim = make_cpp_simulator(self.cpp_lib, self.cpp_hier)

            return self
        except Exception:
            import traceback

            print('Exception caught in "Simulator.setup()": {}'.format(sys.exc_info()))
            print(traceback.format_exc())
            raise ValueError("Error in Simulator.setup(), see previous error")

    def initialize(self):
        try:
            if self.initialized:
                return
            if self.cpp_hier is None:
                self.setup()

            if self.simulation.dry_run:
                return self

            self.cpp_sim.initialize()
            self.initialized = True
            self._auto_dump()  # first dump might be before first advance

            return self
        except Exception:
            print(
                'Exception caught in "Simulator.initialize()": {}'.format(
                    sys.exc_info()[0]
                )
            )
            raise ValueError("Error in Simulator.initialize(), see previous error")

    def _throw(self, e):
        print_rank0(e)
        if exit_on_exception:
            sys.exit(1)
        # or reraise
        raise RuntimeError(e)

    def advance(self, dt=None):
        self._check_init()
        if self.simulation.dry_run:
            return self
        if dt is None:
            dt = self.timeStep()

        self.report = ""
        try:
            self.report = self.cpp_sim.advance(dt)
        except (RuntimeError, TypeError, NameError, ValueError) as e:
            self._throw(f"Exception caught in simulator.py::advance: \n{e}")
        except KeyboardInterrupt as e:
            self._throw(f"KeyboardInterrupt in simulator.py::advance: \n{e}")

        if self._auto_dump() and self.post_advance is not None:
            self.post_advance(self.cpp_sim.currentTime())
        return self

    def times(self):
        return np.arange(
            self.cpp_sim.startTime(),
            self.cpp_sim.endTime() + self.timeStep(),
            self.timeStep(),
        )

    def run(self, plot_times=False, monitoring=None):
        """
        Run the simulation until the end time
        monitoring requires phlop
        """

        self._check_init()

        if monitoring is None:  # check env
            monitoring = SIM_MONITOR
        if self.simulation.dry_run:
            return self
        if monitoring:
            interval = monitoring if isinstance(monitoring, int) else 100  # seconds
            mon.setup_monitoring(interval)
        perf = []
        end_time = self.cpp_sim.endTime()
        t = self.cpp_sim.currentTime()

        tot = 0
        print_rank0("Starting at ", datetime.datetime.now())
        print_rank0("Simulation start time/end time ", t, end_time)
        while t < end_time:
            tick = timem.time()
            self.advance()
            tock = timem.time()
            ticktock = tock - tick
            perf.append(ticktock)
            tot += ticktock
            t = self.cpp_sim.currentTime()
            delta = datetime.timedelta(seconds=tot)
            if cpp.mpi_rank() == 0:
                print(
                    f"t = {t:8.5f} - {ticktock:6.5f}sec - total {delta} {self.report}",
                    end=self.print_eol,
                )

        print_rank0(f"mean advance time = {np.mean(perf)}")
        print_rank0(f"total advance time = {datetime.timedelta(seconds=np.sum(perf))}")
        print_rank0("Finished at ", datetime.datetime.now())

        if plot_times:
            plot_timestep_time(perf)

        mon.monitoring_shutdown()
        return self.reset()

    def _auto_dump(self):
        return self.auto_dump and self.dump()

    def dump(self, *args):
        assert len(args) == 0 or len(args) == 2

        time = self.currentTime() if len(args) == 0 else args[0]
        timestep = self.timeStep() if len(args) == 0 else args[1]

        restarts.dump(self, time, timestep)

        return self.cpp_sim.dump_diagnostics(timestamp=time, timestep=timestep)

    def data_wrangler(self):
        self._check_init()
        if self.cpp_dw is None:
            from pyphare.data.wrangler import DataWrangler

            self.cpp_dw = DataWrangler(self)
        return self.cpp_dw

    def reset(self):
        if self.cpp_sim is not None:
            ph.clearDict()
        if self.cpp_dw is not None:
            self.cpp_dw.kill()
        self.cpp_dw = None
        self.cpp_sim = None
        self.cpp_hier = None
        self.initialized = False
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
        return self.cpp_sim.interp_order  # constexpr static value

    def _check_init(self):
        if not self.initialized:
            self.initialize()

    def _log_to_file(self):
        """
        send C++ std::cout logs to files with env var PHARE_LOG
        Support keys:
            RANK_FILES - logfile per rank
            DATETIME_FILES - logfile with starting datetime timestamp per rank
            CLI  - no logging files, display to cout
            NULL - no logging files, no cout
        """

        logging = os.environ["PHARE_LOG"] = os.environ.get("PHARE_LOG", "RANK_FILES")
        need_log_dir = logging != "CLI" and logging != "NULL"
        if need_log_dir and cpp.mpi_rank() == 0:
            Path(".log").mkdir(exist_ok=True)
        cpp.mpi_barrier()
