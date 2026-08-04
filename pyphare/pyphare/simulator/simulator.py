import atexit
import datetime
import os
import sys
import time as timem
from pathlib import Path

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.cpp import validate
from pyphare.pharein import restarts

from . import monitoring as mon

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

        *  **simulation**: a :class:`~pyphare.pharein.Simulation` object


    **Optional expert arguments**

        These arguments have good default, change them at your own risk.

        *  **print_one_line**: (``bool``), default False, will print simulator info per advance on one line (erasing the previous)
        *  **post_advance**: (``Function``),  default None. A python function called after diagnostics have
           been dumped, letting the user analyze the just-written data before the simulation takes another
           step
        *  **log_to_file**: if True (default), will log prints made from C++ code per MPI rank to the .log directory

    """

    def __init__(self, simulation, **kwargs):
        assert isinstance(simulation, ph.Simulation)  # pylint: disable=no-member
        self.simulation = simulation
        self.cpp_hier = None  # HERE
        self.cpp_sim = None  # BE
        self.cpp_dw = None  # DRAGONS, i.e. use weakrefs if you have to ref these.
        self.post_advance = kwargs.get("post_advance")
        self.initialized = False
        self.print_eol = "\n"
        if kwargs.get("print_one_line", False):
            self.print_eol = "\r"
        self.print_eol = kwargs.get("print_eol", self.print_eol)
        self.log_to_file = kwargs.get("log_to_file", True)

        from pyphare.simulator import _simulator

        _simulator.obj = self

    def __del__(self):
        self.reset()

    def setup(self):
        # mostly to detach C++ class construction/dict parsing from C++ Simulator::init
        try:
            startMPI()

            if all([not self.simulation.dry_run, self.simulation.write_reports]):
                # not necessary during testing
                validate.log_runtime_config()
            validate.check_build_config_is_runtime_compatible()

            if self.log_to_file:
                self._log_to_file()
            ph.populateDict()

            self.cpp_lib = cpp.cpp_lib(self.simulation)
            self.cpp_hier = cpp.cpp_etc_lib().make_hierarchy()
            self.cpp_sim = make_cpp_simulator(self.cpp_lib, self.cpp_hier)
            return self
        except Exception as e:  # see #1434
            import traceback

            print(f'Exception caught in "Simulator.setup()": {sys.exc_info()}')
            print(traceback.format_exc())
            raise ValueError("Error in Simulator.setup(), see previous error") from e

    def initialize(self):
        if self.initialized:
            return self

        self._check_setup()

        if self.simulation.dry_run:  # needs setup
            return self

        try:
            self.cpp_initialize()
            self.dump()  # first dump might be before first advance

            return self
        except Exception as e:  # see #1434
            print(f'Exception caught in "Simulator.initialize()": {sys.exc_info()[0]}')
            raise ValueError(
                "Error in Simulator.initialize(), see previous error"
            ) from e

    def advance(self):
        self._check_init()
        if self.simulation.dry_run:
            return self

        try:
            # dt lives on the C++ side: fixed by config, or recomputed per step when adaptive
            self.cpp_sim.advance()
        except (RuntimeError, TypeError, NameError, ValueError) as e:
            self._throw(f"Exception caught in simulator.py::advance: \n{e}")
        except KeyboardInterrupt as e:
            self._throw(f"KeyboardInterrupt in simulator.py::advance: \n{e}")

        if self.dump() and self.post_advance is not None:
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
        Run the simulation until the end time, calling :meth:`advance` in a loop.

        *  **plot_times**: (``bool``), default False. If True, plots (rank 0 only) the wall-clock time
           taken by each timestep once the run ends.
        *  **monitoring**: (``bool`` or ``int``), default None. Enables periodic resource monitoring
           (CPU/memory, per MPI rank) written to ``.phare/stats/rank.<rank>.yaml``, sampled every N
           seconds if given as an int (default 100 if just True), or falls back to the
           ``PHARE_SIM_MON`` environment variable if left None. Requires the optional ``phlop`` package
           (``python3 -m pip install phlop``); silently does nothing if it isn't installed.
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
        print_rank0("Starting at ", datetime.datetime.now().astimezone())
        print_rank0("Simulation start time/end time ", t, end_time)
        while t < end_time:
            tick = timem.time()
            self.advance()
            tock = timem.time()
            ticktock = tock - tick
            perf.append(ticktock)
            tot += ticktock
            t = self.cpp_sim.currentTime()
            if cpp.mpi_rank() == 0:
                delta = datetime.timedelta(seconds=tot)
                print(
                    f"t = {t:8.5f} - {ticktock:6.5f}sec - total {delta}",
                    end=self.print_eol,
                )

        if perf:
            print_rank0(f"mean advance time = {np.mean(perf)}")
        print_rank0(f"total advance time = {datetime.timedelta(seconds=tot)}")
        print_rank0("Finished at ", datetime.datetime.now().astimezone())

        if plot_times:
            plot_timestep_time(perf)

        mon.monitoring_shutdown()
        return self.reset()

    def dump(self, *args):
        # do diagnostics/restarts dump if active at this time
        if not (len(args) == 0 or len(args) == 2):
            raise ValueError("Simulator::dump: unexpected parameter count")

        time = self.currentTime() if len(args) == 0 else args[0]
        timestep = self.timeStep() if len(args) == 0 else args[1]
        restarts.dump(self, time, timestep)
        return self.cpp_sim.dump_diagnostics(timestamp=time, timestep=timestep)

    def data_wrangler(self):
        self._check_setup()
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
        self._check_setup()
        return self.cpp_sim.timeStep()

    def currentTime(self):
        self._check_setup()
        return self.cpp_sim.currentTime()

    def domain_box(self):
        self._check_setup()
        return self.cpp_sim.domain_box()

    def cell_width(self):
        self._check_setup()
        return self.cpp_sim.cell_width()

    def interp_order(self):
        self._check_setup()
        return self.cpp_sim.interp_order  # constexpr static value

    def cpp_initialize(self):
        # special case to initialize without dumping
        if self.initialized:
            return self
        self._check_setup()
        self.cpp_sim.initialize()
        self.initialized = True
        return self

    def _check_init(self):
        if not self.initialized:
            self.initialize()

    def _check_setup(self):
        if not self.cpp_sim:
            self.setup()

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

    def _throw(self, e):
        print_rank0(e)
        if exit_on_exception:
            sys.exit(1)
        # or reraise
        raise RuntimeError(e)
