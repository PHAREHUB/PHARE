#
#
#

import os
import sys
import datetime
import atexit
import time as timem
from pathlib import Path

import numpy as np
import pyphare.pharein as ph

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


def make_cpp_simulator(
    hier,
    dim,
    interp,
    nbrRefinedPart,
    mhd_timestepper,
    reconstruction,
    limiter,
    riemann,
    hall,
    res,
    hyper_res,
):
    from pyphare.cpp import cpp_lib

    if SCOPE_TIMING:
        mon.timing_setup(cpp_lib())

    nbrRefinedPart_suffix = f"_{nbrRefinedPart}" if nbrRefinedPart else ""
    mhd_timestepper_suffix = f"_{mhd_timestepper}" if mhd_timestepper else ""
    reconstruction_suffix = f"_{reconstruction}" if reconstruction else ""
    limiter_suffix = f"_{limiter}" if limiter else ""
    riemann_suffix = f"_{riemann}" if riemann else ""
    hall_suffix = "_hall" if hall else ""
    res_suffix = "_res" if res else ""
    hyper_res_suffix = "_hyper_res" if hyper_res else ""

    make_sim = f"make_simulator_{dim}_{interp}{nbrRefinedPart_suffix}{mhd_timestepper_suffix}{reconstruction_suffix}{limiter_suffix}{riemann_suffix}{hall_suffix}{res_suffix}{hyper_res_suffix}"

    return getattr(cpp_lib(), make_sim)(hier)


def startMPI():
    if "samrai" not in life_cycles:
        from pyphare.cpp import cpp_lib

        life_cycles["samrai"] = cpp_lib().SamraiLifeCycle()


def print_rank0(*args, **kwargs):
    from pyphare.cpp import cpp_lib

    if cpp_lib().mpi_rank() == 0:
        print(*args, **kwargs)


def plot_timestep_time(timestep_times):
    from pyphare.cpp import cpp_lib

    if cpp_lib().mpi_rank() == 0:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(timestep_times)
        plt.ylabel("timestep time")
        plt.xlabel("timestep")
        fig.savefig("timestep_times.png")

    cpp_lib().mpi_barrier()


class Simulator:
    """

    **Mandatory arguments**

        *  **simulation**: a `Simulation` object


    **Optional expert arguments**

        These arguments have good default, change them at your own risk.

        *  **print_one_line**: (``bool``), default True, will print simulator info per advance on one line (erasing the previous)
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
        if kwargs.get("print_one_line", True):
            self.print_eol = "\r"
        self.print_eol = kwargs.get("print_eol", self.print_eol)
        self.log_to_file = kwargs.get("log_to_file", True)

        self.auto_dump = auto_dump
        import pyphare.simulator._simulator as _simulator

        _simulator.obj = self

    def __del__(self):
        self.reset()

    def setup(self):
        # mostly to detach C++ class construction/dict parsing from C++ Simulator::init
        try:
            import pyphare.cpp.validate as validate_cpp
            from pyphare.cpp import cpp_lib

            startMPI()

            if all([not self.simulation.dry_run, self.simulation.write_reports]):
                # not necessary during testing
                validate_cpp.log_runtime_config()
            validate_cpp.check_build_config_is_runtime_compatible()

            if self.log_to_file:
                self._log_to_file()
            ph.populateDict()
            self.cpp_hier = cpp_lib().make_hierarchy()

            refined_particle_nbr = getattr(
                self.simulation, "refined_particle_nbr", False
            )
            mhd_timestepper = getattr(self.simulation, "mhd_timestepper", False)
            reconstruction = getattr(self.simulation, "reconstruction", False)
            limiter = getattr(self.simulation, "limiter", False)
            riemann = getattr(self.simulation, "riemann", False)
            hall = getattr(self.simulation, "hall", False)
            res = getattr(self.simulation, "res", False)
            hyper_res = getattr(self.simulation, "hyper_res", False)

            self.cpp_sim = make_cpp_simulator(
                self.cpp_hier,
                self.simulation.ndim,
                self.simulation.interp_order,
                refined_particle_nbr,
                mhd_timestepper,
                reconstruction,
                limiter,
                riemann,
                hall,
                res,
                hyper_res,
            )
            return self
        except Exception:
            import traceback

            print('Exception caught in "Simulator.setup()": {}'.format(sys.exc_info()))
            print(traceback.extract_stack())
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

        try:
            self.cpp_sim.advance(dt)
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
        from pyphare.cpp import cpp_lib

        self._check_init()

        if monitoring is None:  # check env
            monitoring = SIM_MONITOR

        if self.simulation.dry_run:
            return self
        if monitoring:
            mon.setup_monitoring(cpp_lib())
        perf = []
        end_time = self.cpp_sim.endTime()
        t = self.cpp_sim.currentTime()

        while t < end_time:
            tick = timem.time()
            self.advance()
            tock = timem.time()
            ticktock = tock - tick
            perf.append(ticktock)
            t = self.cpp_sim.currentTime()
            if cpp_lib().mpi_rank() == 0:
                out = f"t = {t:8.5f}  -  {ticktock:6.5f}sec  - total {np.sum(perf):7.4}sec"
                print(out, end=self.print_eol)

        print_rank0(f"mean advance time = {np.mean(perf)}")
        print_rank0(f"total advance time = {datetime.timedelta(seconds=np.sum(perf))}")

        if plot_times:
            plot_timestep_time(perf)

        mon.monitoring_shutdown(cpp_lib())
        return self.reset()

    def _auto_dump(self):
        return self.auto_dump and self.dump()

    def dump(self, *args):
        assert len(args) == 0 or len(args) == 2

        if len(args) == 0:
            return self.cpp_sim.dump(
                timestamp=self.currentTime(), timestep=self.timeStep()
            )
        return self.cpp_sim.dump(timestamp=args[0], timestep=args[1])

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
        from pyphare.cpp import cpp_lib

        logging = os.environ["PHARE_LOG"] = os.environ.get("PHARE_LOG", "RANK_FILES")
        need_log_dir = logging != "CLI" and logging != "NULL"
        if need_log_dir and cpp_lib().mpi_rank() == 0:
            Path(".log").mkdir(exist_ok=True)
        cpp_lib().mpi_barrier()
