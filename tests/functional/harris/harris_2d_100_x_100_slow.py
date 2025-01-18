#!/usr/bin/env python3

import os
import numpy as np
import matplotlib as mpl
from pathlib import Path

import pyphare.pharein as ph
from pyphare.cpp import cpp_lib
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest
from tools.python3 import plotting as m_plotting

import harris_2d as base

mpl.use("Agg")

SCOPE_TIMING = os.getenv("PHARE_SCOPE_TIMING", "False").lower() in ("true", "1", "t")
"""
  For scope timings to work
  The env var PHARE_SCOPE_TIMING must be == "1" (or "true")
    See src/phare/phare.hpp
  CMake must be configured with: -DwithPhlop=ON
  And a LOG_LEVEL must be defined via compile args: -DPHARE_LOG_LEVEL=1
  Or change the default value in src/core/logger.hpp
  And phlop must be available on PYTHONPATH either from subprojects
   or install phlop via pip
"""

LOAD_BALANCE = os.getenv("LOAD_BALANCE", "True").lower() in ("true", "1", "t")

cpp = cpp_lib()
startMPI()

cells = (100, 100)
final_time = 50
time_step = 0.001
timestamps = np.arange(0, final_time + time_step, final_time / 5)

diag_dir = "phare_outputs/harris_2d_100_x_100_slow"
plot_dir = Path(f"{diag_dir}_plots")
plot_dir.mkdir(parents=True, exist_ok=True)


def config():
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=(0.40, 0.40),
        # refinement="tagging",
        # max_nbr_levels=1,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="1",
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )

    sim = base.config(sim, timestamps)

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    ph.FluidDiagnostics(
        quantity="density", write_timestamps=timestamps, population_name="protons"
    )

    if LOAD_BALANCE:
        ph.LoadBalancer(active=True, auto=True, mode="nppc", tol=0.05)

    return sim


def plot_file_for_qty(qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty("divb", time),
            plot_patches=True,
            vmin=1e-11,
            vmax=2e-10,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty("Ranks", time),
            plot_patches=True,
        )
        run.GetN(time, pop_name="protons").plot(
            filename=plot_file_for_qty("N", time),
            plot_patches=True,
        )
        for c in ["x", "y", "z"]:
            run.GetB(time).plot(
                filename=plot_file_for_qty(f"b{c}", time),
                qty=f"{c}",
                plot_patches=True,
            )
        run.GetJ(time).plot(
            filename=plot_file_for_qty("jz", time),
            qty="z",
            plot_patches=True,
            vmin=-2,
            vmax=2,
        )


class HarrisTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(HarrisTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(HarrisTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()
        if cpp.mpi_rank() == 0:
            plot(diag_dir)
        if SCOPE_TIMING:
            m_plotting.plot_run_timer_data(diag_dir, cpp.mpi_rank())
        cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    HarrisTest().test_run().tearDown()
