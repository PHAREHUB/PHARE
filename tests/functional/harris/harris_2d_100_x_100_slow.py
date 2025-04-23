#!/usr/bin/env python3

import os
import numpy as np
import matplotlib as mpl
from pathlib import Path

import pyphare.pharein as ph

from pyphare.cpp import cpp_lib
from pyphare.simulator.simulator import startMPI

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


class HarrisTest(base.HarrisTest):
    def __init__(self, *args, **kwargs):
        super(HarrisTest, self).__init__(*args, **kwargs)

    def test_run(self):
        super(HarrisTest, self).test_run(diag_dir, config())
        if cpp.mpi_rank() == 0:
            self.plot(timestamps, diag_dir, plot_dir)

        if SCOPE_TIMING:
            self.scope_timing(diag_dir)
        return self


if __name__ == "__main__":
    HarrisTest().test_run().tearDown()
