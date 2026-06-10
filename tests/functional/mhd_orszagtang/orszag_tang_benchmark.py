#!/usr/bin/env python3
"""
Wall-clock benchmark of the Orszag-Tang MHD case, constant vs adaptive time step.
"""
import os

# no scope timing in a benchmark: it must not be part of what we measure
os.environ["PHARE_SCOPE_TIMING"] = "0"

import time

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest

from orszag_tang import config, final_time

ph.NO_GUI()

CONSTANT_DT = 0.0007
CFL_WAVE = 1.0

time_steps = {
    "constant": {"mode": "constant", "value": CONSTANT_DT},
    "adaptive": {"mode": "adaptive", "cfl_wave": CFL_WAVE},
}


def run_one(time_step):
    """Returns (advance wall clock time in seconds, number of time steps)."""
    ph.global_vars.sim = None
    simulator = Simulator(config(time_step=time_step, with_diags=False))
    simulator.initialize()  # initialization is not what we benchmark

    steps = 0
    advance = simulator.advance

    def counted_advance():
        nonlocal steps
        steps += 1
        return advance()

    simulator.advance = counted_advance

    cpp.mpi_barrier()
    tick = time.perf_counter()
    simulator.run()
    tock = time.perf_counter()
    cpp.mpi_barrier()

    simulator.reset()

    ph.global_vars.sim = None
    return tock - tick, steps


class OrszagTangBenchmark(SimulatorTest):
    def test_benchmark(self):
        results = {name: run_one(ts) for name, ts in time_steps.items()}

        if any(steps == 0 for _, steps in results.values()):
            cpp.mpi_barrier()
            return self  # dry run, nothing was advanced

        if cpp.mpi_rank() == 0:
            constant, adaptive = results["constant"], results["adaptive"]
            print(f"\nOrszag-Tang benchmark - final_time = {final_time}")
            print(f"{'mode':<10}{'steps':>8}{'wall clock (s)':>18}{'s/step':>12}")
            for name, (wall, steps) in results.items():
                print(f"{name:<10}{steps:>8}{wall:>18.3f}{wall / steps:>12.4f}")
            print(f"speedup (constant / adaptive) = {constant[0] / adaptive[0]:.2f}x")

        cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    startMPI()
    OrszagTangBenchmark().test_benchmark().tearDown()
