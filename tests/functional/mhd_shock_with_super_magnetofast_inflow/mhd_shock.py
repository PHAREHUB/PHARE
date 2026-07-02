#!/usr/bin/env python3

import os
import numpy as np
from pathlib import Path
from dataclasses import field

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI
from dataclasses import dataclass

from tests.simulator import SimulatorTest

os.environ["PHARE_SCOPE_TIMING"] = "1"  # turn on scope timing

ph.NO_GUI()

@dataclass
class State:
    gamma: float
    rho: float
    p: float
    bx: float
    by: float
    b: float = field(init=False) 
    cs: float = field(init=False)
    va: float = field(init=False)
    cf_fast: float = field(init=False)
    cf_slow: float = field(init=False)

    def __post_init__(self):
        b = np.linalg.norm([self.bx, self.by])
        cs = np.sqrt(self.gamma * self.p / self.rho)
        va = b / np.sqrt(self.rho)
        cos_theta = self.bx / b
        common = np.sqrt((cs**2 + va**2)**2 - (4.0 * cs**2 * va**2 * cos_theta**2))
        cf_fast = np.sqrt((cs**2 + va ** 2 + common) / 2)
        cf_slow = np.sqrt((cs**2 + va ** 2 - common) / 2)
        self.b, self.cs, self.va, self.cf_fast, self.cf_slow = b, cs, va, cf_fast, cf_slow


gamma = 5.0 / 3.0
BX = 0.75
LEFT_INIT = State(gamma, 1.0, 1.0, BX, 1.0)
RIGHT_INIT = State(gamma, 0.125, 0.1, BX, -1.0)
cf_max = max(LEFT_INIT.cf_fast, RIGHT_INIT.cf_fast)
U_INIT = 1.05 * cf_max
ncells = 800
dx = 1.0
cfl = 0.5
time_step = cfl * dx / (U_INIT + cf_max) 
nsafe_cells = 10

nsteps = 1000
final_time = nsteps * time_step
dump_freq = 1
# time_step = 0.2
timestamps = np.arange(0.0, final_time, dump_freq*time_step)
diag_dir = "phare_outputs/shock"


def print_case_info():
    things_to_print = {
        "left state": LEFT_INIT,
        "right state": RIGHT_INIT,
        "cf_max": cf_max,
        "dt": time_step
    }
    if cpp.mpi_rank() == 0:
        for key, thing in things_to_print.items():
            print(f"{key} = {thing}")


def config():
    cells = (ncells,)
    dl = (dx,)

    sim = ph.Simulation(
        smallest_patch_size=15,
        # largest_patch_size=25,
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_mhd_level=1,
        max_nbr_levels=1,
        hyper_resistivity=0.0,
        resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        interp_order=2,
        mhd_timestepper="TVDRK3",
        model_options=["MHDModel"],
        boundary_types="physical",
        boundary_conditions={
            "xlower": {
                "type": "super-magnetofast-inflow",
                "data": {
                    "velocity": U_INIT,
                    "density": LEFT_INIT.rho,
                    "pressure": LEFT_INIT.p,
                    "B" : [BX, LEFT_INIT.by, 0.0],
                },
            },
            "xupper": {"type": "super-magnetofast-outflow"},
        },
    )

    def density(x):
        return np.where(x < nsafe_cells * dx, LEFT_INIT.rho, RIGHT_INIT.rho)

    def vx(x):
        return U_INIT

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def bx(x):
        return BX

    def by(x):
        return np.where(x < nsafe_cells * dx, LEFT_INIT.by, RIGHT_INIT.by)

    def bz(x):
        return 0.0

    def p(x):
        return np.where(x < nsafe_cells * dx, LEFT_INIT.p, RIGHT_INIT.p)

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/shock_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetMHDrho(time).plot(
            filename=plot_file_for_qty(plot_dir, "rho", time), plot_patches=True
        )
        for c in ["x", "y"]:
            run.GetMHDV(time).plot(
                filename=plot_file_for_qty(plot_dir, f"v{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
        run.GetB(time).plot(
            filename=plot_file_for_qty(plot_dir, "by", time),
            plot_patches=True,
            qty="y",
        )
        run.GetMHDP(time).plot(
            filename=plot_file_for_qty(plot_dir, "p", time), plot_patches=True
        )


class ShockTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(ShockTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(ShockTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        # self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()
        print_case_info()
        # if cpp.mpi_rank() == 0:
        #     plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
        #     plot_dir.mkdir(parents=True, exist_ok=True)
        #     plot(diag_dir, plot_dir)
        cpp.mpi_barrier()
        return self


def main():
    Simulator(config()).run()


if __name__ == "__main__":
    startMPI()
    ShockTest().test_run().tearDown()
