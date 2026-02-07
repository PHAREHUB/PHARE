#!/usr/bin/env python3
import os

import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare import cpp 
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest

os.environ["PHARE_SCOPE_TIMING"] = "1"  # turn on scope timing

ph.NO_GUI()

final_time = 25.6
time_step = 0.0005
diag_dir = "phare_outputs/orszag_tang"

time_step_nbr = int(final_time / time_step)
start_dump_time = 0.0
dumpfrequency = 1000
dt = dumpfrequency * time_step
timestamps = (
    dt * np.arange(int((final_time - start_dump_time) / dt) + 1) + start_dump_time
)

hall=True
res=False
hyper_res=True

def config():
    cells = (256, 256)
    dl = (0.2, 0.2)

    sim = ph.Simulation(
        smallest_patch_size=15,
        # largest_patch_size=25,
        time_step_nbr=time_step_nbr,
        time_step=time_step,
        cells=cells,
        dl=dl,
        interp_order=2,
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
        nesting_buffer=1,
        hyper_mode="spatial",
        eta=0.0,
        nu=0.02,
        gamma=5.0 / 3.0,
        reconstruction="Linear",
        limiter="MinMod",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
        hall=hall,
        res=res,
        hyper_res=hyper_res,
        model_options=["MHDModel"],
    )

    B0 = 1.0 / (np.sqrt(4.0 * np.pi))

    def density(x, y):
        return 25.0 / (36.0 * np.pi)

    def vx(x, y):
        Ly = sim.simulation_domain()[1]
        return -np.sin(2.0 * np.pi * y/Ly)

    def vy(x, y):
        Lx = sim.simulation_domain()[0]
        return np.sin(2.0 * np.pi * x/Lx)

    def vz(x, y):
        return 0.0

    def bx(x, y):
        Ly = sim.simulation_domain()[1]
        return -B0 * np.sin(2.0 * np.pi * y/Ly)

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        return B0 * np.sin(4.0 * np.pi * x/Lx)

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 5.0 / (12.0 * np.pi)

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/orszag_tang_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty(plot_dir, "divb", time),
            plot_patches=True,
            vmin=-1e-11,
            vmax=1e-11,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty(plot_dir, "Ranks", time), plot_patches=True
        )
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
                filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
        run.GetMHDP(time).plot(
            filename=plot_file_for_qty(plot_dir, "p", time), plot_patches=True
        )


class OrszagTangTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(OrszagTangTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(OrszagTangTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()
        if cpp.mpi_rank() == 0:
            plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
            plot_dir.mkdir(parents=True, exist_ok=True)
            plot(diag_dir, plot_dir)
        cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    startMPI()
    OrszagTangTest().test_run().tearDown()
