import os
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
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
        common = np.sqrt((cs**2 + va**2) ** 2 - (4.0 * cs**2 * va**2 * cos_theta**2))
        cf_fast = np.sqrt((cs**2 + va**2 + common) / 2)
        cf_slow = np.sqrt((cs**2 + va**2 - common) / 2)
        self.b, self.cs, self.va, self.cf_fast, self.cf_slow = (
            b,
            cs,
            va,
            cf_fast,
            cf_slow,
        )


gamma = 5.0 / 3.0
BX = 0.75
LEFT_INIT = State(gamma, 1.0, 1.0, BX, 1.0)
RIGHT_INIT = State(gamma, 0.125, 0.1, BX, -1.0)
ncells = 800
dx = 1.0


def config(time_step, final_time, timestamps, diag_dir, x0, u, **kwargs):
    sim = ph.Simulation(
        smallest_patch_size=15,
        time_step=time_step,
        final_time=final_time,
        cells=(ncells,),
        dl=(dx,),
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
        gamma=gamma,
        reconstruction="Linear",
        limiter="VanLeer",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
        hall=False,
        model_options=["MHDModel"],
        **kwargs,
    )

    def density(x):
        return np.where(x < x0, LEFT_INIT.rho, RIGHT_INIT.rho)

    def vx(x):
        return u

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def bx(x):
        return BX

    def by(x):
        return np.where(x < x0, LEFT_INIT.by, RIGHT_INIT.by)

    def bz(x):
        return 0.0

    def p(x):
        return np.where(x < x0, LEFT_INIT.p, RIGHT_INIT.p)

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/shock_{qty}_t{time}.png"


def plot(diag_dir, plot_dir, timestamps):
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
        plt.close("all")


class ShockTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super().tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def check(self, diag_dir):
        pass

    def run_case(self, time_step, final_time, timestamps, diag_dir, x0, u, **kwargs):
        self.register_diag_dir_for_cleanup(diag_dir)
        sim = config(time_step, final_time, timestamps, diag_dir, x0, u, **kwargs)
        Simulator(sim).run().reset()
        if not sim.dry_run and cpp.mpi_rank() == 0:
            self.check(diag_dir)
            plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
            plot_dir.mkdir(parents=True, exist_ok=True)
            plot(diag_dir, plot_dir, timestamps)
        cpp.mpi_barrier()
        return self
