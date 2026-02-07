#

import numpy as np
from pathlib import Path

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.run import Run

from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI


from tests.simulator import SimulatorTest


ph.NO_GUI()

# dir(cpp)

start_time = 0
cells = (150, 240, 50)
dl = (0.4, 0.4, 0.4)
diag_dir = "phare_outputs/harris_3d"
time_step = 0.002
final_time = 50
timestamps = np.arange(
    0.0, final_time + time_step, final_time / 100
)  # [] # fix if you want them - np.arange(start_time, 11, .002)

hs = hour_seconds = 3600.0
elapsed_restart_timestamps = [ hs * 1, hs * 2, hs * 4, hs * 9, hs * 12, hs * 15, hs * 23]

hall = True
res = True
hyper_res = True


def config():
    L = 0.5

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        dl=dl,
        cells=cells,
        interp_order=2,
        refinement="tagging",
        max_nbr_levels=1,
        nesting_buffer=1,
        # tagging_threshold=0.5,
        # hyper_resistivity=0.008,
        hyper_mode="spatial",
        resistivity=0.001,
        diag_options={
            "format": "pharevtkhdf",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        restart_options={
            "dir": "checkpoints",
            "mode": "overwrite",
            "elapsed_timestamps": elapsed_restart_timestamps,
            # "restart_time": start_time,
        },
        write_reports=False,
        strict=True,
        eta=0.0,
        nu=0.02,
        gamma=5.0 / 3.0,
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        mhd_timestepper="TVDRK3",
        hall=hall,
        res=res,
        hyper_res=hyper_res,
        model_options=["MHDModel"],
    )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y, z):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1
        x0 = x - 0.5 * Lx
        y1 = y - 0.25 * Ly
        y2 = y - 0.75 * Ly
        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)
        return dBy1 + dBy2

    def bx(x, y, z):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1
        x0 = x - 0.5 * Lx
        y1 = y - 0.25 * Ly
        y2 = y - 0.75 * Ly
        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)
        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.25, L) - S(y, Ly * 0.75, L)) + dBx1 + dBx2

    def bz(x, y, z):
        return 0.0

    def b2(x, y, z):
        return bx(x, y, z) ** 2 + by(x, y, z) ** 2 + bz(x, y, z) ** 2

    # maybe use T to compute this
    def p(x, y, z):
        return 1.0 - b2(x, y, z) / 2.0

    def density(x, y, z):
        return p(x, y, z) / 2.0

    def T(x, y, z):
        K = 0.7
        temp = 1.0 / density(x, y, z) * (K - b2(x, y, z) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y, z):
        return 0.0

    def vy(x, y, z):
        return 0.0

    def vz(x, y, z):
        return 0.0

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    pop_name = "protons"
    for time in timestamps:
        # run.GetDivB(time).plot(
        #     filename=plot_file_for_qty(plot_dir, "divb", time),
        #     plot_patches=True,
        #     vmin=1e-11,
        #     vmax=2e-10,
        # )
        # run.GetRanks(time).plot(
        #     filename=plot_file_for_qty(plot_dir, "Ranks", time), plot_patches=True
        # )
        run.GetN(time, pop_name=pop_name).plot(
            filename=plot_file_for_qty(plot_dir, "N", time), plot_patches=True
        )

        for c in ["x", "y", "z"]:
            run.GetB(time, all_primal=False).plot(
                filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                plot_patches=True,
                qty=f"B{c}",
            )

            run.GetE(time, all_primal=False).plot(
                filename=plot_file_for_qty(plot_dir, f"e{c}", time),
                plot_patches=True,
                qty=f"E{c}",
            )
        # run.GetJ(time).plot(
        #     filename=plot_file_for_qty(plot_dir, "jz", time),
        #     qty="z",
        #     plot_patches=True,
        #     vmin=-2,
        #     vmax=2,
        # )
        # run.GetPressure(time, pop_name=pop_name).plot(
        #     filename=plot_file_for_qty(plot_dir, "Pxx", time),
        #     qty=pop_name + "_Pxx",
        #     plot_patches=True,
        #     vmin=0,
        #     vmax=2.7,
        # )
        # run.GetPressure(time, pop_name=pop_name).plot(
        #     filename=plot_file_for_qty(plot_dir, "Pzz", time),
        #     qty=pop_name + "_Pzz",
        #     plot_patches=True,
        #     vmin=0,
        #     vmax=1.5,
        # )


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
        # self.register_diag_dir_for_cleanup(diag_dir)

        Simulator(config()).run(monitoring=0).reset()

        # if cpp.mpi_rank() == 0:
        #     plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
        #     plot_dir.mkdir(parents=True, exist_ok=True)
        #     plot(diag_dir, plot_dir)
        # cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    startMPI()
    HarrisTest().test_run().tearDown()


