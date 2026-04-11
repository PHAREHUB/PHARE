#
#
#


import sys
import numpy as np
from pathlib import Path

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy import hierarchy_utils as hootils

from tests.simulator import SimulatorTest

ph.NO_GUI()


time_step = 0.005
final_time = 100.0
dt = time_step * 100
timestamps = np.arange(0, final_time + dt, dt)


def density(x):
    return 1.0


def bx(x):
    return 1.0


def by(x):
    return 0.0


def bz(x):
    return 0.0


def T(x):
    return 0.125**2


def vWeak(x):
    from pyphare.pharein.global_vars import sim

    L = sim.simulation_domain()[0]
    x0 = 0.5 * L
    sigma = 2.0
    bubble = 0.08 * np.exp(-((x - x0) ** 2) / (2 * sigma**2))
    return bubble


def vNull(x):
    return 0.0


def vth(x):
    return np.sqrt(T(x))


vvv = {
    "vbulkx": vWeak,
    "vbulky": vNull,
    "vbulkz": vNull,
    "vthx": vth,
    "vthy": vth,
    "vthz": vth,
}


def config(**kwargs):
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        boundary_types="periodic",
        hyper_resistivity=0.001,
        cells=512,
        dl=0.25,
        diag_options={
            "format": "phareh5",
            "options": {"dir": kwargs["diagdir"], "mode": "overwrite"},
        },
    )

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, "nbr_part_per_cell": 200, **vvv},
    )

    ph.ElectronModel(closure="isothermal", Te=kwargs.get("Te", 0.0))

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["domain"]:
        ph.ParticleDiagnostics(
            quantity=quantity, write_timestamps=timestamps, population_name="protons"
        )
    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(test, diag_dir, Te):
    run = Run(diag_dir)
    plot_dir = Path(f"{diag_dir}_plots")
    plot_dir.mkdir(parents=True, exist_ok=True)

    time = 0
    V_finest = run.GetVi(time)[:]
    V_gaussian = V_finest.gaussian()
    for c in ["x", "y", "z"]:
        V_gaussian.plot(
            filename=plot_file_for_qty(plot_dir, f"v{c}_gaussian", time), qty=f"{c}"
        )

    times = np.asarray((0, 20, 40, 60, 80, 100))
    hootils.plot_velocity_peaks_over_time(
        run, times, Te, plot_dir / "velocity_peaks.png", sigma=6
    )


class TestWeakPerbabtion1D(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(TestWeakPerbabtion1D, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(TestWeakPerbabtion1D, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def run(self, Te):
        diag_dir = f"phare_outputs/test_weak_perturbation1d/{cpp.mpi_size()}"
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config(diagdir=diag_dir, Te=Te)).run().reset()

        if cpp.mpi_rank() == 0:
            plot(self, diag_dir, Te)
        return self


if __name__ == "__main__":
    startMPI()

    Te = 0.2
    if len(sys.argv) > 1:
        Te = float(sys.argv[1])

    TestWeakPerbabtion1D().run(Te).tearDown()
