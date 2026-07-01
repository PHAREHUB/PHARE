#

import numpy as np
from pathlib import Path

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI

from tests.simulator import SimulatorTest


ph.NO_GUI()


start_time = 0
cells = (105, 105, 105)
# cells = (22, 22, 22)
dl = (0.4, 0.4, 0.4)
diag_dir = "phare_outputs/harris_3d"
time_step = 0.001
final_time = 0.001
timestamps = []

hs = hour_seconds = 3600.0
elapsed_restart_timestamps = [hs * 1, hs * 3, hs * 6]
ppc = 3


def config():
    L = 0.5

    sim = ph.Simulation(
        interp_order=1,
        time_step=time_step,
        final_time=final_time,
        dl=dl,
        cells=cells,
        # refinement="tagging",
        # max_nbr_levels=1,
        refinement_boxes={"L0": {"B0": [[5] * 3, [99] * 3]}},
        # refinement_boxes={"L0": {"B0": [[3] * 3, [8] * 3]}},
        nesting_buffer=1,
        tagging_threshold=0.5,
        hyper_resistivity=0.008,
        hyper_mode="spatial",
        resistivity=0.001,
        # diag_options={
        #     "format": "pharevtkhdf",
        #     "options": {"dir": diag_dir, "mode": "overwrite"},
        # },
        restart_options={
            # "dir": "checkpoints",
            "mode": "overwrite",
            # "elapsed_timestamps": [0],
            "restart_time": "auto",
        },
        write_reports=False,
        strict=False,
        tag_buffer=3,
    )

    def density(x, y, z):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def by(x, y, z):
        return 0

    def bx(x, y, z):
        return 0

    def bz(x, y, z):
        return 0.0

    def b2(x, y, z):
        return bx(x, y, z) ** 2 + by(x, y, z) ** 2 + bz(x, y, z) ** 2

    def T(x, y, z):
        K = 0.7
        temp = 1.0 / density(x, y, z) * (K - b2(x, y, z) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y, z):
        return 1e-15

    def vy(x, y, z):
        return 1e-15

    def vz(x, y, z):
        return 1e-15

    def vthx(x, y, z):
        return 1e-15

    def vthy(x, y, z):
        return 1e-15

    def vthz(x, y, z):
        return 1e-15

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": ppc,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
            "init": {"seed": cpp.mpi_rank() + 12},
        },
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.LoadBalancer(active=True, mode="nppc", tol=0.05, every=1000)

    return sim


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
        sim = Simulator(config()).run(monitoring=0)
        sim.print_summary()
        sim.reset()
        return self


if ph.PHARE_EXE:
    config()

elif __name__ == "__main__":
    startMPI()
    HarrisTest().test_run().tearDown()
