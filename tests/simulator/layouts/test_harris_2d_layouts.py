import numpy as np

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.simulator.simulator import startMPI

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags
from tests.simulator.layouts import compare_hierarchies

ph.NO_GUI()
cells = (100, 50)
time_step = 0.005
final_time = time_step * 2
ppc = 55


def config(layout):
    L = 0.5

    diag_dir = f"phare_outputs/test/layouts/harris_2d_{layout}/{cpp.mpi_size()}"
    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=3,
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite", "fine_dump_lvl_max": 10},
        },
        nesting_buffer=1,
        particle_layout=layout,
    )

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        return dBy1 + dBy2

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 0.7
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": ppc,
    }

    model = ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 12334}},
    )
    dump_all_diags(model.populations)
    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.LoadBalancer(active=True, every=1, mode="nppc", tol=0.001)
    return sim, diag_dir


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

    def _run(self, layout):
        ph.global_vars.sim = None
        sim, diag_dir = config(layout)
        # self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(sim).run().reset()
        return diag_dir

    def test_run(self):
        diag_dir0 = self._run("AoSMapped")
        diag_dir1 = self._run("AoSPCTS")

        # if cpp.mpi_rank() == 0:
        #     compare_hierarchies(
        #         self,
        #         Run(diag_dir0),
        #         Run(diag_dir1),
        #         atol=dict(b=5e-15, e=3e-14, moments=2e-14, particles=1e-14),
        #     )
        return self


if __name__ == "__main__":
    startMPI()
    HarrisTest().test_run().tearDown()
