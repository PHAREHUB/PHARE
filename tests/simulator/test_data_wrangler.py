#!/usr/bin/env python3

import unittest
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator
from tests.simulator import populate_simulation

from pyphare.pharesee.hierarchy.fromsim import hierarchy_from_sim


def populate_simulation_mhd(dim, interp):
    from pyphare.pharein.simulation import check_patch_size

    cells = [20 for _ in range(dim)]
    _, smallest_patch_size = check_patch_size(dim, interp_order=interp, cells=cells)
    dl = [1.0 / v for v in cells]

    ph.global_vars.sim = None
    sim = ph.Simulation(
        interp_order=interp,
        smallest_patch_size=smallest_patch_size,
        largest_patch_size=[20] * dim,
        time_step_nbr=1,
        final_time=0.001,
        boundary_types=["periodic"] * dim,
        cells=cells,
        dl=dl,
        diag_options={},
        nesting_buffer=0,
        strict=True,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="Linear",
        limiter="VanLeer",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
        hall=False,
        model_options=["MHDModel"],
        max_mhd_level=1,
    )

    ph.MHDModel(
        density=lambda *xyz: 1.0,
        bx=lambda *xyz: 1.0,
        by=lambda *xyz: 0.0,
        bz=lambda *xyz: 0.0,
        vx=lambda *xyz: 0.0,
        vy=lambda *xyz: 0.0,
        vz=lambda *xyz: 0.0,
        p=lambda *xyz: 1.0,
    )

    return sim


class DataWranglerTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(DataWranglerTest, self).__init__(*args, **kwargs)
        self.dw = None
        self.simulator = None

    def test_hybrid(self):
        ndim, interp = 2, 1
        self.simulator = Simulator(populate_simulation(ndim, interp))
        self.simulator.initialize()

        hier = None
        hier = hierarchy_from_sim(self.simulator, qty="EM_B_x", hier=hier, sync=True)
        hier = hierarchy_from_sim(self.simulator, qty="EM_B_y", hier=hier, sync=True)
        hier = hierarchy_from_sim(self.simulator, qty="density", hier=hier, sync=True)

        if cpp.mpi_rank() == 0:
            hier.plot(
                filename="data_wrangler.Bx.png",
                qty="EM_B_x",
                plot_patches=True,
                levels=(0,),
            )
            hier.plot(
                filename="data_wrangler.By.png",
                qty="EM_B_y",
                plot_patches=True,
                levels=(0,),
            )

        hier = hierarchy_from_sim(self.simulator, qty="particles", pop="protons")
        ppc = 100  # set in makeBasicModel / defaultPopulationSettings
        expected = ppc * int(np.prod([20] * ndim))
        total = sum(p.patch_datas["particles"].size() for p in hier.level(0).patches)
        self.assertGreaterEqual(total, (expected * 0.95) / cpp.mpi_size())

    def test_mhd(self):
        ndim, interp = 2, 1
        self.simulator = Simulator(populate_simulation_mhd(ndim, interp))
        self.simulator.initialize()

        hier = None
        hier = hierarchy_from_sim(self.simulator, qty="rho", hier=hier)
        hier = hierarchy_from_sim(self.simulator, qty="EM_B_x", hier=hier)
        hier = hierarchy_from_sim(self.simulator, qty="V_x", hier=hier)
        hier = hierarchy_from_sim(self.simulator, qty="P", hier=hier)
        hier = hierarchy_from_sim(self.simulator, qty="Etot", hier=hier)

        if cpp.mpi_rank() == 0:
            hier.plot(
                filename="data_wrangler_mhd.rho.png", qty="rho", plot_patches=True
            )
            hier.plot(
                filename="data_wrangler_mhd.Bx.png", qty="EM_B_x", plot_patches=True
            )
            hier.plot(filename="data_wrangler_mhd.Vx.png", qty="V_x", plot_patches=True)
            hier.plot(filename="data_wrangler_mhd.P.png", qty="P", plot_patches=True)
            hier.plot(
                filename="data_wrangler_mhd.Etot.png", qty="Etot", plot_patches=True
            )

    def tearDown(self):
        del self.dw
        if self.simulator is not None:
            self.simulator.reset()


if __name__ == "__main__":
    unittest.main()
