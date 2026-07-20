#!/usr/bin/env python3
"""Inner-boundary functional smoke test: sphere in a fully periodic box.

A uniform flow is initialized around a spherical inner boundary at the centre
of a periodic domain. The run exercises the full inner-boundary stack
(classification, safe-state pinning, ghost-cell moment BCs, degraded fluxes on
under-resolved levels) without requiring outer physical boundary conditions.

Asserts after the run:
  - the run completes without NaN/exception,
  - rho and P are finite and positive everywhere on the dumped hierarchy,
  - cells deep inside the body hold the configured safe-state density.
"""

import os

import numpy as np

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest

ph.NO_GUI()


cells = (64, 64)
dl = (0.15625, 0.15625)  # 10 x 10 box
domain_size = (cells[0] * dl[0], cells[1] * dl[1])
center = (domain_size[0] / 2, domain_size[1] / 2)  # (5, 5)
radius = 1.0

safe_density = 2.0
safe_pressure = 0.5

time_step = 1e-3
time_step_nbr = 50
final_time = time_step * time_step_nbr

diag_dir = "phare_outputs/mhd_inner_boundary_sphere_periodic"


def config():
    sim = ph.Simulation(
        smallest_patch_size=10,
        time_step=time_step,
        time_step_nbr=time_step_nbr,
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
        reconstruction="Linear",
        limiter="VanLeer",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
        hall=False,
        res=False,
        hyper_res=False,
        model_options=["MHDModel"],
        inner_boundary={
            "name": "sphere",
            "shape": "sphere",
            "center": list(center),
            "radius": radius,
            "inactive_safe_state": {
                "density": safe_density,
                "pressure": safe_pressure,
                "velocity": [0.0, 0.0, 0.0],
                "B": 0.0,
            },
        },
    )

    def density(x, y):
        return 1.0 + 0.0 * x

    def vx(x, y):
        return 0.5 + 0.0 * x

    def vy(x, y):
        return 0.0 * x

    def vz(x, y):
        return 0.0 * x

    def bx(x, y):
        return 0.0 * x

    def by(x, y):
        return 0.0 * x

    def bz(x, y):
        return 1e-2 + 0.0 * x

    def p(x, y):
        return 1.0 + 0.0 * x

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    for quantity in ["rho", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=[0.0, final_time])

    return sim


class MHDInnerBoundarySphereTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(MHDInnerBoundarySphereTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(MHDInnerBoundarySphereTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()

        if cpp.mpi_rank() == 0:
            self._check_hierarchy()
        cpp.mpi_barrier()
        return self

    def _check_hierarchy(self):
        import h5py  # lazy import (convention)

        rho_h5 = os.path.join(diag_dir, "mhd_rho.h5")
        p_h5 = os.path.join(diag_dir, "mhd_P.h5")

        rho_hier = hierarchy_from(h5_filename=rho_h5, times=final_time)
        p_hier = hierarchy_from(h5_filename=p_h5, times=final_time)

        checked_fluid = 0
        checked_body = 0
        for ilvl, level in rho_hier.levels(final_time).items():
            for patch in level.patches:
                pd = patch.patch_datas["mhdRho"]
                nbr_ghosts = pd.ghosts_nbr
                data = pd.dataset[:]
                inner = data[tuple(slice(g, -g) for g in nbr_ghosts)]

                self.assertTrue(
                    np.all(np.isfinite(inner)), f"NaN/inf rho on L{ilvl} {patch.id}"
                )
                self.assertTrue(np.all(inner > 0.0), f"rho <= 0 on L{ilvl} {patch.id}")

                # classify cell centres against the sphere
                x = pd.x[nbr_ghosts[0] : -nbr_ghosts[0]]
                y = pd.y[nbr_ghosts[1] : -nbr_ghosts[1]]
                X, Y = np.meshgrid(x, y, indexing="ij")
                r = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

                # deep-inside cells (a couple of cells clear of the ghost shell)
                # must hold the configured safe-state density
                deep = r < radius - 3 * max(dl)
                if np.any(deep):
                    checked_body += deep.sum()
                    np.testing.assert_allclose(
                        inner[deep],
                        safe_density,
                        rtol=1e-12,
                        err_msg=f"safe-state density not held on L{ilvl} {patch.id}",
                    )

                fluid = r > radius + 3 * max(dl)
                checked_fluid += fluid.sum()

        for ilvl, level in p_hier.levels(final_time).items():
            for patch in level.patches:
                pd = patch.patch_datas["mhdP"]
                nbr_ghosts = pd.ghosts_nbr
                inner = pd.dataset[:][tuple(slice(g, -g) for g in nbr_ghosts)]
                self.assertTrue(
                    np.all(np.isfinite(inner)), f"NaN/inf P on L{ilvl} {patch.id}"
                )
                self.assertTrue(np.all(inner > 0.0), f"P <= 0 on L{ilvl} {patch.id}")

        self.assertGreater(checked_body, 0, "no in-body cells were checked")
        self.assertGreater(checked_fluid, 0, "no fluid cells were checked")


def main():
    startMPI()
    MHDInnerBoundarySphereTest().test_run().tearDown()


if __name__ == "__main__":
    main()
