#!/usr/bin/env python3
"""Visualization case for inner boundary mesh data at initialization.

A sphere is placed at the centre of the domain.  The signed-distance field
and cell-status map are written as VTK diagnostics at t=0 so that they can
be inspected in ParaView or VisIt right after running this script.
"""

import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

cells = (50, 50)
dl = (0.2, 0.2)
domain_size = (cells[0] * dl[0], cells[1] * dl[1])  # 10 x 10
center = (domain_size[0] / 2, domain_size[1] / 2)   # (5, 5)
radius = 2.0

time_step = 0.001
diag_dir = "phare_outputs/mhd_inner_boundary_sphere"


def config():
    sim = ph.Simulation(
        smallest_patch_size=10,
        time_step=time_step,
        final_time=time_step,  # one step — diagnostics are captured at t=0 (init)
        cells=cells,
        dl=dl,
        interp_order=2,
        refinement="tagging",
        max_mhd_level=1,
        max_nbr_levels=1,
        hall=True,
        res=False,
        hyper_res=True,
        hyper_resistivity=0.0,
        resistivity=0.0,
        diag_options={
            "format": "pharevtkhdf",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        mhd_timestepper="TVDRK3",
        model_options=["MHDModel"],
        inner_boundary={
            "name": "sphere",
            "shape": "sphere",
            "center": list(center),
            "radius": radius,
        },
    )

    def density(x, y):
        return 1.0

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def bx(x, y):
        return 1.0

    def by(x, y):
        return 0.0

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 1.0

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    # Inner boundary mesh data — main purpose of this script.
    ph.MHDDiagnostics(quantity="IBSignedDistance", write_timestamps=[0.0])
    ph.MHDDiagnostics(quantity="IBCellStatus", write_timestamps=[0.0])

    # Standard MHD quantities for context.
    ph.MHDDiagnostics(quantity="rho", write_timestamps=[0.0])
    ph.MHDDiagnostics(quantity="P", write_timestamps=[0.0])
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=[0.0])

    return sim


def main():
    Simulator(config()).run().reset()
    ph.global_vars.sim = None


if __name__ == "__main__":
    startMPI()
    main()
