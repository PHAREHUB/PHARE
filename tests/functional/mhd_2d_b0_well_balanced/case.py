#!/usr/bin/env python3
"""
Well-balanced test for the B = B0 + B1 magnetic-field splitting in the MHD solver.

The exact solution is REST: uniform rho, uniform P, V = 0, B1 = 0, and a static
curl-free / divergence-free background field B0. With j0 = curl(B0) = 0 and B1 = 0
there is no force (j x B = 0) and no EMF, so a well-balanced scheme must keep the
state at rest forever.

Any growth of |V| is therefore attributable to the way B0 is handled in the flux /
CT EMF (B0 reconstructed or interpolated -> spurious jump, non-zero discrete
divergence of the B0 Maxwell stress). With B0 sampled analytically everywhere and
its self-stress removed, max|V| stays at machine zero.

Two B0 modes, selected by env var PHARE_B0_MODE:
  - "dipole"  (default): smooth 2D dipole with its singularity OUTSIDE the domain,
                         so B0 is smooth but has a real gradient -> probes well-balancing.
  - "uniform": constant B0 vector (grad B0 = 0) -> control, |V| ~ machine zero trivially.

Env knobs: PHARE_NCELLS (default 128), PHARE_RECON (default "WENOZ"),
           PHARE_LIMITER (default "None"), PHARE_B1_MODE ("zero"|"uniform"),
           PHARE_V0 (uniform bulk Vx, default 0).
"""

import os
import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

small = 1e-12

B0_MODE = os.environ.get("PHARE_B0_MODE", "dipole")
B1_MODE = os.environ.get("PHARE_B1_MODE", "zero")  # "zero" or "uniform"
NCELLS = int(os.environ.get("PHARE_NCELLS", "128"))
RECON = os.environ.get("PHARE_RECON", "WENOZ")
LIMITER = os.environ.get("PHARE_LIMITER", "None")
V0 = float(os.environ.get("PHARE_V0", "0.0"))

gamma = 5.0 / 3.0

L = 4.0                       # domain side length
dipole_center = (-1.5, -1.5)  # outside the domain
m_dipole = 1.0
B0_UNIFORM = (0.0, 0.3, 0.0)
B1_UNIFORM = (0.0, 0.2, 0.0)


def r(x, y):
    return np.sqrt(x**2 + y**2)


def dipole_2d_x(x, y):
    xr = x - dipole_center[0]
    yr = y - dipole_center[1]
    return (m_dipole / (2.0 * np.pi)) * (2.0 * xr * yr / (r(xr, yr) ** 4 + small))


def dipole_2d_y(x, y):
    xr = x - dipole_center[0]
    yr = y - dipole_center[1]
    return (m_dipole / (2.0 * np.pi)) * ((yr**2 - xr**2) / (r(xr, yr) ** 4 + small))


if B0_MODE == "dipole":

    def b0x(x, y):
        return dipole_2d_x(x, y)

    def b0y(x, y):
        return dipole_2d_y(x, y)

elif B0_MODE == "uniform":

    def b0x(x, y):
        return B0_UNIFORM[0] + 0.0 * x

    def b0y(x, y):
        return B0_UNIFORM[1] + 0.0 * x

else:
    raise ValueError(f"unknown PHARE_B0_MODE={B0_MODE!r}")


def b0z(x, y):
    return 0.0 * x


def density(x, y):
    return 1.0 + 0.0 * x


def p(x, y):
    return 1.0 + 0.0 * x


def vzero(x, y):
    return 0.0 * x


def vx0(x, y):
    return V0 + 0.0 * x


cells = np.array([NCELLS, NCELLS])
domain_size = np.array([L, L])
dl = domain_size / cells

B0_scale = 2.0
cs = np.sqrt(gamma * 1.0 / 1.0)
cfast = cs + B0_scale / np.sqrt(1.0)
cfl = 0.2
time_step = cfl * dl.min() / cfast
n_steps = 100
final_time = n_steps * time_step
dump_every = 5
timestamps = np.arange(0, n_steps + 1, dump_every) * time_step

diag_dir = f"phare_outputs_{B0_MODE}_b1{B1_MODE}_n{NCELLS}_{RECON}"


def config():
    sim = ph.Simulation(
        smallest_patch_size=20,
        time_step=time_step,
        final_time=final_time,
        cells=cells.tolist(),
        dl=dl.tolist(),
        interp_order=2,
        refinement="tagging",
        max_mhd_level=1,
        max_nbr_levels=1,
        hall=False,
        res=False,
        hyper_res=False,
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
        reconstruction=RECON,
        limiter=LIMITER,
        riemann="Rusanov",
        mhd_timestepper="TVDRK3",
        model_options=["MHDModel"],
        boundary_types=("physical", "physical"),
        boundary_conditions={
            "xlower": {"type": "open"},
            "xupper": {"type": "open"},
            "ylower": {"type": "open"},
            "yupper": {"type": "open"},
        },
    )

    model_kwargs = dict(
        density=density,
        vx=vx0,
        vy=vzero,
        vz=vzero,
        b0x=b0x,
        b0y=b0y,
        b0z=b0z,
        p=p,
    )
    if B1_MODE == "uniform":
        # uniform perturbation field; total B = B0 + B1. Equilibrium (j1 = curl B1 = 0)
        # so the magnetized fluid stays at rest -> tests cross-term well-balancing.
        model_kwargs.update(
            b1x=lambda x, y: B1_UNIFORM[0] + 0.0 * x,
            b1y=lambda x, y: B1_UNIFORM[1] + 0.0 * x,
            b1z=lambda x, y: B1_UNIFORM[2] + 0.0 * x,
        )
    ph.MHDModel(**model_kwargs)

    ph.MHDDiagnostics(quantity="rho", write_timestamps=timestamps)
    ph.MHDDiagnostics(quantity="P", write_timestamps=timestamps)
    ph.MHDDiagnostics(quantity="V", write_timestamps=timestamps)
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    return sim


def main():
    print(
        f"[b0_well_balanced] mode={B0_MODE} ncells={NCELLS} recon={RECON} "
        f"limiter={LIMITER} dt={time_step:.3e} nsteps={n_steps} dir={diag_dir}"
    )
    Simulator(config()).run().reset()
    ph.global_vars.sim = None


if __name__ == "__main__":
    startMPI()
    main()
