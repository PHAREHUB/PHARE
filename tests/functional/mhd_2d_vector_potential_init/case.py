#!/usr/bin/env python3
"""
Vector-potential init for the MHD B = B0 + B1 splitting (2D).

Both B0 (a dipole) and B1 (a uniform IMF) are divergence-free analytically. When prescribed
component-wise (b0x/b0y, b1x/b1y) and sampled at the staggered face centres, the *discrete*
divergence div B = dx Bx + dy By is generically nonzero. When prescribed through an out-of-plane
vector potential A_z (a0z / a1z), B = curl(A_z z_hat) is built with the same discrete curl as
Faraday, so the discrete div B = 0 to machine precision (div . curl = 0).

The initial state is REST (V = 0). Selected by env var PHARE_B_INIT:
  - "potential" (default): B0 from a0z, B1 from a1z          -> discrete div B ~ machine zero
  - "components":          B0 from b0x/b0y, B1 from b1x/b1y  -> discrete div B != 0

Knobs: PHARE_NCELLS (default 64), PHARE_RECON (WENOZ), PHARE_LIMITER (None).
"""

import os
import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

small = 1e-12

B_INIT = os.environ.get("PHARE_B_INIT", "potential")
NCELLS = int(os.environ.get("PHARE_NCELLS", "64"))
RECON = os.environ.get("PHARE_RECON", "WENOZ")
LIMITER = os.environ.get("PHARE_LIMITER", "None")

gamma = 5.0 / 3.0

# --- geometry: dipole singularity OUTSIDE the [0,L]^2 domain (smooth B0 inside) ---
L = 4.0
dipole_center = (-1.5, -1.5)
m_dipole = 1.0
b1_imf = 0.2  # uniform perturbation field along -y (IMF in B1)

density = lambda x, y: 1.0 + 0.0 * x
p = lambda x, y: 1.0 + 0.0 * x
vzero = lambda x, y: 0.0 * x


def r(x, y):
    return np.sqrt(x**2 + y**2)


def _xr(x, y):
    return x - dipole_center[0]


def _yr(x, y):
    return y - dipole_center[1]


# Analytic dipole B0 = curl(psi0 z_hat), psi0 = (m/2pi) xr / r^2 (xr,yr relative to singularity)
def dipole_b0x(x, y):
    xr, yr = _xr(x, y), _yr(x, y)
    return -(m_dipole / (2.0 * np.pi)) * (2.0 * xr * yr / (r(xr, yr) ** 4 + small))


def dipole_b0y(x, y):
    xr, yr = _xr(x, y), _yr(x, y)
    return -(m_dipole / (2.0 * np.pi)) * ((yr**2 - xr**2) / (r(xr, yr) ** 4 + small))


def psi0(x, y):
    # B0 = curl(psi0 z_hat): Bx = d psi0/dy = dipole_b0x, By = -d psi0/dx = dipole_b0y
    xr, yr = _xr(x, y), _yr(x, y)
    return (m_dipole / (2.0 * np.pi)) * xr / (r(xr, yr) ** 2 + small)


# Uniform IMF B1 = (0, -b1_imf, 0) = curl(a1z z_hat) with a1z = b1_imf * x
def imf_b1x(x, y):
    return 0.0 * x


def imf_b1y(x, y):
    return -b1_imf + 0.0 * x


def a1z(x, y):
    # Bx = d a1z/dy = 0 ; By = -d a1z/dx = -b1_imf  -> a1z = b1_imf * x
    return b1_imf * x


cells = np.array([NCELLS, NCELLS])
dl = np.array([L / NCELLS, L / NCELLS])
time_step = 0.2 * dl.min()
n_steps = 4
final_time = n_steps * time_step
timestamps = np.arange(0, n_steps + 1, 2) * time_step

diag_dir = f"phare_outputs_{B_INIT}_n{NCELLS}_{RECON}"


def config():
    sim = ph.Simulation(
        smallest_patch_size=20,
        time_step=time_step,
        final_time=final_time,
        cells=cells.tolist(),
        dl=dl.tolist(),
        interp_order=1,
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
    )

    common = dict(density=density, vx=vzero, vy=vzero, vz=vzero, p=p)
    if B_INIT == "potential":
        ph.MHDModel(a0z=psi0, a1z=a1z, **common)
    elif B_INIT == "components":
        ph.MHDModel(b0x=dipole_b0x, b0y=dipole_b0y, b1x=imf_b1x, b1y=imf_b1y, **common)
    else:
        raise ValueError(f"unknown PHARE_B_INIT={B_INIT}")

    ph.MHDDiagnostics(quantity="rho", write_timestamps=timestamps)
    ph.MHDDiagnostics(quantity="V", write_timestamps=timestamps)

    return sim


def main():
    print(f"[vector_potential_init] B_INIT={B_INIT} ncells={NCELLS} dir={diag_dir}")
    Simulator(config()).run().reset()
    ph.global_vars.sim = None


if __name__ == "__main__":
    startMPI()
    main()
