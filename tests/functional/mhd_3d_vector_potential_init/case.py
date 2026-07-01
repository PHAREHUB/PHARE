#!/usr/bin/env python3
"""
Vector-potential init for the MHD B = B0 + B1 splitting (3D).

3D analogue of mhd_2d_vector_potential_init: exercises the FULL discrete curl B = curl(A) with a
3D vector potential A = (Ax, Ay, Az) (the 2D test only uses the out-of-plane Az). The potentials
are chosen so B0 and B1 are UNIFORM (hence force-free) and the resulting state stays at rest:

  a0 = (0, -B0x * z, 0)  ->  B0 = curl(a0) = (B0x, 0, 0)   via  Bx = -dAy/dz   (tests deriv<Z>)
  a1 = (0,  B1z * x, 0)  ->  B1 = curl(a1) = (0, 0, B1z)   via  Bz =  dAy/dx   (tests deriv<X>)

Both fields are divergence-free analytically AND discretely (div . curl = 0). The initial state
is REST (V = 0); uniform B exerts no Lorentz force so |V| must stay ~0 over the short run.

Selected by env var PHARE_B_INIT:
  - "potential" (default): B0 from a0*, B1 from a1*
  - "components":          B0 from b0x, B1 from b1z

Knobs: PHARE_NCELLS (default 16), PHARE_RECON (WENOZ), PHARE_LIMITER (None).
"""

import os
import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()

B_INIT = os.environ.get("PHARE_B_INIT", "potential")
NCELLS = int(os.environ.get("PHARE_NCELLS", "16"))
RECON = os.environ.get("PHARE_RECON", "WENOZ")
LIMITER = os.environ.get("PHARE_LIMITER", "None")

gamma = 5.0 / 3.0

L = 1.0
B0x = 0.3  # uniform background field along +x
B1z = 0.2  # uniform perturbation field along +z

density = lambda x, y, z: 1.0 + 0.0 * x
p = lambda x, y, z: 1.0 + 0.0 * x
vzero = lambda x, y, z: 0.0 * x
zero = lambda x, y, z: 0.0 * x

# --- vector potentials (full 3D): only the y-component is nonzero ---
# a0 = (0, -B0x z, 0) -> B0 = (B0x, 0, 0)
a0y = lambda x, y, z: -B0x * z
# a1 = (0, B1z x, 0) -> B1 = (0, 0, B1z)
a1y = lambda x, y, z: B1z * x

# --- component-wise counterparts (uniform) ---
b0x = lambda x, y, z: B0x + 0.0 * x
b1z = lambda x, y, z: B1z + 0.0 * x

cells = np.array([NCELLS, NCELLS, NCELLS])
dl = np.array([L / NCELLS, L / NCELLS, L / NCELLS])
time_step = 0.2 * dl.min()
n_steps = 4
final_time = n_steps * time_step
timestamps = np.arange(0, n_steps + 1, 2) * time_step

diag_dir = f"phare_outputs_{B_INIT}_n{NCELLS}_{RECON}"


def config():
    sim = ph.Simulation(
        smallest_patch_size=8,
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
        # SSPRK4_5 (not TVDRK3): the built ideal (no Hall/res/hyper) 3D WENOZ permutation.
        mhd_timestepper="SSPRK4_5",
        model_options=["MHDModel"],
    )

    common = dict(density=density, vx=vzero, vy=vzero, vz=vzero, p=p)
    if B_INIT == "potential":
        ph.MHDModel(
            a0x=zero, a0y=a0y, a0z=zero, a1x=zero, a1y=a1y, a1z=zero, **common
        )
    elif B_INIT == "components":
        ph.MHDModel(b0x=b0x, b1z=b1z, **common)
    else:
        raise ValueError(f"unknown PHARE_B_INIT={B_INIT}")

    ph.MHDDiagnostics(quantity="rho", write_timestamps=timestamps)
    ph.MHDDiagnostics(quantity="V", write_timestamps=timestamps)

    return sim


def main():
    print(f"[vector_potential_init_3d] B_INIT={B_INIT} ncells={NCELLS} dir={diag_dir}")
    Simulator(config()).run().reset()
    ph.global_vars.sim = None


if __name__ == "__main__":
    startMPI()
    main()
