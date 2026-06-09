#!/usr/bin/env python3
import os

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import pyphare.pharein as ph
from pyphare import cpp 
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest

# Note: this test does not handle mpi yet (would require gathering several patches on rank 0).
# It scans every reconstruction in one run, each compiled as its own permutation in res/sim/all.txt.

os.environ["PHARE_SCOPE_TIMING"] = "1"

ph.NO_GUI()

mode_speed = {
    "Alfven": 1.0,
    "Fast": 2.0,
    "Slow": 0.5,
    "Entropy": 1.0,
}

mode = "Alfven"
time_step = 0.0007
final_time = 1.0 / mode_speed[mode]
timestamps = [0.0, final_time]
diag_dir = "phare_outputs/convergence"

# The 3D scheme is globally 2nd order (midpoint flux quadrature, face projections),
# so every reconstruction converges at order 2 here regardless of its 1D order.
expected_orders = {
    "Linear": 2.0,
    "WENO3": 2.0,
    "WENOZ": 2.0,
    "MP5": 2.0,
}

# Limiter per reconstruction (limiters are only valid with Linear).
limiters = {
    "Linear": "VanLeer",
    "WENO3": "None",
    "WENOZ": "None",
    "MP5": "None",
}

# SSPRK4_5 so the temporal error never dominates the spatial convergence.
mhd_timestepper = "SSPRK4_5"
ghosts = 4

tolerance = 0.15


def config(nx, reconstruction, limiter):
    sim = ph.Simulation(
        # smallest_patch_size=15,
        # largest_patch_size=25,
        time_step=time_step,
        final_time=final_time,
        cells=(2*nx,nx,nx),
        dl=(3.0 / (2 * nx), 1.5 / nx, 1.5 / nx),        refinement="tagging",
        max_mhd_level=1,
        max_nbr_levels=1,
        hyper_resistivity=0.0,
        resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
        nesting_buffer=1,
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction=reconstruction,
        limiter=limiter,
        riemann="Rusanov",
        mhd_timestepper=mhd_timestepper,
        model_options=["MHDModel"],
    )

    eps = 1.e-6
    background = np.asarray([1., 1. if mode == "Entropy" else 0., 0., 0., 1., 3./2., 0., 1./sim.gamma])
    Rfast_p = 1./(2.*np.sqrt(5.)) * np.asarray([2., 4., 2., 0., 0., 4., 0., 9.])
    Rfast_m = 1./(2.*np.sqrt(5.)) * np.asarray([2., -4., -2., 0., 0., 4., 0., 9.])
    Ralven_p = np.asarray([0., 0., 0., 1., 0., 0., 1., 0.])
    Ralven_m = np.asarray([0., 0., 0., -1., 0., 0., 1., 0.])
    Rslow_p = 1./(2.*np.sqrt(5.)) * np.asarray([4., 2., 4., 0., 0., -2., 0., 3.])
    Rslow_m = 1./(2.*np.sqrt(5.)) * np.asarray([4., -2., -4., 0., 0., -2., 0., 3.])
    Rentrop = 1./2. * np.asarray([2., 2., 0., 0., 0., 0., 0., 1.])

    R_mode = {
        "Alfven": Ralven_p,
        "Fast":   Rfast_p,
        "Slow":   Rslow_p,
        "Entropy": Rentrop,
    }

    R = R_mode[mode]

    sin_a = 2/3
    cos_a = np.sqrt(1 - sin_a**2)

    sin_b = 2/np.sqrt(5)
    cos_b = np.sqrt(1 - sin_b**2)

    e1 = np.array([cos_a*cos_b, cos_a*sin_b, sin_a])
    e2 = np.array([-sin_b, cos_b, 0.0])
    e3 = np.array([-sin_a*cos_b, -sin_a*sin_b, cos_a])

    T = np.vstack([e1, e2, e3]).T

    def rotate(wave_frame):
        return T @ wave_frame

    v_bg = rotate(background[1:4])
    b_bg = rotate(background[4:7])

    Rv = rotate(np.array(R[1:4]))
    Rb = rotate(np.array(R[4:7]))

    e_bg = background[7]/(sim.gamma - 1) + 0.5*(background[0] * np.dot(v_bg, v_bg) + np.dot(b_bg, b_bg))

    def phase(x,y,z):
        return np.cos(2.*np.pi*x1(x,y,z))

    def x1(x,y,z):
        return x*cos_a*cos_b + y*cos_a*sin_b + z*sin_a

    def density(x,y,z):
        return background[0] + eps * R[0] * phase(x,y,z)

    def rhovx(x,y,z):
        return background[0]*v_bg[0] + eps * Rv[0] * phase(x,y,z)

    def rhovy(x,y,z):
        return background[0]*v_bg[1] + eps * Rv[1] * phase(x,y,z) 

    def rhovz(x,y,z):
        return background[0]*v_bg[2] + eps * Rv[2] * phase(x,y,z) 

    def vx(x,y,z):
        return rhovx(x,y,z)/density(x,y,z)

    def vy(x,y,z):
        return rhovy(x,y,z)/density(x,y,z)

    def vz(x,y,z):
        return rhovz(x,y,z)/density(x,y,z)

    def bx(x,y,z):
        return b_bg[0] + eps * Rb[0] * phase(x,y,z) 

    def by(x,y,z):
        return b_bg[1] + eps * Rb[1] * phase(x,y,z) 

    def bz(x,y,z):
        return b_bg[2] + eps * Rb[2] * phase(x,y,z) 

    def E(x,y,z):
        return e_bg + eps * R[7] * phase(x,y,z) 

    def p(x,y,z):
        return (E(x,y,z) - 0.5*(density(x,y,z) * (vx(x,y,z)*vx(x,y,z) + vy(x,y,z)*vy(x,y,z) + vz(x,y,z)*vz(x,y,z)) + (bx(x,y,z)*bx(x,y,z) + by(x,y,z)*by(x,y,z) + bz(x,y,z)*bz(x,y,z)))) * (sim.gamma - 1)

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    return sim


# using by error is arbitrary now, add the error for everyone now
def compute_error(run, final_time, Nx, Dx, ghosts=0):
    coords = np.arange(Nx + 2 * ghosts) * Dx + 0.5 * Dx
    from pyphare.pharesee.hierarchy.hierarchy_utils import single_patch_for_LO
    computed_by = single_patch_for_LO(run.GetB(final_time, all_primal=False).By).levels()[0].patches[0].patch_datas["By"].dataset[:]

    expected_by = single_patch_for_LO(run.GetB(0., all_primal=False).By).levels()[0].patches[0].patch_datas["By"].dataset[:]

    # expected_by = 1e-6 * np.cos(2 * np.pi * (coords - final_time))
    return np.sum(np.abs(computed_by - expected_by)) / len(computed_by)


def run_convergence(reconstruction, limiter):
    N_base = 16
    dx_values, errors, N_values = [], [], []

    while N_base <= 128:
        Nx, Ny, Nz = 2*N_base, N_base, N_base
        Dx, Dy, Dz = 3.0/Nx, 1.5/Ny, 1.5/Nz

        ph.global_vars.sim = None
        Simulator(config(N_base, reconstruction, limiter)).run().reset()

        run = Run(diag_dir)
        error = compute_error(run, final_time, Nx, Dx)

        dx_values.append(Dx)
        N_values.append(N_base)
        errors.append(error)

        N_base *= 2

    dx_values = np.array(dx_values)
    slope, intercept = np.polyfit(np.log(dx_values), np.log(errors), 1)
    expected = expected_orders[reconstruction]

    fitted_line = np.exp(intercept) * dx_values**slope
    plt.figure(figsize=(10, 6))
    plt.loglog(dx_values, errors, "o-", label=f"Data (Slope: {slope:.2f})")
    plt.loglog(dx_values, fitted_line, "--", label="Fitted Line")
    plt.xlabel("Δx", fontsize=16)
    plt.ylabel("Error (L1 Norm)", fontsize=16)
    plt.title(f"{mode} - {reconstruction}", fontsize=20)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend(fontsize=20)
    plt.savefig(f"{diag_dir}/convergence_{reconstruction}.png", dpi=200)
    plt.close()

    relative_error = abs(slope - expected) / abs(expected)
    assert relative_error < tolerance, f"{reconstruction}: got {slope}, expected {expected}"


def main():
    for reconstruction, limiter in limiters.items():
        run_convergence(reconstruction, limiter)


if __name__ == "__main__":
    main()
