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

time_step = 5e-4
final_time = 1.0
timestamps = [0.0, final_time]
diag_dir = "phare_outputs/convergence"

# Expected orders for different reconstructions
expected_orders = {
    "Constant": 1.0,
    "Linear": 2.0,
    "WENO3": 3.0,
    "WENOZ": 5.0,
    "MP5": 5.0,
}

# Limiter per reconstruction (limiters are only valid with Linear).
limiters = {
    "Constant": "None",
    "Linear": "VanLeer",
    "WENO3": "None",
    "WENOZ": "None",
    "MP5": "None",
}

# SSPRK4_5 so the temporal error never dominates the spatial convergence.
mhd_timestepper = "SSPRK4_5"
ghosts = 2

tolerance = 0.15


def config(nx, dx, reconstruction, limiter):
    sim = ph.Simulation(
        smallest_patch_size=15,
        # largest_patch_size=25,
        time_step=time_step,
        final_time=final_time,
        cells=(nx,),
        dl=(dx,),        refinement="tagging",
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

    def density(x):
        return 1.0

    def vx(x):
        return 0.0

    def vy(x):
        return -1e-6 * np.cos(2 * np.pi * x)

    def vz(x):
        return 0.0

    def bx(x):
        return 1.0

    def by(x):
        return 1e-6 * np.cos(2 * np.pi * x)

    def bz(x):
        return 0.0

    def p(x):
        return 0.1

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    return sim


def compute_error(run, final_time, Nx, Dx, ghosts=0):
    coords = np.arange(Nx + 2 * ghosts) * Dx + 0.5 * Dx
    from pyphare.pharesee.hierarchy.hierarchy_utils import single_patch_for_LO
    computed_by = single_patch_for_LO(run.GetB(final_time, all_primal=False).By).levels()[0].patches[0].patch_datas["By"].dataset[:]

    expected_by = single_patch_for_LO(run.GetB(0., all_primal=False).By).levels()[0].patches[0].patch_datas["By"].dataset[:]

    # expected_by = 1e-6 * np.cos(2 * np.pi * (coords - final_time))
    return np.sum(np.abs(computed_by - expected_by)) / len(computed_by)


def run_convergence(reconstruction, limiter):
    Nx0 = 50
    Dx0 = 1.0 / Nx0
    Nx, Dx = Nx0, Dx0

    dx_values, errors = [], []

    while Dx > Dx0 / 32.0 and Nx < 1600:
        ph.global_vars.sim = None
        sim = config(Nx, Dx, reconstruction, limiter)
        Simulator(sim).run().reset()
        if sim.dry_run:
            return
        run = Run(diag_dir)
        error = compute_error(run, final_time, Nx, Dx, ghosts)
        dx_values.append(Dx)
        errors.append(error)
        Dx /= 2.0
        Nx *= 2

    dx_values = np.array(dx_values)
    log_dx = np.log(dx_values)
    log_errors = np.log(errors)
    slope, intercept = np.polyfit(log_dx, log_errors, 1)
    expected = expected_orders[reconstruction]

    fitted_line = np.exp(intercept) * dx_values**slope
    plt.figure(figsize=(10, 6))
    plt.loglog(dx_values, errors, "o-", label=f"Data (Slope: {slope:.2f})")
    plt.loglog(dx_values, fitted_line, "--", label="Fitted Line")
    plt.xlabel("Δx", fontsize=16)
    plt.ylabel("Error (L1 Norm)", fontsize=16)
    plt.title(f"Convergence Plot - {reconstruction}", fontsize=20)
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

