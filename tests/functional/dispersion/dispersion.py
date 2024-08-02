#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator


import matplotlib as mpl
import numpy as np

mpl.use("Agg")


def fromNoise():
    # in this configuration there are no prescribed waves
    # and only eigen modes existing in the simulation noise
    # will be visible

    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=20,
        # the following time step number
        # and final time mean that the
        # smallest frequency will be 2/100
        # and the largest 2/dt  = 2e3
        time_step_nbr=100000,
        final_time=100.0,
        boundary_types="periodic",
        # smallest wavelength will be 2*0.2=0.4
        # and largest 50
        cells=500,
        dl=0.2,
        diag_options={
            "format": "phareh5",
            "options": {"dir": "dispersion", "mode": "overwrite"},
        },
    )

    def density(x):
        return 1.0

    def by(x):
        return 0.0

    def bz(x):
        return 0.0

    def bx(x):
        return 1.0

    def vx(x):
        return 0.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vthx(x):
        return 0.01

    def vthy(x):
        return 0.01

    def vthz(x):
        return 0.01

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": 1, "density": density, **vvv}
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time + sim.time_step, sim.time_step)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def prescribedModes():
    # in this configuration there user prescribed
    # wavelength, at which more energy is thus expected
    # than in modes naturally present in the simulation noise

    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=20,
        # the following time step number
        # and final time mean that the
        # smallest frequency will be 2/100
        # and the largest 2/dt  = 2e3
        time_step_nbr=100000,
        final_time=100.0,
        boundary_types="periodic",
        # smallest wavelength will be 2*0.2=0.4
        # and largest 50
        cells=500,
        dl=0.2,
        diag_options={
            "format": "phareh5",
            "options": {"dir": "dispersion", "mode": "overwrite"},
        },
    )

    def density(x):
        return 1.0

    def by(x):
        L = sim.simulation_domain()
        return 0.1 * np.cos(2 * np.pi * x / L[0])

    def bz(x):
        L = sim.simulation_domain()
        return -0.1 * np.sin(2 * np.pi * x / L[0])

    def bx(x):
        return 1.0

    def vx(x):
        return 0.0

    def vy(x):
        L = sim.simulation_domain()
        return 0.1 * np.cos(2 * np.pi * x / L[0])

    def vz(x):
        L = sim.simulation_domain()
        return 0.1 * np.sin(2 * np.pi * x / L[0])

    def vthx(x):
        return 0.01

    def vthy(x):
        return 0.01

    def vthz(x):
        return 0.01

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": 1, "density": density, **vvv}
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    timestamps = np.arange(0, sim.final_time + sim.time_step, sim.time_step)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def main():
    Simulator(fromNoise()).run()


if __name__ == "__main__":
    main()
