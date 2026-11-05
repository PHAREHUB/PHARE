#!/usr/bin/env python3
# Copyright 2026 Finnish Meteorological Institute

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator


import matplotlib as mpl
import numpy as np

mpl.use("Agg")


def shock_tube():
    species_mass = 1.6726217770000001e-27
    particle_temp_nrj_ratio = 1.38064852e-23
    temperature = 7.243e4
    elementary_charge = 1.602176634e-19

    sim = ph.Simulation(
        time_step_nbr=1000,
        final_time=0.1,
        boundary_types="periodic",
        cells=5000,
        dl=20,
        diag_options={
            "format": "phareh5",
            "options": {"dir": "shock_tube", "mode": "overwrite"},
        },
    )

    def density(x):
        return np.where(np.abs(x - 1e5/2) < 1e5/4, 3e6, 1e6)

    def bx(x):
        return 1.5e-9

    def by(x):
        return np.where(np.abs(x - 1e5/2) < 1e5/4, -1e-9, +1e-9)

    def bz(x):
        return 0.0

    def vx(x):
        return 0.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vthx(x):
        return np.sqrt(particle_temp_nrj_ratio * temperature / species_mass)

    def vthy(x):
        return np.sqrt(particle_temp_nrj_ratio * temperature / species_mass)

    def vthz(x):
        return np.sqrt(particle_temp_nrj_ratio * temperature / species_mass)

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": elementary_charge, "density": density, **vvv}
    )

    ph.ElectronModel(closure = "isothermal", Te = temperature)

    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time + sim.time_step, sim.time_step)[::100]

    ph.ElectromagDiagnostics(quantity = 'B', write_timestamps = timestamps)
    ph.ElectromagDiagnostics(quantity = 'E', write_timestamps = timestamps)
    ph.FluidDiagnostics(quantity = 'mass_density', write_timestamps = timestamps)
    ph.FluidDiagnostics(quantity = 'bulkVelocity', write_timestamps = timestamps)

    return sim


def main():
    Simulator(shock_tube()).run()


if __name__ == "__main__":
    main()
