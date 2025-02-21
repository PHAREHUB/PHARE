#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.use("Agg")


def config():
    """Configure the simulation

    This function defines the Simulation object,
    user initialization model and diagnostics.
    """
    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=20,
        time_step_nbr=2000,  # number of time steps (not specified if time_step and final_time provided)
        final_time=20.0,  # simulation final time (not specified if time_step and time_step_nbr provided)
        boundary_types="periodic",  # boundary condition, string or tuple, length == len(cell) == len(dl)
        cells=500,  # integer or tuple length == dimension
        dl=1.0,  # mesh size of the root level, float or tuple
        refinement_boxes={
            "L0": {"B0": [(80,), (180,)]},
            "L1": {"B0": [(200,), (300,)]},
        },
        diag_options={
            "format": "phareh5",
            "options": {"dir": "td_noflow", "mode": "overwrite"},
        },
    )

    def density(x):
        return 1.0

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def bx(x):
        return 0.0

    def by(x):
        L = sim.simulation_domain()[0]
        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))

    def bz(x):
        return 0.5

    def b2(x):
        return bx(x) ** 2 + by(x) ** 2 + bz(x) ** 2

    def T(x):
        K = 1
        return 1 / density(x) * (K - b2(x) * 0.5)

    def vx(x):
        return 0.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vthx(x):
        return T(x)

    def vthy(x):
        return T(x)

    def vthz(x):
        return T(x)

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

    ph.ElectronModel(closure="isothermal", Te=0.12)

    dt_dump = 0.1
    n_dump = int(sim.final_time / dt_dump) + 1
    timestamps = np.linspace(0, sim.final_time, n_dump)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for quantity in ["charge_density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
        )
    return sim


def main():
    Simulator(config()).run()


if __name__ == "__main__":
    main()
