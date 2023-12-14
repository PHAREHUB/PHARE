#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("Agg")


from tests.diagnostic import all_timestamps


def config_uni(**kwargs):
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
            "L0": {"B0": [(100,), (200,)]},
            "L1": {"B0": [(300,), (350,)]},
        },
        diag_options={
            "format": "phareh5",
            "options": {"dir": kwargs["diagdir"], "mode": "overwrite"},
        },
    )

    def density(x):
        return 1.0

    def bx(x):
        return 0.0

    def by(x):
        return 1.0

    def bz(x):
        return 0.5

    def vx(x):
        return kwargs["vx"]

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def vthx(x):
        return 0.1

    def vthy(x):
        return 0.1

    def vthz(x):
        return 0.1

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

    timestamps = all_timestamps(sim)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    return sim


def config_td(**kwargs):
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
        cells=200,  # integer or tuple length == dimension
        dl=1.0,  # mesh size of the root level, float or tuple
        refinement_boxes={
            "L0": {"B0": [(50,), (150,)]},
            "L1": {"B0": [(125,), (175,)]},
        },
        diag_options={
            "format": "phareh5",
            "options": {"dir": kwargs["diagdir"], "mode": "overwrite"},
        },
    )

    def density(x):
        return 1.0

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def bx(x):
        return 0.0

    def by(x):
        from pyphare.pharein.global_vars import sim

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
        return kwargs["vx"]

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

    timestamps = all_timestamps(sim)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    return sim


def main():
    for name, config in zip(("uni", "td"), (config_uni, config_td)):
        params = [
            {"vx": -1, "diagdir": name + "_vxm2"},
            {"vx": 2, "diagdir": name + "_vx2"},
        ]
        for param in params:
            if param["vx"] > -1:
                continue
            print("-----------------------------------")
            print(param)
            print("-----------------------------------")

            Simulator(config(**param)).run()
            ph.global_vars.sim = None


if __name__ == "__main__":
    main()
