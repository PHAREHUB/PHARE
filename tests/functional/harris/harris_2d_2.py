#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI

import numpy as np
import matplotlib as mpl

mpl.use("Agg")


def config():
    L = 0.5
    sim = ph.Simulation(
        time_step=0.005,
        time_step_nbr=5,
        cells=(800, 400),
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=1,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="10",
        hyper_resistivity=0.002,
        resistivity=0.001,
        strict=True,
    )

    def density(x, y):
        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        Lx, Ly = sim.simulation_domain()
        sigma = 1.0
        dB = 0.1
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)
        return dBy1 + dBy2

    def bx(x, y):
        Lx, Ly = sim.simulation_domain()
        sigma = 1.0
        dB = 0.1
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)
        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + dBx1 + dBx2

    def bz(x, y):
        return 0.0

    def b2(x, y):
        return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

    def T(x, y):
        K = 0.7
        temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
        assert np.all(temp > 0)
        return temp

    def vxyz(x, y):
        return 0.0

    def vthxyz(x, y):
        return np.sqrt(T(x, y))

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": vxyz,
            "vbulky": vxyz,
            "vbulkz": vxyz,
            "vthx": vthxyz,
            "vthy": vthxyz,
            "vthz": vthxyz,
            "nbr_part_per_cell": 100,
            "init": {"seed": 1333337},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    return sim


if __name__ == "__main__":
    Simulator(config()).run()
