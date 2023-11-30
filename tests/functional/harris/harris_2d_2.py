#!/usr/bin/env python3

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI

import numpy as np

# import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use("Agg")

from pyphare.cpp import cpp_lib
cpp = cpp_lib()
# startMPI()

diag_outputs = "phare_outputs/test/harris/2d"
# from datetime import datetime


def config():
    L = 0.5
    sim = ph.Simulation(
        time_step=0.005,
        final_time=30.0,
        # boundary_types="periodic",
        cells=(800, 400),
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=3,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="10",
        hyper_resistivity=0.002,
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": ".", "mode": "overwrite"},
        },
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
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly

        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        return dBy1 + dBy2

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
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

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return np.sqrt(T(x, y))

    def vthy(x, y):
        return np.sqrt(T(x, y))

    def vthz(x, y):
        return np.sqrt(T(x, y))

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
        "nbr_part_per_cell": 100,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
        },  # , "init":{"seed":cpp.mpi_rank()+int(sys.argv[1])}}
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    dt = 10 * sim.time_step
    nt = sim.final_time / dt + 1
    timestamps = dt * np.arange(nt)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

    return sim


def get_time(path, time, datahier=None):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    datahier = hierarchy_from(h5_filename=path + "/EM_E.h5", time=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path + "/EM_B.h5", time=time, hier=datahier)
    return datahier


def post_advance(new_time):
    if cpp.mpi_rank() == 0:
        print(f"running tests at time {new_time}")
        from tests.simulator.test_advance import AdvanceTestBase

        test = AdvanceTestBase()
        test.base_test_overlaped_fields_are_equal(
            get_time(diag_outputs, new_time), new_time
        )
        print(f"tests passed")


def main():
    sim = config()
    s = Simulator(sim)  # , post_advance=post_advance)
    s.initialize()
    post_advance(0)
    s.run()


if __name__ == "__main__":
    main()
