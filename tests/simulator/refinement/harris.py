#!/usr/bin/env python3

import pyphare.pharein as ph  # lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import (
    ElectromagDiagnostics,
    FluidDiagnostics,
    ParticleDiagnostics,
    InfoDiagnostics,
)
from pyphare.pharein import MetaDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharein import global_vars as gv
from pyphare.pharein import LoadBalancer
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tests.simulator.test_advance import AdvanceTestBase

mpl.use("Agg")
from pyphare.cpp import cpp_lib

cpp = cpp_lib()
startMPI()

test = AdvanceTestBase()


def config():

    start_time = 0.0
    L = 0.5
    Simulation(
        time_step=0.010,
        final_time=0.15,
        # boundary_types="periodic",
        cells=(80, 40),
        dl=(0.40, 0.40),
        refinement="tagging",
        max_nbr_levels=1,
        nesting_buffer=1,
        clustering="tile",
        tag_buffer="10",
        tagging_threshold=0.4,
        hyper_resistivity=0.000,
        # hyper_mode="constant",
        resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": "diags_master", "mode": "overwrite"},
        },
        restart_options={
            "dir": "checkpoints",
            "mode": "overwrite",
            "elapsed_timestamps": [36000, 79000],
        },
    )  # ,"restart_time":start_time }

    def density(x, y):
        from pyphare.pharein.global_vars import sim

        Ly = sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        from pyphare.pharein.global_vars import sim

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
        from pyphare.pharein.global_vars import sim

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

    MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            **vvv,
            "init": {"seed": cpp.mpi_rank() + 11},
        },
    )

    ElectronModel(closure="isothermal", Te=0.0)

    LoadBalancer(active=True, mode="nppc", tol=0.05, every=1000)

    sim = ph.global_vars.sim
    dt = 1.0 * sim.time_step
    nt = (sim.final_time - start_time) / dt
    timestamps = start_time + dt * np.arange(nt)
    print(timestamps)

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
        )

    for quantity in ["density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
        )

    InfoDiagnostics(quantity="particle_count", write_timestamps=timestamps)


def check_hier(hier):
    for ilvl, lvl in hier.levels().items():
        for patch in lvl.patches:
            assert all(patch.box.shape > 5)
    return hier


def get_time(path, time=None, datahier=None):
    if time is not None:
        time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    datahier = hierarchy_from(h5_filename=path + "/EM_E.h5", times=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path + "/EM_B.h5", times=time, hier=datahier)
    return datahier


def get_hier(path):
    return get_time(path)


def post_advance(new_time):
    datahier = check_hier(
        get_hier("/home/aunai/Documents/code/phare/phare_jobs/run104/diags_master")
    )
    errors = test.base_test_overlaped_fields_are_equal(datahier, new_time)


def main():

    config()
    simulator = Simulator(gv.sim, print_one_line=False, post_advance=post_advance)
    simulator.initialize()
    simulator.run()


if __name__ == "__main__":
    main()
