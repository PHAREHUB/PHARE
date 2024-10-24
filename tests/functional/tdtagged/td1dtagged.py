#!/usr/bin/env python3


import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.run import Run


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from tests.simulator.test_advance import AdvanceTestBase
from tests.diagnostic import all_timestamps
from pyphare.cpp import cpp_lib

cpp = cpp_lib()

mpl.use("Agg")


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
    return 2.0


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


# used to only test on the early particle diagnostic files
particle_diagnostics = {"count": 10, "idx": 0}


def simulation_params(diagdir, **extra):
    params = {
        "interp_order": 1,
        "time_step_nbr": 500,
        "time_step": 0.04,
        "boundary_types": "periodic",
        "cells": 200,
        "hyper_resistivity": 0.01,
        "dl": 1.0,
        "diag_options": {
            "format": "phareh5",
            "options": {"dir": diagdir, "mode": "overwrite"},
        },
    }
    params.update(**extra)
    return params


def config(**options):
    sim = ph.Simulation(**options)
    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz, protons={"charge": 1, "density": density, **vvv}
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)

    timestamps = all_timestamps(sim)

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for pop in sim.model.populations:
        for quantity in ["domain"]:
            ph.ParticleDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps[: particle_diagnostics["count"] + 1],
                population_name=pop,
            )

    return sim


def withTagging(diagdir):
    return config(**simulation_params(diagdir, refinement="tagging", max_nbr_levels=3))


def noRefinement(diagdir):
    return config(**simulation_params(diagdir))


def make_figure():
    from scipy.optimize import curve_fit

    rwT = Run("./withTagging")
    rNoRef = Run("./noRefinement")

    plot_time = 11
    v = 2

    BH = rwT.GetB(plot_time)
    BwT = rwT.GetB(plot_time, merged=True, interp="linear")
    BNoRef = rNoRef.GetB(plot_time, merged=True, interp="linear")
    JwT = rwT.GetJ(plot_time, merged=True, interp="linear")
    JNoRef = rNoRef.GetJ(plot_time, merged=True, interp="linear")

    xbywT = BwT["By"][1][0]
    bywT = BwT["By"][0](xbywT)
    xbyNoRef = BNoRef["By"][1][0]
    byNoRef = BNoRef["By"][0](xbyNoRef)
    xjzwT = JwT["Jz"][1][0]
    jzwT = JwT["Jz"][0](xjzwT)
    xjzNoRef = JNoRef["Jz"][1][0]
    jzNoRef = JNoRef["Jz"][0](xjzNoRef)

    fig, axarr = plt.subplots(nrows=3, figsize=(8, 8))

    def S(x, x0, l):
        return 0.5 * (1 + np.tanh((x - x0) / l))

    def by(x):
        L = 200
        v1 = -1
        v2 = 1
        return v1 + (v2 - v1) * (S(x, L * 0.25, 1) - S(x, L * 0.75, 1))

    wT0 = 150.0

    ax0, ax1, ax2 = axarr

    ax0.plot(xbyNoRef, byNoRef, color="k", ls="-")
    ax0.plot(xbywT, bywT, color="royalblue", ls="-")
    ax0.plot(xbyNoRef, by(xbyNoRef), color="darkorange", ls="--")

    ax1.set_xlim((wT0, 195))
    ax1.set_ylim((-1.5, 2))
    ax1.plot(xbyNoRef, byNoRef, color="k", ls="-")
    ax1.plot(xbywT, bywT, color="royalblue", ls="-")

    ax2.plot(xjzwT, jzwT)
    ax2.plot(xjzNoRef, jzNoRef, color="k")
    ax2.set_xlim((wT0, 195))
    ax2.set_ylim((-1.5, 0.5))

    # draw level patches
    for ilvl, level in BH.levels().items():
        for patch in level.patches:
            dx = patch.layout.dl[0]
            x0 = patch.origin[0]
            x1 = (patch.box.upper[0] + 1) * patch.layout.dl[0]
            for ax in (ax1, ax2, ax0):
                ax.axvspan(
                    x0,
                    x1,
                    color="royalblue",
                    ec="k",
                    alpha=0.2,
                    ymin=ilvl / 4,
                    ymax=(ilvl + 1) / 4,
                )

    from pyphare.pharesee.plotting import zoom_effect

    zoom_effect(ax0, ax1, wT0, 195)

    for ax in (ax0, ax1, ax2):
        ax.axvline(wT0 + plot_time * v, color="r")

    fig.savefig("tdtagged1d.png")

    # select data around the rightward TD
    idx = np.where((xbywT > 150) & (xbywT < 190))
    xx = xbywT[idx]
    bby = bywT[idx]

    # now we will fit by_fit to the data
    # and we expect to find x0=172 and L=1
    # or close enough
    def by_fit(x, x0, L):
        v1 = 1
        v2 = -1
        return v1 + (v2 - v1) * S(x, x0, L)

    popt, pcov = curve_fit(by_fit, xx, bby, p0=(150, 1))
    x0, L = popt

    if np.abs(L - 1) > 0.5:
        raise RuntimeError(f"L (={L}) too far from 1.O")
    if np.abs(x0 - (150 + plot_time * v)) > 0.5:
        raise RuntimeError(f"x0 (={x0}) too far from 172")


test = AdvanceTestBase()


def get_time(path, time):
    if time is not None:
        time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from

    return hierarchy_from(h5_filename=path + "/ions_pop_protons_domain.h5", times=time)


def post_advance(new_time):
    if (
        particle_diagnostics["idx"] < particle_diagnostics["count"]
        and cpp.mpi_rank() == 0
    ):
        particle_diagnostics["idx"] += 1
        datahier = get_time(ph.global_vars.sim.diag_options["options"]["dir"], new_time)
        test.base_test_domain_particles_on_refined_level(datahier, new_time)


def main():
    Simulator(noRefinement(diagdir="noRefinement")).run().reset()
    ph.global_vars.sim = None

    Simulator(
        withTagging(diagdir="withTagging"), post_advance=post_advance
    ).run().reset()
    ph.global_vars.sim = None

    if cpp.mpi_rank() == 0:
        make_figure()


if __name__ == "__main__":
    main()
