#!/usr/bin/env python3


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest

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
    return np.sqrt(T(x))


def vthy(x):
    return np.sqrt(T(x))


def vthz(x):
    return np.sqrt(T(x))


vvv = {
    "vbulkx": vx,
    "vbulky": vy,
    "vbulkz": vz,
    "vthx": vthx,
    "vthy": vthy,
    "vthz": vthz,
}


# the only dump time, at which all checks are done
final_time = 11.0
nppc = 100


def simulation_params(diagdir, **extra):
    params = {
        "interp_order": 1,
        "final_time": final_time,
        "time_step": 0.04,
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
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            "nbr_part_per_cell": nppc,
            **vvv,
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)

    timestamps = [final_time]

    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["charge_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for pop in sim.model.populations:
        for quantity in ["domain"]:
            ph.ParticleDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                population_name=pop,
            )

    ph.LoadBalancer(auto=True, mode="nppc", tol=0.05)

    return sim


def withTagging(diagdir):
    return config(**simulation_params(diagdir, refinement="tagging", max_nbr_levels=3))


def noRefinement(diagdir):
    return config(**simulation_params(diagdir))


def make_figure():
    from scipy.optimize import curve_fit

    from pyphare.pharesee.plotting import zoom_effect

    rwT = Run("./withTagging")
    rNoRef = Run("./noRefinement")

    params = simulation_params("")
    L = params["cells"] * params["dl"]
    plot_time = final_time
    v = vx(0)

    BH = rwT.GetB(plot_time)
    BywT = rwT.GetB(plot_time, all_primal=False).finest(plot_time, "By", "linear")
    ByNoRef = rNoRef.GetB(plot_time, all_primal=False).finest(plot_time, "By", "linear")
    JzwT = rwT.GetJ(plot_time, all_primal=False).finest(plot_time, "Jz", "linear")
    JzNoRef = rNoRef.GetJ(plot_time, all_primal=False).finest(plot_time, "Jz", "linear")

    xbywT, bywT = BywT.x, BywT[:]
    xbyNoRef, byNoRef = ByNoRef.x, ByNoRef[:]
    xjzwT, jzwT = JzwT.x, JzwT[:]
    xjzNoRef, jzNoRef = JzNoRef.x, JzNoRef[:]

    def by_expected(x):
        # initial profile carried at v without deformation, periodic domain
        x0 = np.mod(x - v * plot_time, L)
        return -1 + 2 * (S(x0, L * 0.25, 1) - S(x0, L * 0.75, 1))

    # sheets start at L/4 and 3L/4, zoom windows follow them
    sheets = [np.mod(x0 + v * plot_time, L) for x0 in (0.25 * L, 0.75 * L)]
    windows = [(xs - 22, xs + 23) for xs in sheets]

    fig, axd = plt.subplot_mosaic(
        [["full", "full"], ["byL", "byR"], ["jzL", "jzR"]],
        figsize=(11, 9),
        layout="constrained",
    )
    fig.suptitle(
        f"1D tangential discontinuities advected at V = {v:g}, t = {plot_time:g}"
    )

    ax0 = axd["full"]
    ax0.plot(xbyNoRef, byNoRef, color="k", label="no refinement")
    ax0.plot(xbywT, bywT, color="royalblue", label="tagging, 3 levels")
    ax0.plot(
        xbyNoRef,
        by_expected(xbyNoRef),
        color="darkorange",
        ls="--",
        label="initial profile moved by V t",
    )
    ax0.set_ylabel("By")
    ax0.set_xlim(0, L)

    for side, (xmin, xmax), xs in zip("LR", windows, sheets):
        axby, axjz = axd["by" + side], axd["jz" + side]
        axby.plot(xbyNoRef, byNoRef, color="k")
        axby.plot(xbywT, bywT, color="royalblue")
        axby.plot(xbyNoRef, by_expected(xbyNoRef), color="darkorange", ls="--")
        axby.set_ylim((-1.5, 2))
        axby.set_ylabel("By")
        axjz.plot(xjzNoRef, jzNoRef, color="k")
        axjz.plot(xjzwT, jzwT, color="royalblue")
        axjz.set_ylabel("Jz")
        axjz.set_xlabel("x")
        for ax in (axby, axjz):
            ax.set_xlim((xmin, xmax))
        zoom_effect(ax0, axby, xmin, xmax)
        for ax in (ax0, axby, axjz):
            ax.axvline(xs, color="r", lw=1)
    axd["jzL"].set_ylim((-0.5, 1.5))
    axd["jzR"].set_ylim((-1.5, 0.5))

    # refined patches of the tagged run, level i drawn in the band [i/4, (i+1)/4]
    for ilvl, level in BH.levels().items():
        for patch in level.patches:
            x0 = patch.origin[0]
            x1 = (patch.box.upper[0] + 1) * patch.layout.dl[0]
            for ax in axd.values():
                ax.axvspan(
                    x0,
                    x1,
                    color="royalblue",
                    ec="k",
                    alpha=0.2,
                    ymin=ilvl / 4,
                    ymax=(ilvl + 1) / 4,
                )

    handles, labels = ax0.get_legend_handles_labels()
    handles += [
        plt.Line2D([], [], color="r", lw=1),
        mpl.patches.Patch(color="royalblue", alpha=0.2, ec="k"),
    ]
    labels += ["expected sheet position", "patches of levels 0, 1, 2 (bottom to top)"]
    fig.legend(handles, labels, loc="outside lower center", ncol=5, fontsize="small")

    fig.savefig("tdtagged1d.png")

    # select data around the rightward TD
    xmin, xmax = windows[1]
    idx = np.where((xbywT > xmin) & (xbywT < xmax))
    xx = xbywT[idx]
    bby = bywT[idx]

    # fit the right sheet: expect thickness 1 at x = 3L/4 + V t
    def by_fit(x, x0, l):
        v1 = 1
        v2 = -1
        return v1 + (v2 - v1) * S(x, x0, l)

    popt, _pcov = curve_fit(by_fit, xx, bby, p0=(xmin, 1))
    x0, l = popt

    if np.abs(l - 1) > 0.5:
        raise RuntimeError(f"L (={l}) too far from 1.0")
    if np.abs(x0 - sheets[1]) > 0.5:
        raise RuntimeError(f"x0 (={x0}) too far from {sheets[1]}")


def check_refined_particles(diagdir):
    from pyphare.pharesee.hierarchy import hierarchy_from

    datahier = hierarchy_from(
        h5_filename=diagdir + "/ions_pop_protons_domain.h5", times=f"{final_time:.10f}"
    )
    HybridAdvanceTest().base_test_domain_particles_on_refined_level(
        datahier, final_time
    )


def check_invariants(diagdir):
    """exact on the periodic L0: particle count, mean By (= 0) and mean Bz (= 0.5)"""
    run = Run(diagdir)

    cells = simulation_params(diagdir)["cells"]
    n_L0 = sum(
        patch.patch_datas["protons_domain"].size()
        for patch in run.GetParticles(final_time, "protons").level(0).patches
    )
    if n_L0 != cells * nppc:
        raise RuntimeError(f"L0 particle count {n_L0} != {cells * nppc}")

    B = run.GetB(final_time, all_primal=False)  # native grid, no duplicated nodes
    for qty, expected in (("By", 0.0), ("Bz", 0.5)):
        values = np.concatenate(
            [
                patch.patch_datas[qty][patch.patch_datas[qty].box]
                for patch in B.level(0).patches
            ]
        )
        # conserved to round-off, but dumps may be single precision
        tol = max(np.finfo(values.dtype).eps, 1e-12)
        mean = np.mean(values, dtype=np.float64)
        if np.abs(mean - expected) > tol:
            raise RuntimeError(f"L0 mean {qty} = {mean} != {expected}")


def main():
    sim = noRefinement(diagdir="noRefinement")
    Simulator(sim).run().reset()
    ph.global_vars.sim = None

    Simulator(withTagging(diagdir="withTagging")).run().reset()
    ph.global_vars.sim = None

    if not sim.dry_run and cpp.mpi_rank() == 0:
        make_figure()
        check_refined_particles("withTagging")
        for diagdir in ("noRefinement", "withTagging"):
            check_invariants(diagdir)


if __name__ == "__main__":
    main()
