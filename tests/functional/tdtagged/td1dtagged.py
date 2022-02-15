#!/usr/bin/env python3


import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
from pyphare.pharesee.run import Run


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')


from tests.diagnostic import all_timestamps


def density(x):
    return 1.


def S(x,x0,l):
    return 0.5*(1+np.tanh((x-x0)/l))


def bx(x):
    return 0.


def by(x):
    from pyphare.pharein.global_vars import sim
    L = sim.simulation_domain()[0]
    v1=-1
    v2=1.
    return v1 + (v2-v1)*(S(x,L*0.25,1) -S(x, L*0.75, 1))


def bz(x):
    return 0.5


def b2(x):
    return bx(x)**2 + by(x)**2 + bz(x)**2


def T(x):
    K = 1
    return 1/density(x)*(K - b2(x)*0.5)


def vx(x):
    return 2.


def vy(x):
    return 0.


def vz(x):
    return 0.


def vthx(x):
    return T(x)


def vthy(x):
    return T(x)


def vthz(x):
    return T(x)

vvv = {"vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz }


# used to only test on the early particle diagnostic files
particle_diagnostics = {"count":10, "idx":0}

def simulation_params(diagdir, **extra):
    params = {
        "interp_order": 1,
        "time_step_nbr":2000,
        "time_step":.01,
        "boundary_types":"periodic",
        "cells":200,
        "dl":1.0,
        "diag_options":{"format": "phareh5",
                        "options": {
                              "dir": diagdir,
                              "mode":"overwrite"}}
    }
    params.update(**extra)
    return params


def config(**options):

    ph.Simulation(**options)
    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, **vvv}
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)

    sim = ph.global_vars.sim

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

    for pop in sim.model.populations:
        for quantity in ['domain']:
            ph.ParticleDiagnostics(quantity=quantity,
                                compute_timestamps=timestamps[:particle_diagnostics["count"]+1],
                                write_timestamps=timestamps[:particle_diagnostics["count"]+1],
                                population_name=pop)



def withTagging(diagdir):
    config(**simulation_params(diagdir,
                               refinement="tagging",
                               max_nbr_levels=3))


def staticRefinement(diagdir):
    config(**simulation_params(diagdir,
                               refinement_boxes={"L0": [Box1D(46, 52), Box1D(146, 152)]}))


def noRefinement(diagdir):
    config(**simulation_params(diagdir))




def make_figure():
    from scipy.optimize import curve_fit
    from pyphare.pharesee.hierarchy import flat_finest_field

    rwT    = Run("./withTagging")
    rNoRef = Run("./noRefinement")

    plot_time = 11
    v = 2

    BwT = rwT.GetB(plot_time)
    BNoRef = rNoRef.GetB(plot_time)
    JwT = rwT.GetJ(plot_time)
    JNoRef = rNoRef.GetJ(plot_time)

    bywT, xbywT  = flat_finest_field(BwT, "By")
    byNoRef, xbyNoRef  = flat_finest_field(BNoRef, "By")
    jzwT, xjzwT  = flat_finest_field(JwT, "Jz")
    jzNoRef, xjzNoRef  = flat_finest_field(JNoRef, "Jz")

    fig, axarr = plt.subplots(nrows=3, figsize=(8,8))

    def S(x, x0, l):
        return 0.5*(1+np.tanh((x-x0)/l))

    def by(x):
        L=200
        v1=-1
        v2=1
        return v1 + (v2-v1)*(S(x, L*0.25, 1)-S(x, L*0.75, 1))

    wT0 = 150.
    ax0, ax1, ax2 = axarr
    ax0.plot(xbywT, bywT)
    ax0.plot(xbyNoRef, byNoRef, color="k", alpha=0.6)
    ax0.plot(xbyNoRef, by(xbyNoRef), ls='--')

    ax1.plot(xbywT, bywT)
    ax1.plot(xbyNoRef, byNoRef, color='k')
    ax1.set_xlim((wT0,195))
    ax1.set_ylim((-1.5, 2))

    ax2.plot(xjzwT, jzwT)
    ax2.plot(xjzNoRef, jzNoRef, color='k')
    ax2.set_xlim((wT0,195))
    ax2.set_ylim((-1.5, 0.5))


    # draw level patches
    for ilvl,level in BwT.levels().items():
        for patch in level.patches:
            dx = patch.layout.dl[0]
            x0 = patch.origin[0]
            x1 = (patch.box.upper[0]+1)*patch.layout.dl[0]
            for ax in (ax1, ax2, ax0):
                ax.axvspan(x0, x1, color='b',ec='k', alpha=0.2,
                            ymin=ilvl/4, ymax=(ilvl+1)/4)


    from pyphare.pharesee.plotting import zoom_effect
    zoom_effect(ax0, ax1, wT0, 195)

    for ax in (ax0, ax1, ax2):
        ax.axvline(wT0+plot_time*v, color="r")

    fig.savefig("td1dtagged.png")

    # select data around the rightward TD
    idx = np.where((xbywT>150) & (xbywT<190))
    xx = xbywT[idx]
    bby = bywT[idx]


    # now we will fit by_fit to the data
    # and we expect to find x0=172 and L=1
    # or close enough
    def by_fit(x, x0,L):
        v1=1
        v2=-1
        return v1 + (v2-v1)*S(x, x0, L)

    popt,pcov = curve_fit(by_fit, xx, bby, p0=(150,1))
    x0, L = popt

    if np.abs(L-1)>0.5:
        raise RuntimeError(f"L (={L}) too far from 1.O")
    if np.abs(x0-(150+plot_time*v))>0.5:
        raise RuntimeError(f"x0 (={x0}) too far from 172")


from tests.simulator.test_advance import AdvanceTestBase
from pyphare.cpp import cpp_lib
cpp = cpp_lib()


test = AdvanceTestBase()

def get_time(path, time):
    if time is not None:
        time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from
    return hierarchy_from(h5_filename=path+"/ions_pop_protons_domain.h5", time=time)



def post_advance(new_time):
    if particle_diagnostics["idx"] < particle_diagnostics["count"] and cpp.mpi_rank() == 0:
        particle_diagnostics["idx"] += 1
        datahier = get_time(gv.sim.diag_options["options"]["dir"], new_time)
        test.base_test_domain_particles_on_refined_level(datahier, new_time)



def main():

    noRefinement(diagdir="noRefinement")
    Simulator(gv.sim).run()
    gv.sim = None

    withTagging(diagdir="withTagging")
    Simulator(gv.sim, post_advance=post_advance).run()
    gv.sim = None

    if cpp.mpi_rank() == 0:
        make_figure()




if __name__=="__main__":
    main()
