#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics,FluidDiagnostics, ParticleDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

from pyphare.cpp import cpp_lib
cpp = cpp_lib()
diag_outputs="phare_outputs/test/harris/2d"
from datetime import datetime

def config():

    Simulation(
        smallest_patch_size=15,
        largest_patch_size=25,
        time_step_nbr=1000,
        time_step=0.001,
        boundary_types="periodic",
        cells=(100,100),
        dl=(0.2, 0.2),
        refinement_boxes={},
        hyper_resistivity=0.001,
        resistivity=0.001,
        diag_options={"format": "phareh5",
                      "options": {"dir": diag_outputs,
                                  "mode":"overwrite"}},
        strict=True,
    )

    def density(x, y):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()[1]
        return 0.2 + 1./np.cosh((y-L*0.3)/0.5)**2 + 1./np.cosh((y-L*0.7)/0.5)**2


    def S(y, y0, l):
        return 0.5*(1. + np.tanh((y-y0)/l))


    def by(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = (x - 0.5 * Lx)
        y1 = (y - 0.3 * Ly)
        y2 = (y - 0.7 * Ly)
        w3 = np.exp(-(x0*x0 + y1*y1) / (w2*w2))
        w4 = np.exp(-(x0*x0 + y2*y2) / (w2*w2))
        w5 = 2.0*w1/w2
        return (w5 * x0 * w3) + ( -w5 * x0 * w4)


    def bx(x, y):
        from pyphare.pharein.global_vars import sim
        Lx = sim.simulation_domain()[0]
        Ly = sim.simulation_domain()[1]
        w1 = 0.2
        w2 = 1.0
        x0 = (x - 0.5 * Lx)
        y1 = (y - 0.3 * Ly)
        y2 = (y - 0.7 * Ly)
        w3 = np.exp(-(x0*x0 + y1*y1) / (w2*w2))
        w4 = np.exp(-(x0*x0 + y2*y2) / (w2*w2))
        w5 = 2.0*w1/w2
        v1=-1
        v2=1.
        return v1 + (v2-v1)*(S(y,Ly*0.3,0.5) -S(y, Ly*0.7, 0.5)) + (-w5*y1*w3) + (+w5*y2*w4)


    def bz(x, y):
        return 0.


    def b2(x, y):
        return bx(x,y)**2 + by(x, y)**2 + bz(x, y)**2


    def T(x, y):
        K = 1
        temp = 1./density(x, y)*(K - b2(x, y)*0.5)
        assert np.all(temp >0)
        return temp

    def vx(x, y):
        return 0.


    def vy(x, y):
        return 0.


    def vz(x, y):
        return 0.


    def vthx(x, y):
        return np.sqrt(T(x, y))


    def vthy(x, y):
        return np.sqrt(T(x, y))


    def vthz(x, y):
        return np.sqrt(T(x, y))


    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz,
        "nbr_part_per_cell":100
    }

    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density,  **vvv, "init":{"seed": 12334}},

    )

    ElectronModel(closure="isothermal", Te=0.0)


    sim = ph.global_vars.sim
    dt = 100 * sim.time_step
    nt = sim.final_time/dt+1
    timestamps = (dt * np.arange(nt))
    print(timestamps)


    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )


    for quantity in ["density", "bulkVelocity"]:
        FluidDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
            )

   #for popname in ("protons",):
   #    for name in ["domain", "levelGhost", "patchGhost"]:
   #        ParticleDiagnostics(quantity=name,
   #                            compute_timestamps=timestamps,
   #                            write_timestamps=timestamps,
   #                            population_name=popname)




def get_time(path, time):
    time = "{:.10f}".format(time)
    from pyphare.pharesee.hierarchy import hierarchy_from
    datahier = None
    datahier = hierarchy_from(h5_filename=path+"/EM_E.h5", time=time, hier=datahier)
    datahier = hierarchy_from(h5_filename=path+"/EM_B.h5", time=time, hier=datahier)
    # datahier = hierarchy_from(h5_filename=path+"/ions_pop_protons_domain.h5", time=time, hier=datahier)
    # datahier = hierarchy_from(h5_filename=path+"/ions_pop_protons_levelGhost.h5", time=time, hier=datahier)
    # datahier = hierarchy_from(h5_filename=path+"/ions_pop_protons_patchGhost.h5", time=time, hier=datahier)
    # from pyphare.pharesee.hierarchy import merge_particles
    # merge_particles(datahier)
    return datahier

def post_advance(new_time):
    if cpp.mpi_rank() == 0:
        from tests.simulator.test_advance import AdvanceTestBase
        test = AdvanceTestBase()
        print("opening diags")
        now = datetime.now()
        datahier = get_time(diag_outputs, new_time)
        print("diags opened in", datetime.now() - now)
        test.base_test_overlaped_fields_are_equal(datahier, new_time)
        # test.base_test_overlapped_particledatas_have_identical_particles(datahier, new_time)
        # n_particles = np.array(cells).prod() * ppc
        # n_particles_at_t = 0
        # for patch in datahier.level(0, new_time).patches:
        #     n_particles_at_t += patch.patch_datas["protons_particles"].dataset[patch.box].size()
        # test.assertEqual(n_particles, n_particles_at_t)
        # print("coarsest_time", new_time, n_particles_at_t)


import numpy as np
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pyphare
from scipy.ndimage import gaussian_filter as gf
from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.geometry import hierarchy_overlaps
from pyphare.core import box as boxm

def pdata(ip, hier, name):
    x  = hier.patch_levels[0].patches[ip].patch_datas[name].x[:]
    y  = hier.patch_levels[0].patches[ip].patch_datas[name].y[:]
    nx = x.size
    ny = y.size
    qty = hier.patch_levels[0].patches[ip].patch_datas[name].dataset[:].reshape((nx,ny))
    return qty, x, y
def getData(path, time=None):
    times = get_times_from_h5(path+"EM_B.h5")
    times.max()
    if time is None:
        time = times.max()
    B = hierarchy_from(h5_filename = path+"EM_B.h5", time = "{:0.10f}".format(time))
    E = hierarchy_from(h5_filename = path+"EM_E.h5", time = "{:0.10f}".format(time))
    Ni= hierarchy_from(h5_filename = path+"ions_density.h5", time = "{:0.10f}".format(time))
    return B,E,Ni
def plot2D(hier,time, qty, vmin, vmax, **kwargs):
    fig,ax = plt.subplots(figsize=(10,5))
    for ilvl, plvl in hier.levels(time).items():
        for ip, patch in enumerate(plvl.patches):
            print("patch info", patch.id, patch.box, patch.box.shape)
            pdat = patch.patch_datas[qty]
            nbr_ghost = 5
            val = pdat.dataset[:]
            x = pdat.x
            y = pdat.y
            nx = x.size
            ny = y.size
            val2d = val.reshape((nx,ny))[nbr_ghost:-nbr_ghost, nbr_ghost:-nbr_ghost]
            x = x[nbr_ghost:-nbr_ghost]
            y = y[nbr_ghost:-nbr_ghost]
            dx,dy = pdat.layout.dl
            im = ax.pcolormesh(x, y, val2d.T, cmap="Spectral_r", vmin=vmin, vmax=vmax)
            r = Rectangle((patch.layout.box.lower[0]*dx, patch.layout.box.lower[1]*dy),
                          patch.layout.box.shape[0]*dx, patch.layout.box.shape[1]*dy,
                         fc="none", ec="k", alpha=0.4, lw=0.8)
            ax.add_patch(r)
            ax.text(x.min() + (x.max()-x.min())/5,
                    y.min() + (y.max()-y.min())/5, f"{patch.id}",
                   fontsize=10)
    ax.set_aspect('auto')
    ax.set_title(f" qty = {qty} t = {time}")
    fig.colorbar(im, ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if "filename" in kwargs:
        fig.savefig(kwargs["filename"])


def main():

    config()
    # startMPI()
    s = Simulator(gv.sim, post_advance=post_advance)
    s.initialize()
    if cpp.mpi_rank() == 0:
       plot2D(get_time(diag_outputs, 0), 0, "Ez", -1, 1, filename="thing.png")
    post_advance(0)
    s.run()


if __name__=="__main__":
    main()