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



def config(**kwargs):

    Simulation(
        smallest_patch_size=20,
        largest_patch_size=20,
        time_step_nbr=10000,
        final_time=10.,
        cells=(40, 40),
        dl=(0.2, 0.2),
        hyper_resistivity=0.001,
        resistivity=0.000,
        diag_options={"format": "phareh5",
                      "options": {"dir": kwargs["diagdir"],
                                  "mode":"overwrite"}}
    )

    def density(x, y):
        return 1.

    def bx(x, y):
        return kwargs.get("bx",1.)

    def by(x, y):
        return kwargs.get("by", 0.)

    def bz(x, y):
        return kwargs.get("bz",0.)

    def vx(x, y):
        return kwargs.get("vx",0.)

    def vy(x, y):
        return kwargs.get("vy",0.)

    def vz(x, y):
        return 0.

    def vthx(x, y):
        return 1.

    def vthy(x, y):
        return 1.

    def vthz(x, y):
        return 1.

    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz,
        "nbr_part_per_cell":100
    }

    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density,  **vvv}
    )

    ElectronModel(closure="isothermal", Te=kwargs.get("Te",0.))



    sim = ph.global_vars.sim
    dt = 10*sim.time_step
    nt = sim.final_time/dt+1
    timestamps = dt * np.arange(nt)
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

from pyphare.pharesee.hierarchy import get_times_from_h5
from pyphare.pharesee.hierarchy import hierarchy_from


def getData(path, time=None):
    times = get_times_from_h5(path+"EM_B.h5")
    times.max()
    if time is None:
        time = times.max()

    B = hierarchy_from(h5_filename = path+"EM_B.h5", time = "{:0.10f}".format(time))
    E = hierarchy_from(h5_filename = path+"EM_E.h5", time = "{:0.10f}".format(time))
    Ni= hierarchy_from(h5_filename = path+"ions_density.h5", time = "{:0.10f}".format(time))

    return B,E,Ni



def getEmeans(path):
    """
    returns the average electric field for each component over the whole domain
    for each time of the run at the given path
    """

    # this assumes all E,B,Ni have the same timestamps, which they do here..
    times = get_times_from_h5(path+"EM_B.h5")

    Emeans = {qty:np.zeros_like(times) for qty in ("Ex", "Ey", "Ez")}
    for it,t in enumerate(times):
        B,E, Ni = getData(path, t)
        ncells = {qty:[] for qty in ("Ex", "Ey", "Ez")}
        for component in ("Ex", "Ey", "Ez"):
            for ip, patch in enumerate(E.patch_levels[0].patches):
                pdat = patch.patch_datas[component]
                val = pdat.dataset[:]
                x = pdat.x
                nbr_ghost = 5
                y = pdat.y
                nx = x.size
                ny = y.size
                val2d = val.reshape((nx,ny))[nbr_ghost:-nbr_ghost, nbr_ghost:-nbr_ghost]
                Emeans[component][it]+=val2d.sum()
                ncells[component]+= [val2d.size]
            Emeans[component][it] /= sum(ncells[component])
    return Emeans

import glob

def makeFigs():

    runs = glob.glob("run*")

    expecteds = ((0,0,0),
                (0,0,0),
                (0,0,1),
                (-1,1,0),
                (0,1,0),
                (0,0,-1),
                )

    runid = 0
    for run, expected in zip(runs, expecteds):
        Emean = getEmeans(run)

        fig, axrr = plt.subplots(ncols=3, figsize=(9,4))
        for ax,component,exp in zip(axrr, ("Ex", "Ey", "Ez"), expected):
            ax.plot(np.arange(0,10,0.01), Emean[component], label="avg")
            ax.axhline(Emean[component].mean(), label="time avg")
            ax.axhline(exp ,ls="--", label="theory", color="r")
            ax.legend()
        fig.tight_layout()
        fig.savefig(f"Emean_{runid}.png")




def main():

    params={
        0:{},
        1:{"vx":1.}, # V//B
        2:{"vy":1.}, # Ez = 1
        3:{"vx":1.,"vy":1.,"by":1.,"bz":1.}, #Ex=-1, Ey=1, Ez=-1
        4:{"bz":1.,"bx":0.,"vx":1.}, # Ey = 1
        5:{"by":1.,"bx":0.,"vx":1.}, # Ez = -1
    }
    base_diag_dir = "run{:03d}"
    for irun, param in params.items():
        param.update({"diagdir":base_diag_dir.format(irun)})
        config(**param)
        simulator = Simulator(gv.sim).initialize().run()
        gv.sim = None



if __name__=="__main__":
    main()
