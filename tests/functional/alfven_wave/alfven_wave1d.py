#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
from pyphare.pharesee.hierarchy import get_times_from_h5
from tests.diagnostic import all_timestamps
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import finest_field


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')



####################################################################
#
#                     Simulation configuration
#
####################################################################
def config():

    # configure the simulation

    Simulation(
        smallest_patch_size=50,
        largest_patch_size=50,
        time_step_nbr=100000,        # number of time steps (not specified if time_step and final_time provided)
        final_time=1000,             # simulation final time (not specified if time_step and time_step_nbr provided)
        boundary_types="periodic", # boundary condition, string or tuple, length == len(cell) == len(dl)
        cells=1000,                # integer or tuple length == dimension
        dl=1,                  # mesh size of the root level, float or tuple
        hyper_resistivity = 0.001,
        refinement_boxes={"L0": {"B0": [(450, ), (550, )]}},
        diag_options={"format": "phareh5",
                      "options": {"dir": ".","mode":"overwrite"}}
    )


    def density(x):
        return 1.


    def by(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.01*np.cos(2*np.pi*x/L[0])


    def bz(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.01*np.sin(2*np.pi*x/L[0])


    def bx(x):
        return 1.


    def vx(x):
        return 0.


    def vy(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.01*np.cos(2*np.pi*x/L[0])

    def vz(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.01*np.sin(2*np.pi*x/L[0])


    def vthx(x):
        return 0.01


    def vthy(x):
        return 0.01


    def vthz(x):
        return 0.01


    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthx, "vthy": vthy, "vthz": vthz
    }

    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, **vvv}
    )

    ElectronModel(closure="isothermal", Te=0.0)

    timestamps = all_timestamps(gv.sim)

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




####################################################################
#                      post processing
####################################################################

def wave(x, a0, k, phi):
    return a0*np.cos(k*x + phi)


def phase_speed(run_path, ampl, xmax):
    from scipy.signal import medfilt
    from scipy.optimize import curve_fit
    import os
    time = get_times_from_h5(os.path.join(run_path,"EM_B.h5"))
    r = Run(run_path)
    phase = np.zeros_like(time)
    amplitude = np.zeros_like(time)
    wave_vec = np.zeros_like(time)

    for it, t in enumerate(time):
        B = r.GetB(t)
        by, xby = finest_field(B, "By")
        a, k, phi = curve_fit(wave, xby, by,p0=(ampl, 2*np.pi/xmax, 0))[0]
        phase[it] = phi
        amplitude[it] = a
        wave_vec[it]  = k

    vphi = medfilt(np.gradient(phase, time)/wave_vec, kernel_size=7)
    return vphi, time, phase, amplitude, wave_vec



def main():
    from pybindlibs.cpp import mpi_rank
    from pyphare.pharesee.run import Run
    from pyphare.pharesee.hierarchy import finest_field
    from pyphare.cpp import cpp_lib
    mpi_rank = cpp_lib().mpi_rank


    config()
    Simulator(gv.sim).initialize().run()

    if mpi_rank() == 0:

        vphi, t, phi, a, k = phase_speed(".", 0.01, 1000)

        r = Run(".")
        t = get_times_from_h5("EM_B.h5")
        fig, ax = plt.subplots(figsize=(9,5), nrows=1)

        B = r.GetB(t[int(len(t)/2)])
        by, xby = finest_field(B, "By")
        ax.plot(xby, by, label="t = 500", alpha=0.6)

        sorted_patches=sorted(B.patch_levels[1].patches,
                              key=lambda p:p.box.lower[0])

        x0 = sorted_patches[0].patch_datas["By"].x[0]
        x1 = sorted_patches[-1].patch_datas["By"].x[-1]

        B = r.GetB(t[-1])
        by, xby = finest_field(B, "By")
        ax.plot(xby, by, label="t = 1000", alpha=0.6)
        ax.plot(xby, wave(xby, 0.01, 2*np.pi/1000., 2*np.pi/1000*500),
                color="k", ls="--", label="T=500 (theory)")

        B = r.GetB(t[0])
        by, xby = finest_field(B, "By")
        ax.plot(xby, by, label="t = 0", color="k")

        ax.set_xlabel("x")
        ax.set_ylabel(r"$B_y$")
        ax.legend(ncol=4, loc="upper center")
        ax.set_ylim((-0.012, 0.013))
        ax.set_title(r"$V_\phi = {:6.4f}$".format(vphi.mean()))

        ax.axvspan(x0, x1, alpha=0.2)
        fig.tight_layout()

        fig.savefig("alfven_wave.png", dpi=200)

        assert np.mean(np.abs(vphi-1) < 5e-2)


if __name__=="__main__":
    main()
