#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv
from pyphare.pharesee.hppierarchy import get_times_from_h5
from tests.diagnostic import all_timestamps
from pyphare.pharesee.run import Run
import os
from pyphare.pharein.global_vars import sim



import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')



def config():

    # most unstable mode at k=0.19, that is lambda = 33
    # hence the length of the box is 33, and the fourier mode will be 1

    Simulation(
        smallest_patch_size=20,
        largest_patch_size=60,
        time_step=0.0005, 
        final_time=40,
        boundary_types="periodic",
        cells=165,
        dl=0.2,
        hyper_resistivity = 0.01,
        refinement_boxes={"L0": {"B0": [( 50, ), (110, )]},
                          "L1": {"B0": [(140, ), (180, )]} },
        diag_options={"format": "phareh5",
                      "options": {"dir": "ion_ion_beam", "mode": "overwrite"}}
    )

    def densityMain(x):
        return 1.
   
    def densityBeam(x):
        return .01
   
    def bx(x):
        return 1.
   
    def by(x):
        return 0.
   
    def bz(x):
        return 0.
   
    def vB(x):
        return 5.
   
    def v0(x):
        return 0.
   
    def vth(x):
        return np.sqrt(0.1)
   
   
    vMain = {
        "vbulkx": v0, "vbulky": v0, "vbulkz": v0,
        "vthx": vth, "vthy": vth, "vthz": vth
    }
   
   
    vBulk = {
        "vbulkx": vB, "vbulky": v0, "vbulkz": v0,
        "vthx": vth, "vthy": vth, "vthz": vth
    }
   
   
    MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        main={"charge": 1, "density": densityMain, **vMain},
        beam={"charge": 1, "density": densityBeam, **vBulk}
    )

    ElectronModel(closure="isothermal", Te=0.1)

    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time, 1.)

    for quantity in ["B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )





def getMode(t, m):
    # return the mode 'm' of the FFT : we the use m=1 for 1 wavelength in the whole simulation domain
    return np.absolute(np.fft.fft(t)[m])

def yaebx(x, a, b):
    return a*np.exp(np.multiply(b, x))

def growth_b_right_hand(run_path):
    file = os.path.join(run_path, "EM_B.h5")
    times = get_times_from_h5(file)

    r = Run(run_path)
    byz = np.array([])

    from scipy.optimize import curve_fit

    for time in times:
        interp_by, x = r.GetB(time, merged=True, interp='nearest')['By']
        interp_bz, x = r.GetB(time, merged=True, interp='nearest')['Bz']

        by = interp_by(x[0])
        bz = interp_bz(x[0])

        mm = getMode(by-1j*bz, 1) # last arg is the mode

        byz = np.append(byz, mm)

    popt, pcov = curve_fit(yaebx, times, byz, p0=[0.1, 0.1]) # last args are initial guess

    return times, byz, popt, pcov


def main():
    from pybindlibs.cpp import mpi_rank
    from matplotlib import rc

    config()
    Simulator(gv.sim).initialize().run()

    if mpi_rank() == 0:

        times, byz, popt, pcov = growth_b_right_hand(os.path.join(os.curdir, "ion_ion_beam"))

        rc('text', usetex = True)

        fig, ax = plt.subplots(figsize=(6,4), nrows=1)

        ax.plot(times, byz, color='k', linestyle=' ', marker='o')
        ax.plot(times, yaebx(times, popt[0], popt[1]), color='red')

        ax.set_xlabel("$t \Omega_p$")
        ax.set_ylabel("$B_y - \imath B_z \ $ (m=1)")
        ax.set_title("Ion/Ion beam instability : Right-Hand-Resonand mode")
        ax.text(30, 0, "$\gamma / \Omega_p$ = {:6.4f}".format(popt[1]))

        fig.tight_layout()
        fig.savefig("ion_ion_beam.pdf", dpi=200)

        cov = np.fabs(np.array(pcov).flatten())

        assert np.fabs(popt[1]-0.09) < 1e-2
        assert cov.max() < 1e-2


if __name__=="__main__":
    main()

