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
        final_time=50,
        time_step=0.0005,
        boundary_types="periodic",
        cells=165,
        dl=0.2,
        hyper_resistivity = 0.01,
        refinement_boxes={"L0": {"B0": [( 50, ), (110, )]},
                          "L1": {"B0": [(140, ), (180, )]} },
        diag_options={"format": "phareh5",
                      "options": {"dir": "ion_ion_beam1d", "mode": "overwrite"}}
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

    ElectronModel(closure="isothermal", Te=0.0)

    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time, 1.)

    for quantity in ["B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )





def yaebx(x, a, b):
    return a*np.exp(np.multiply(b, x))


def growth_b_right_hand(run_path):
    file = os.path.join(run_path, "EM_B.h5")
    times = get_times_from_h5(file)

    r = Run(run_path)
    first_mode = np.array([])

    from scipy.optimize import curve_fit

    for time in times:
        B_hier = r.GetB(time, merged=True, interp="linear")

        by_interpolator, xyz_finest = B_hier["By"]
        bz_interpolator, xyz_finest = B_hier["Bz"]

        # remove the last point so that "x" is periodic wo. last point = first point
        x = xyz_finest[0][:-1]

        by = by_interpolator(x)
        bz = by_interpolator(x)

        # get the mode 1, as it is the most unstable in a box of length 33
        mode1 = np.absolute(np.fft.fft(by-1j*bz)[1])
        first_mode = np.append(first_mode, mode1)

    popt, pcov = curve_fit(yaebx, times, first_mode, p0=[0.08, 0.09])

    # now the signal is stripped from its exponential part
    damped_mode=first_mode*yaebx(times, 1/popt[0], -popt[1])

    # find the omega for which "damped_mode" is the largest :
    # this term is twice the one it should be because "mode1" resulting from
    # an absolute value, this (cosx)^2 = cos(2x) then appears at the 2nd
    # harmonoic (hence the factor 0.5 to get "omega")
    # the factor "+1" is because we remove the DC component, so the value
    # given by argmax has also to miss this value
    omegas = np.fabs(np.fft.fft(damped_mode).real)
    omega = 0.5*(omegas[1:omegas.size//2].argmax()+1)*2*np.pi/times[-1]

    print(omegas[0:6])
    print(omegas[1:omegas.size//2].argmax())
    print(2*np.pi/times[-1])
    print(0.5*(omegas[1:omegas.size//2].argmax()+1)*2*np.pi/times[-1])

    return times, first_mode, popt[0], popt[1], damped_mode, omega


def main():
    from pybindlibs.cpp import mpi_rank

    config()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

    if mpi_rank() == 0:

        times, first_mode, ampl, gamma, damped_mode, omega \
                = growth_b_right_hand(os.path.join(os.curdir, "ion_ion_beam1d"))

        fig, (ax1, ax2) = plt.subplots(2, 1)

        ax1.set_title("Right Hand Resonant mode (Beam instability)")
        ax1.stem(times, first_mode, linefmt='-k', basefmt=' ', use_line_collection=True)
        ax1.plot(times, yaebx(times, ampl, gamma), color='r', linestyle='-', marker='')
        ax1.text(0.04, 0.80, "From Gary et al., 1985 (ApJ : 10.1086/162797)", transform=ax1.transAxes)
        ax1.set_ylabel("Most unstable mode")
        ax1.set_title("Right Hand Resonant mode (Beam instability)")
        ax1.text(0.30, 0.50, "gamma = {:5.3f}... expected 0.09".format(gamma), transform=ax1.transAxes)

        ax2.plot(times, damped_mode, color='g', linestyle='', marker='o')
        ax2.set_xlabel("Time")
        ax2.set_ylabel("Real mode")
        ax2.text(0.48, 0.30, "~ 3 periods until t=50", transform=ax2.transAxes)
        ax2.text(0.40, 0.20, "omega (real) = {:5.3f}... expected 0.19".format(omega), transform=ax2.transAxes)

        fig.savefig("ion_ion_beam1d.png")

        # compare with the values given gary et al. 1985
        assert np.fabs(gamma-0.09) < 2e-2
        #assert np.fabs(omega-0.19) < 8e-2


if __name__=="__main__":
    main()

