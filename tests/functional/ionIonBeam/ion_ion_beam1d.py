import os

import numpy as np
import pyphare.pharein as ph
import matplotlib.pyplot as plt

from pyphare.cpp import cpp_lib
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5
from pyphare.pharesee.run import Run

ph.NO_GUI()

cpp = cpp_lib()


def config():
    # most unstable mode at k=0.19, that is lambda = 33
    # hence the length of the box is 33, and the fourier mode will be 1

    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=60,
        final_time=100,
        time_step=0.001,
        boundary_types="periodic",
        cells=165,
        dl=0.2,
        hyper_resistivity=0.01,
        refinement="tagging",
        max_nbr_levels=3,
        # refinement_boxes={
        #     "L0": {"B0": [(50,), (110,)]},
        #     "L1": {"B0": [(140,), (180,)]},
        # },
        diag_options={
            "format": "phareh5",
            "options": {"dir": "ion_ion_beam1d", "mode": "overwrite"},
        },
    )

    def densityMain(x):
        return 1.0

    def densityBeam(x):
        return 0.01

    def bx(x):
        return 1.0

    def by(x):
        return 0.0

    def bz(x):
        return 0.0

    def vB(x):
        return 5.0

    def v0(x):
        return 0.0

    def vth(x):
        return np.sqrt(0.1)

    vMain = {
        "vbulkx": v0,
        "vbulky": v0,
        "vbulkz": v0,
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
    }

    vBulk = {
        "vbulkx": vB,
        "vbulky": v0,
        "vbulkz": v0,
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
    }

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        main={"charge": 1, "density": densityMain, **vMain},
        beam={"charge": 1, "density": densityBeam, **vBulk},
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    timestamps = np.arange(0, sim.final_time, 0.1)

    for quantity in ["B", "E"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

    for pop_name in ["main", "beam"]:
        ph.ParticleDiagnostics(
            quantity="domain", population_name=pop_name, write_timestamps=timestamps
        )
    return sim


def yaebx(x, a, b):
    return a * np.exp(np.multiply(b, x))


def growth_b_right_hand(run_path, time_offset):
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks

    file = os.path.join(run_path, "EM_B.h5")
    times = get_times_from_h5(file)
    dt = times[1] - times[0]
    r = Run(run_path)
    first_mode = np.array([])

    for time in times:
        B_hier = r.GetB(time, merged=True, interp="linear")

        by_interpolator, xyz_finest = B_hier["By"]
        bz_interpolator, xyz_finest = B_hier["Bz"]

        # remove the last point so that "x" is periodic wo. last point = first point
        x = xyz_finest[0][:-1]

        by = by_interpolator(x)
        bz = bz_interpolator(x)

        # get the mode 1, as it is the most unstable in a box of length 33
        mode1 = np.absolute(np.fft.fft(by - 1j * bz)[1])
        first_mode = np.append(first_mode, mode1)

    ioffset = int(time_offset / dt)
    imax = find_peaks(first_mode, width=ioffset)[0][0]

    # the curve_fit is performed from time index 0 to imax-ioffset as this offset prevent to use
    # the final part of the curve which is no more exponential as this is the end of the linear mode
    popt, pcov = curve_fit(
        yaebx, times[: imax - ioffset], first_mode[: imax - ioffset], p0=[0.08, 0.09]
    )

    # now the signal is stripped from its exponential part
    damped_mode = first_mode[: imax - ioffset] * yaebx(
        times[: imax - ioffset], 1 / popt[0], -popt[1]
    )

    # find the omega for which "damped_mode" is the largest :
    # this term is twice the one it should be because "mode1" resulting from
    # an absolute value, this (cosx)^2 = cos(2x) then appears at the 2nd
    # harmonoic (hence the factor 0.5 to get "omega")
    # the factor "+1" is because we remove the DC component, so the value
    # given by argmax has also to miss this value
    omegas = np.fabs(np.fft.fft(damped_mode).real)
    omega = (
        0.5
        * (omegas[1 : omegas.size // 2].argmax() + 1)
        * 2
        * np.pi
        / times[imax - 1 - ioffset]
    )

    return times, first_mode, popt[0], popt[1], damped_mode, omega


def main():
    from scipy.signal import find_peaks

    time_offset = 10.0
    # this is an offset so the exponential fit associated to the linear phase is not performed
    # until the time at which the B value gets the greater... but is rather before with an offset
    # given by time_offset : the exponential fit is performed on [0, times[imax] - time_offset]

    Simulator(config()).run()

    if cpp.mpi_rank() == 0:
        times, first_mode, ampl, gamma, damped_mode, omega = growth_b_right_hand(
            os.path.join(os.curdir, "ion_ion_beam1d"), time_offset
        )

        dt = times[1] - times[0]
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(
            4, 1, figsize=(6.5, 12), constrained_layout=True
        )
        imax = find_peaks(first_mode, width=int(10 / dt))[0][
            0
        ]  # the smallest width of the peak is 10

        ax1.set_title("Time evolution of the first right-hand circular mode amplitude")
        ax1.plot(
            times,
            first_mode,
            color="k",
            label="|$\hat{B}_y (m=1,t)-i\hat{B}_z (m=1,t)$|",
        )
        ax1.plot(
            times[:imax],
            yaebx(times[:imax], ampl, gamma),
            color="r",
            linestyle="-",
            label="$B_0. \exp(\gamma t), \ with\ \gamma =${:5.5f} (expected 0.09)".format(
                gamma
            ),
        )
        ax1.axvline(0, 0, yaebx(times[imax], ampl, gamma), color="red", linestyle="--")
        ax1.axvline(
            times[imax] - time_offset,
            0,
            yaebx(times[imax], ampl, gamma),
            color="red",
            linestyle="--",
        )
        ax1.legend()
        ax1.set_xlabel("t - Time")

        r = Run(os.path.join(os.curdir, "ion_ion_beam1d"))
        dt = times[1] - times[0]
        ax_t2, ax_t3, ax_t4 = ax2.twinx(), ax3.twinx(), ax4.twinx()
        vmin, vmax = -1.3, 5.6
        for t in [0, times[imax], times[-1]]:
            if t == 0:
                ax, ax_t = ax2, ax_t2
            elif t == times[-1]:
                ax, ax_t = ax4, ax_t4
            else:
                ax, ax_t = ax3, ax_t3

            ions = r.GetParticles(t, ["main", "beam"])

            ions.dist_plot(
                axis=("x", "Vx"),
                ax=ax,
                norm=0.4,
                finest=True,
                gaussian_filter_sigma=(1, 1),
                vmin=vmin,
                vmax=vmax,
                dv=0.05,
                title="t = {:.1f}".format(t),
                xlabel="",
                ylabel="",
            )
            ax.set_xlim((0, 33))
            ax.axvline(10, vmin, vmax, color="k")
            ax.axvline(22, vmin, vmax, color="k")
            ax.axvline(14, vmin, vmax, color="red")
            ax.axvline(18, vmin, vmax, color="red")

            E_hier = r.GetE(time=t, merged=True, interp="linear")
            ey_interpolator, xyz_finest = E_hier["Ey"]
            ax_t.plot(
                xyz_finest[0],
                ey_interpolator(xyz_finest[0]),
                linewidth=2,
                color="dimgray",
            )
            ax_t.set_ylim((-0.4, 0.4))

        ax_t3.set_ylabel("Ey(x)")
        ax3.set_ylabel("Vx - Velocity")
        ax4.set_xlabel("X - Position")
        fig.savefig("ion_ion_beam1d.png")

        assert np.fabs(gamma - 0.09) < 2e-2


if __name__ == "__main__":
    main()
