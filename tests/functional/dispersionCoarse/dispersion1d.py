#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharein import global_vars as gv

#from pyphare.pharesee.hierarchy import finest_field
import os
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import get_times_from_h5

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')



from pyphare.cpp import cpp_lib
cpp = cpp_lib()
startMPI()


def fromNoise():

    # in this configuration there are no prescribed waves
    # and only eigen modes existing in the simulation noise
    # will be visible

    Simulation(
        smallest_patch_size=20,
        largest_patch_size=20,

        # the following time step number
        # and final time mean that the
        # smallest frequency will be 2/100
        # and the largest 2/dt  = 2e3
        time_step_nbr=100000,
        final_time=100.,

        boundary_types="periodic",

        # smallest wavelength will be 2*0.2=0.4
        # and largest 50
        cells=500,
        dl=0.2,
        diag_options={"format": "phareh5",
                      "options": {"dir": "fromNoise1d",
                                  "mode":"overwrite"}}
    )

    def density(x):
        return 1.

    def by(x):
        return 0.

    def bz(x):
        return 0.

    def bx(x):
        return 1.

    def vx(x):
        return 0.

    def vy(x):
        return 0.

    def vz(x):
        return 0.

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

    ElectronModel(closure="isothermal", Te=0.)


    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time +sim.time_step, sim.time_step)

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )




def omega(k, p):
    k2 = k*k
    return 0.5*k2*(np.sqrt(1+4/k2)+p)




def setOfModes(polarization, modes, b_amplitudes, seed):

    Simulation(
        smallest_patch_size=20,
        largest_patch_size=50,
        time_step_nbr=300000,
        final_time=200.,
        boundary_types="periodic",
        cells=2000,
        dl=0.2,
        diag_options={"format": "phareh5",
                      "options": {"dir": "setOfModes1d",
                                  "mode":"overwrite"}}
    )

    assert(len(modes) == len(b_amplitudes))

    # list of wave_numbers for the given box
    from pyphare.pharein.global_vars import sim
    L = sim.simulation_domain()[0]
    wave_numbers = [2*np.pi*m/L for m in modes]

    # using faraday : v1 = -w b1 / (k . B0)
    #v_amplitudes = [-b*omega(k, p)/k for (k, b, p) in zip(wave_numbers, b_amplitudes, polarizations)]


    def density(x):
        # no density fluctuations as whistler and AIC are not compressional
        return 1.

    def by(x):
        modes = 0.0
        for (k, b) in zip(wave_numbers, b_amplitudes):
            modes += b*np.cos(k*x)
        return modes

    def bz(x):
        modes = 0.0
        for (k, b) in zip(wave_numbers, b_amplitudes):
            modes += b*np.sin(k*x)*polarization
        return modes

    def bx(x):
        return 1.

    def vx(x):
        return 0.

    def vy(x):
        #modes = 0.0
        #for (k, v, f) in zip(wave_numbers, v_amplitudes, phases):
        #    modes += v*np.cos(k*x+f)
        #return modes
        return 0.0

    def vz(x):
        #modes = 0.0
        #for (k, v, f) in zip(wave_numbers, b_amplitudes, phases):
        #    modes += v*np.sin(k*x+f)
        #return modes
        return 0.0

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
        main={"charge": 1, "density": density, **vvv, "init":{"seed": seed}}
    )

    ElectronModel(closure="isothermal", Te=0.)

    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time+sim.time_step, 10*sim.time_step)


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

    return wave_numbers, b_amplitudes



# ___ post-processing functions
def get_all_w(run_path, wave_numbers, polarization):
    file = os.path.join(run_path, "EM_B.h5")
    times = get_times_from_h5(file)

    nm = len(wave_numbers)
    print('number of modes : {}'.format(nm))

    r = Run(run_path)
    byz = np.array([])

    for time in times:

        interp_by, x = r.GetB(time, merged=True, interp='nearest')['By']
        interp_bz, x = r.GetB(time, merged=True, interp='nearest')['Bz']
        by = interp_by(x[0])
        bz = interp_bz(x[0])

        # polarization = +1 for R mode, -1 for L mode
        byz = np.concatenate((byz, by+polarization*1j*bz))

    nx = x[0].shape[0]
    nt = times.shape[0]
    byz = np.reshape(byz, (nt, nx))

    BYZ = np.absolute(np.fft.fft2(byz)[:(nt+1)//2, :(nx+1)//2])
    BYZ_4_all_W = np.sum(BYZ, axis=0)

    idx = np.argsort(BYZ_4_all_W)
    kmodes = idx[-nm:]

    wmodes = np.array([])
    for i in range(nm):
        wmodes = np.append(wmodes, np.argmax(BYZ[:,idx[-nm+i]]))
        #wmodes.append(np.argmax(BYZ[:,idx[-nm+i]]))

    idx = np.argsort(kmodes)

    print(kmodes[idx], wmodes[idx])

    return kmodes[idx], wmodes[idx], BYZ



def main():
    # list of modes : m = 1 is for 1 wavelength in the whole domain
    modes = [4, 8, 16, 32, 64, 128, 256, 512]

    # lists of amplitudes of the magnetic field amplitudes
    b_amplitudes = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]

    # polarization : -1 for L mode
    wave_nums, b1 = setOfModes(-1, modes, b_amplitudes, cpp.mpi_rank()+1)
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

    from pybindlibs.cpp import mpi_rank
    from matplotlib import rc

    if mpi_rank() == 0:
        sim = ph.global_vars.sim

        L = sim.simulation_domain()[0]
        T = sim.final_time

        #for the left mode
        ki, wi, byz = get_all_w(os.path.join(os.curdir, "setOfModes1d"), wave_nums, -1)

        np.save('left1d.npy', byz)

        k_numL = 2*np.pi*ki/L
        w_numL = 2*np.pi*wi/T

    ph.global_vars.sim = None

    # list of modes : m = 1 is for 1 wavelength in the whole domain
    modes = [4, 8, 16, 32, 64, 128, 256, 512]

    # lists of amplitudes of the magnetic field amplitudes
    b_amplitudes = [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005]

    # polarization : +1 for R mode
    wave_nums, b1 = setOfModes(+1, modes, b_amplitudes, cpp.mpi_rank()+1)
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

    if mpi_rank() == 0:
        #sim = ph.global_vars.sim

        L = sim.simulation_domain()[0]
        T = sim.final_time

        #for the riht mode
        ki, wi, byz = get_all_w(os.path.join(os.curdir, "setOfModes1d"), wave_nums, +1)

        np.save('right1d.npy', byz)

        k_numR = 2*np.pi*ki/L
        w_numR = 2*np.pi*wi/T

        rc('text', usetex = True)

        fig, ax = plt.subplots(figsize=(4,3), nrows=1)

        k_the = np.arange(0.04, 10, 0.001)
        w_thR = omega(k_the, +1)
        w_thL = omega(k_the, -1)

        ax.plot(k_the, w_thR, '-k')
        ax.plot(k_the, w_thL, '-k')
        ax.plot(k_numR, w_numR, 'b+', label='$R \ mode$', markersize=8)
        ax.plot(k_numL, w_numL, 'rx', label='$L \ mode$', markersize=8)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('$k_{\parallel} c / \omega_p$')
        ax.set_ylabel('$\omega / \Omega_p$')

        ax.legend(loc='upper left', frameon=False)

        fig.tight_layout()
        fig.savefig("dispersion1d.pdf", dpi=200)

        w_theR = omega(k_numR, +1)
        w_theL = omega(k_numL, -1)

        errorL = 100*np.fabs(w_numL-w_theL)/w_theL
        errorR = 100*np.fabs(w_numR-w_theR)/w_theR

        with open('dispersion1d.txt', 'w') as f:
            print(*('error Left ... k = {:.4f}   w_the = {:.4f}   w_num = {:.4f}   err = {:.4f}'.\
                    format(k, W, w, e) for (k, W, w, e) in zip(k_numL, w_theL, w_numL, errorL)), sep="\n", file=f)
            print(*('error Right... k = {:.4f}   w_the = {:.4f}   w_num = {:.4f}   err = {:.4f}'.\
                    format(k, W, w, e) for (k, W, w, e) in zip(k_numR, w_theR, w_numR, errorR)), sep="\n", file=f)

        targetL = np.array([ 3.,  6.,  1.,  4.,  1.,  2.,  4.,  1.])
        targetR = np.array([ 3.,  6.,  1.,  3.,  0.,  3.,  6., 20.])

        np.testing.assert_allclose(errorL, targetL, rtol=1e-2, atol=8)
        np.testing.assert_allclose(errorR, targetR, rtol=1e-2, atol=8)

        #assert errorL.max() < 10.0
        #assert errorR.max() < 30.0


if __name__=="__main__":
    main()

