#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharein import global_vars as gv

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


def omega(k, p):
    k2 = k*k
    return 0.5*k2*(np.sqrt(1+4/k2)+p)




def setOfModes(polarization, modes, b_amplitudes, theta, seed):

    Simulation(
        smallest_patch_size=10,
        largest_patch_size=10,
        time_step_nbr=800000,   # 40000
        final_time=200.,        # 40.
        boundary_types=("periodic", "periodic"),
        cells=(400, 20),
        dl=(0.1, 0.1),
        diag_options={"format": "phareh5",
                      "options": {"dir": "setOfModes2d",
                                  "mode":"overwrite"}}
    )

    assert(len(modes) == len(b_amplitudes))

    # list of wave_numbers for the given box
    from pyphare.pharein.global_vars import sim
    L = sim.simulation_domain()[0]
    wave_numbers = [2*np.pi*m/L for m in modes]

    wave_num_x = [k*np.cos(theta) for k in wave_numbers]
    wave_num_y = [k*np.sin(theta) for k in wave_numbers]


    def density(x, y):
        # no density fluctuations as whistler and AIC are not compressional
        return 1.

    def bx(x, y):
        modes = np.cos(theta) # DC magnetic field of amplitude 1
        for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes):
            modes -= b*np.cos(kx*x+ky*y)*np.sin(theta)
        return modes

    def by(x, y):
        modes = np.sin(theta) # DC magnetic field of amplitude 1
        for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes):
            modes += b*np.cos(kx*x+ky*y)*np.cos(theta)
        return modes

    def bz(x, y):
        modes = 0.0
        for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes):
            modes += b*np.sin(kx*x+ky*y)*polarization
        return modes

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def vthx(x, y):
        return 0.01

    def vthy(x, y):
        return 0.01

    def vthz(x, y):
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
def get_all_w(run_path, wave_numbers, polarization, theta):
    file = os.path.join(run_path, "EM_B.h5")
    times = get_times_from_h5(file)

    nm = len(wave_numbers)
    print('number of modes : {}'.format(nm))

    r = Run(run_path)
    blz = np.array([])

    for time in times:
        interp_bx, xy = r.GetB(time, merged=True, interp='bilinear')['Bx']
        interp_by, xy = r.GetB(time, merged=True, interp='bilinear')['By']
        interp_bz, xy = r.GetB(time, merged=True, interp='bilinear')['Bz']

        X, Y = (xy[0], xy[0]*np.tan(theta))

        # the 1d space sampling is then of size nx
        # the (l, t) image will then be of size (nx, nt)
        bx = interp_bx(X, Y)
        by = interp_by(X, Y)
        bz = interp_bz(X, Y)

        # the "l" direction is then at +pi/2 from the theta, normal to B0
        bl = by*np.cos(theta)-bx*np.sin(theta)

        # polarization = +1 for R mode, -1 for L mode
        blz = np.concatenate((blz, bl+polarization*1j*bz))

    nx = xy[0].shape[0]
    nt = times.shape[0]
    blz = np.reshape(blz, (nt, nx))

    BLZ = np.absolute(np.fft.fft2(blz)[:(nt+1)//2, :(nx+1)//2])
    BLZ_4_all_W = np.sum(BLZ, axis=0)

    idx = np.argsort(BLZ_4_all_W)
    kmodes = idx[-nm:]

    #wmodes = []
    wmodes = np.array([])
    for i in range(nm):
        #wmodes.append(np.argmax(BLZ[1:,idx[-nm+i]]))
        wmodes = np.append(wmodes, np.argmax(BLZ[1:,idx[-nm+i]]))

    idx = np.argsort(kmodes)

    #print(kmodes, wmodes)
    print(kmodes[idx], wmodes[idx])

    #return kmodes, np.array(wmodes), BLZ
    return kmodes[idx], wmodes[idx], BLZ



def main():

    # angle of the oblique mode (in radians) : has to be arctan2(L[1], L[0])
    theta = np.arctan2(10, 200)

    # list of modes : m = 1 is for 1 wavelength in the whole domain
    modes = [4, 8, 16, 32, 64]

    # lists of amplitudes of the magnetic field amplitudes
    b_amplitudes = [0.002, 0.002, 0.002, 0.002, 0.002]

    # polarization : -1 for L mode
    wave_nums, b1 = setOfModes(-1, modes, b_amplitudes, theta, cpp.mpi_rank()+1)
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

    from pybindlibs.cpp import mpi_rank
    from matplotlib import rc

    if mpi_rank() == 0:
        sim = ph.global_vars.sim

        # hmmm... because theta has to be set before simulator !
        np.testing.assert_allclose(np.tan(theta) , sim.simulation_domain()[1]/sim.simulation_domain()[0], atol=1e-15)

        L = sim.simulation_domain()[0]
        T = sim.final_time

        #for the left mode
        ki, wi, blz = get_all_w(os.path.join(os.curdir, "setOfModes2d"), wave_nums, -1, theta)

        np.save('left2d.npy', blz)

        k_numL = 2*np.pi*ki*np.cos(theta)/L
        w_numL = 2*np.pi*wi/T

    #because the simulation is already set
    ph.global_vars.sim = None

    # list of modes : m = 1 is for 1 wavelength in the whole domain
    modes = [4, 8, 16, 32, 64]

    # lists of amplitudes of the magnetic field amplitudes
    b_amplitudes = [0.005, 0.005, 0.005, 0.005, 0.005]

    # polarization : +1 for R mode
    wave_nums, b1 = setOfModes(+1, modes, b_amplitudes, theta, cpp.mpi_rank()+1)
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

    if mpi_rank() == 0:
        #sim = ph.global_vars.sim

        L = sim.simulation_domain()[0]
        T = sim.final_time

        #for the riht mode
        ki, wi, blz = get_all_w(os.path.join(os.curdir, "setOfModes2d"), wave_nums, +1, theta)

        np.save('right2d.npy', blz)

        k_numR = 2*np.pi*ki*np.cos(theta)/L
        w_numR = 2*np.pi*wi/T

        rc('text', usetex = True)

        fig, ax = plt.subplots(figsize=(4,3), nrows=1)

        k_the = np.arange(0.2, 20, 0.001)
        w_thR = omega(k_the, +1)
        w_thL = omega(k_the, -1)

        #ax.imshow(zobi, origin='lower', cmap='viridis_r')
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
        fig.savefig("dispersion2d.pdf", dpi=200)

        w_theR = omega(k_numR, +1)
        w_theL = omega(k_numL, -1)

        errorL = 100*np.fabs(w_numL-w_theL)/w_theL
        errorR = 100*np.fabs(w_numR-w_theR)/w_theR

        with open('dispersion2d.txt', 'w') as f:
            print(*('error Left ... k = {:.4f}   w_the = {:.4f}   w_num = {:.4f}   err = {:.4f}'.\
                    format(k, W, w, e) for (k, W, w, e) in zip(k_numL, w_theL, w_numL, errorL)), sep="\n", file=f)
            print(*('error Right... k = {:.4f}   w_the = {:.4f}   w_num = {:.4f}   err = {:.4f}'.\
                    format(k, W, w, e) for (k, W, w, e) in zip(k_numR, w_theR, w_numR, errorR)), sep="\n", file=f)



if __name__=="__main__":
    main()



