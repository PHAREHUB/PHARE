
import sys
import unittest
import numpy as np
from ddt import ddt, data, unpack
from tests.diagnostic import all_timestamps
import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest, parse_cli_args

from pyphare.cpp import cpp_lib
cpp = cpp_lib()

from pyphare.pharein.simulation import supported_dimensions

permutations = [
    [dim, res] for dim in supported_dimensions() for res in ["low", "high"]
]

configs = {
    1 : {
        "modes_left" : [4, 8, 16, 32, 64, 128, 256, 512],
        "modes_right" : [4, 8, 16, 32, 64, 128, 256, 512],
        "b_amplitudes_left" : [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
        "b_amplitudes_right" : [0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005],
        "low" :{
            "time_step_nbr": 300000,
            "cells"        : 2000,
            "dl"           : 0.2,
        },
        "high":{
            "time_step_nbr": 800000,
            "cells"        : 4000,
            "dl"           : 0.1,
        },
        "base":{ # this is overridden by resolution specific keys if they exist
            "smallest_patch_size":20,
            "largest_patch_size":50,
            "final_time":200.,
            "boundary_types":"periodic",
            "diag_options":{"format": "phareh5",
                          "options": {"dir": "setOfModes1d",
                                      "mode":"overwrite"}}
        },
    },
    2 : {
        "modes_left" : [4, 8, 16, 32, 64],
        "modes_right" : [4, 8, 16, 32, 64],
        "b_amplitudes_left" : [0.002, 0.002, 0.002, 0.002, 0.002],
        "b_amplitudes_right" : [0.005, 0.005, 0.005, 0.005, 0.005],
        "low" :{
            "time_step_nbr": 150000,
            "cells"        : (400, 20),
            "dl"           : (0.2, 0.2),
        },
        "high":{
            "time_step_nbr": 800000,
            "cells"        : (200, 10),
            "dl"           : (0.1, 0.1),
        },
        "base":{ # this is overridden by resolution specific keys if they exist
            "smallest_patch_size":10,
            "largest_patch_size":10,
            "final_time":200.,
            "boundary_types":("periodic", "periodic"),
            "diag_options":{"format": "phareh5",
                          "options": {"dir": "setOfModes2d",
                                      "mode":"overwrite"}}
        },
    },
}

def get_theta():
    """
    define the angle of the DC magnetic field so that it is exactly
    along the diagonal of the 2D box
    """
    from pyphare.pharein.global_vars import sim
    L = sim.cells
    return np.arctan2(L[1], L[0])


def setup(dim, resolution, modes, b_amplitudes, polarization, seed=12345):
    assert(len(modes) == len(b_amplitudes))

    ph.global_vars.sim = None

    config = configs[dim]["base"].copy()
    config.update(configs[dim][resolution])
    sim = ph.Simulation(**config)

    # list of wave_numbers for the given box
    from pyphare.pharein.global_vars import sim
    L = sim.simulation_domain()
    wave_numbers = [2*np.pi*m/L[0] for m in modes]

    if dim == 2:
        theta = get_theta()
        wave_num_x = [k*np.cos(theta) for k in wave_numbers]
        wave_num_y = [k*np.sin(theta) for k in wave_numbers]


    all_funcs = {
        "base" :{ # this is overridden by dimension specific keys if they exist
            "density": lambda *xyz: 1,
            "vthx": lambda *xyz: .01,
            "vthy": lambda *xyz: .01,
            "vthz": lambda *xyz: .01,
            "vbulkx": lambda *xyz: 0,
            "vbulky": lambda *xyz: 0,
            "vbulkz": lambda *xyz: 0,
        },
        1: {
            "bx": lambda *xyz: 1,
            "by": lambda x: np.sum([b*np.cos(k*x) for (k, b) in zip(wave_numbers, b_amplitudes)]),
            "bz": lambda x: np.sum([b*np.sin(k*x)*polarization for (k, b) in zip(wave_numbers, b_amplitudes)]),
        },
        2: { # could be moved to a function see, https://github.com/PHAREHUB/PHARE/blob/master/tests/simulator/__init__.py#L93 # rm comment
            "bx": lambda x, y: np.cos(theta) - np.sum([b*np.cos(kx*x+ky*y)*np.sin(theta) for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes)]),
            "by": lambda x, y: np.sin(theta) + np.sum([b*np.cos(kx*x+ky*y)*np.cos(theta) for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes)]),
            "bz": lambda x, y: np.sum([b*np.sin(kx*x+ky*y)*polarization for (kx, ky, b) in zip(wave_num_x, wave_num_y, b_amplitudes)]),
        },
    }

    funcs = all_funcs["base"].copy()
    funcs.update(all_funcs[dim])

    vvv = {
        "vbulkx": funcs["vbulkx"], "vbulky": funcs["vbulky"], "vbulkz": funcs["vbulkz"],
        "vthx": funcs["vthx"], "vthy": funcs["vthy"], "vthz": funcs["vthz"]
    }
    ph.MaxwellianFluidModel(
        bx=funcs["bx"], by=funcs["by"], bz=funcs["bz"],
        main={"charge": 1, "density": funcs["density"], **vvv, "init":{"seed": seed}}
    )
    ph.ElectronModel(closure="isothermal", Te=0.)

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

    return sim, wave_numbers, b_amplitudes





def get_all_w_1d(run_path, wave_numbers, polarization):
    file = os.path.join(run_path, "EM_B.h5")
    times = get_times_from_h5(file)

    nm = len(wave_numbers)

    r = Run(run_path)
    byz = np.array([])

    for time in times:
        B_hier = r.GetB(time, merged=True, interp="linear")

        by_interpolator, xyz_finest = B_hier["By"]
        bz_interpolator, xyz_finest = B_hier["Bz"]

        # remove the last point so that "x" is periodic wo. last point = first point
        x = xyz_finest[0][:-1]

        by = by_interpolator(x)
        bz = bz_interpolator(x)

        byz = np.concatenate((byz, by+polarization*1j*bz))

    ### nx = x[0].shape[0]
    nx = x.shape[0]
    nt = times.shape[0]
    byz = np.reshape(byz, (nt, nx))

    BYZ = np.absolute(np.fft.fft2(byz)[:(nt+1)//2, :(nx+1)//2])
    # for each k modes, the summation is over all angular frequencies
    # this is because one can have some blur on the adjacent modes
    # so one can not identify the appropriate w mode
    BYZ_4_all_W = np.sum(BYZ, axis=0)

    idx = np.argsort(BYZ_4_all_W)
    kmodes = idx[-nm:]

    wmodes = np.array([])
    for i in range(nm):
        wmodes = np.append(wmodes, np.argmax(BYZ[:,idx[-nm+i]]))

    idx = np.argsort(kmodes)

    return kmodes[idx], wmodes[idx]



def get_all_w_2d(run_path, wave_numbers, polarization, theta):
    file = os.path.join(run_path, "EM_B.h5")
    times = get_times_from_h5(file)

    nm = len(wave_numbers)

    r = Run(run_path)
    blz = np.array([])

    for time in times:
        B_hier = r.GetB(time, merged=True, interp="linear")

        bx_interpolator, xyz_finest = B_hier["Bx"]
        by_interpolator, xyz_finest = B_hier["By"]
        bz_interpolator, xyz_finest = B_hier["Bz"]

        # remove the last point so that "x" is periodic wo. last point = first point
        X, Y = (xyz_finest[0][:-1], xyz_finest[0][:-1]*np.tan(theta))

        # the 1d space sampling is then of size nx
        # the (l, t) image will then be of size (nx, nt)
        bx = bx_interpolator(X, Y)
        by = by_interpolator(X, Y)
        bz = bz_interpolator(X, Y)

        # the "l" direction is then at +pi/2 from the theta, normal to B0
        bl = by*np.cos(theta)-bx*np.sin(theta)

        # polarization = +1 for R mode, -1 for L mode
        blz = np.concatenate((blz, bl+polarization*1j*bz))

    nx = X.shape[0]
    nt = times.shape[0]
    blz = np.reshape(blz, (nt, nx))

    BLZ = np.absolute(np.fft.fft2(blz)[:(nt+1)//2, :(nx+1)//2])
    # for each k modes, the summation is over all angular frequencies
    # this is because one can have some blur on the adjacent modes
    # so one can not identify the appropriate w mode
    BLZ_4_all_W = np.sum(BLZ, axis=0)

    idx = np.argsort(BLZ_4_all_W)
    kmodes = idx[-nm:]

    wmodes = np.array([])
    for i in range(nm):
        wmodes = np.append(wmodes, np.argmax(BLZ[1:,idx[-nm+i]]))

    idx = np.argsort(kmodes)

    return kmodes[idx], wmodes[idx]



def post_sim_left(sim, **kwargs):
    if cpp.mpi_rank() == 0 and dim == 2:
        theta = get_theta()
    if cpp.mpi_rank() == 0:

        sim = ph.global_vars.sim

        L = sim.simulation_domain()[0]
        T = sim.final_time

        ki, wi = get_all_w_1d(os.path.join(os.curdir, "setOfModes1d"), wave_nums, -1)

        k_num_left_1d = 2*np.pi*ki/L
        w_num_left_1d = 2*np.pi*wi/T


        np.testing.assert_allclose(np.tan(theta) , sim.simulation_domain()[1]/sim.simulation_domain()[0], atol=1e-15)

        L = sim.simulation_domain()[0]
        T = sim.final_time

        ki, wi = get_all_w_2d(os.path.join(os.curdir, "setOfModes2d"), wave_nums, -1, theta)

        k_num_left_2d = 2*np.pi*ki*np.cos(theta)/L
        w_num_left_2d = 2*np.pi*wi/T

    cpp.mpi_barrier() # KEEP THIS!

    return k_num_left_1d, w_num_left_1d, k_num_left_2d, w_num_left_2d



def post_sim_right(sim, **kwargs):
    if cpp.mpi_rank() == 0 and dim == 2:
        theta = get_theta()
    if cpp.mpi_rank() == 0:
        # do things

        ### ........................... TODO

        pass
    cpp.mpi_barrier() # KEEP THIS!



def omega(k, p):
    k2 = k*k
    return 0.5*k2*(np.sqrt(1+4/k2)+p)



def dispersion(dim, resolution):
    seed = cpp.mpi_rank()+1
    config = configs[dim][resolution]

    # we need to set, run and post process 1 run for the left mode
    modes = configs[dim]["modes_left"]
    b_amplitudes = configs[dim]["b_amplitudes_left"]
    polarization = -1
    sim, wave_nums, b1 = setup(dim, resolution, modes, b_amplitudes, polarization, seed)
    Simulator(sim).run().reset()
    k_num_left_1d, w_num_left_1d, k_num_left_2d, w_num_left_2d = post_sim_left(sim)

    # ... and set, run and post process 1 run for the right mode
    modes = configs[dim]["modes_right"]
    b_amplitudes = configs[dim]["b_amplitudes_right"]
    polarization = 1
    sim, wave_nums, b1 = setup(dim, resolution, modes, b_amplitudes, polarization, seed)
    Simulator(sim).run().reset()
    k_num_right_1d, w_num_right_1d, k_num_right_2d, w_num_right_2d = post_sim_right(dim)


    # plot the dispersion function
    k_analytic = np.arange(0.2, 20, 0.001)
    w_analytic_right = omega(k_the, +1)
    w_analytic_left = omega(k_the, -1)

    fig, ax = plt.subplots(figsize=(4,3), nrows=1)

    ax.plot(k_analytic, w_analytic_right, '-k')
    ax.plot(k_analytic, w_analytic_left, '-k')

    ax.plot(k_num_left_1d, w_num_left_1d, 'rx', label='AIC', markersize=8)
    ax.plot(k_num_left_2d, w_num_left_2d, 'rx', markersize=5)

    ax.plot(k_num_right_1d, w_num_right_1d, 'b+', label='Whistler', markersize=8)
    ax.plot(k_num_right_2d, w_num_right_2d, 'b+', markersize=5)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Wave Number')
    ax.set_ylabel('Angular Frequency')

    ax.legend(loc='upper left', frameon=False)

    fig.tight_layout()
    fig.savefig("dispersion.png", dpi=200)


    error_left_1d  = 100*np.fabs(w_num_left_1d -omega(w_num_left_1d,  -1))/omega(w_num_left_1d,  -1)
    error_right_1d = 100*np.fabs(w_num_right_1d-omega(w_num_right_1d, -1))/omega(w_num_right_1d, +1)

    target_left_1d  = np.array([ 3.,  6.,  1.,  4.,  1.,  2.,  4.,  1.])
    target_right_1d = np.array([ 3.,  6.,  1.,  3.,  0.,  3.,  6., 20.])

    np.testing.assert_allclose(error_left_1d,  target_left_1d,  rtol=1e-2, atol=8)
    np.testing.assert_allclose(error_right_1d, target_right_1d, rtol=1e-2, atol=8)


    error_left_2d  = 100*np.fabs(w_num_left_2d -omega(w_num_left_2d,  -1))/omega(w_num_left_2d,  -1)
    error_right_2d = 100*np.fabs(w_num_right_2d-omega(w_num_right_2d, -1))/omega(w_num_right_2d, +1)

    target_left_2d  = np.array([18.,  9.,  7.,  9.,  6.])
    target_right_2d = np.array([ 4.,  3.,  3.,  9., 29.])

    np.testing.assert_allclose(error_left_2d,  target_left_2d,  rtol=1e-2, atol=8)
    np.testing.assert_allclose(error_right_2d, target_right_2d, rtol=1e-2, atol=8)










class DispersionTest(SimulatorTest):
    # ddt mangles method names so direct lookup isn't easy
    def test_dispersion(self, dim, resolution):
        dispersion(dim, resolution)

@ddt
class DDTDispersionTest(SimulatorTest):
    @data(*permutations)
    @unpack
    def test_dispersion(self, dim, resolution):
        dispersion(dim, resolution)


if __name__=="__main__":

    if len(sys.argv) == 3:
        dim, resolution = parse_cli_args()
        if resolution not in ['low', 'high']:
            raise ValueError('arg should be "low" or "high"')

        test = DispersionTest()
        test.setUp()
        test.test_dispersion(int(dim), resolution)
        test.tearDown()

    elif len(sys.argv) == 1:

        loader = unittest.TestLoader()
        suites = []
        for suite in loader.loadTestsFromTestCase(DDTDispersionTest):
            suites += [suite]
        tests = unittest.TestSuite(suites)
        unittest.TextTestRunner(verbosity=2).run(tests)

    else:
        print('example usage: $script $dim $resolution=[low/high]')

