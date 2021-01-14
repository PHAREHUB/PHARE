#!/usr/bin/env python3

import pyphare.pharein as ph #lgtm [py/import-and-import-from]
from pyphare.pharein import Simulation
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectromagDiagnostics, FluidDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.simulator.simulator import Simulator
from pyphare.pharein import global_vars as gv


import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
mpl.use('Agg')



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
                      "options": {"dir": "dispersion",
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

    ElectronModel(closure="isothermal", Te=0.12)


    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time +sim.time_step, sim.time_step)
    print(any(timestamps[1:]-timestamps[:-1]< 0.001))
    print(sim.time_step)
    print((timestamps[1:]-timestamps[:-1]).min())
    exit()

    for quantity in ["E", "B"]:
        ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )




def prescribedModes(): #TODOit

    # in this configuration there user prescribed
    # wavelength, at which more energy is thus expected
    # than in modes naturally present in the simulation noise

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
        diag_options={"format": "phareh5", "options": {"dir": "dispersion",
                                                       "mode":"overwrite"}}
    )

    def density(x):
        return 1.


    def by(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.1*np.cos(2*np.pi*x/L[0])


    def bz(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.1*np.sin(2*np.pi*x/L[0])


    def bx(x):
        return 1.


    def vx(x):
        return 0.


    def vy(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.1*np.cos(2*np.pi*x/L[0])

    def vz(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()
        return 0.1*np.sin(2*np.pi*x/L[0])


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
        protons={"charge": 1, "density": density, **vvv, "init": {"seed": 1337}}
    )

    ElectronModel(closure="isothermal", Te=0.12)



    sim = ph.global_vars.sim

    timestamps = np.arange(0, sim.final_time +sim.time_step, sim.time_step)



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








def main():
    fromNoise()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()

if __name__=="__main__":
    main()
