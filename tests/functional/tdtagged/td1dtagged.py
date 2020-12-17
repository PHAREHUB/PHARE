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





def config():
    """ Configure the simulation

    This function defines the Simulation object,
    user initialization model and diagnostics.
    """
    Simulation(
        smallest_patch_size=10,
        largest_patch_size=20,
        time_step_nbr=100,        # number of time steps (not specified if time_step and final_time provided)
        final_time=25,             # simulation final time (not specified if time_step and time_step_nbr provided)
        boundary_types="periodic", # boundary condition, string or tuple, length == len(cell) == len(dl)
        cells=500,                # integer or tuple length == dimension
        dl=1,                  # mesh size of the root level, float or tuple
        max_nbr_levels=3,          # (default=1) max nbr of levels in the AMR hierarchy
        nesting_buffer=2,
        refinement = "tagging",
        diag_options={"format": "phareh5", "options": {"dir": "phare_outputs","mode":"overwrite"}}
    )


    def density(x):
        return 1.


    def S(x,x0,l):
        return 0.5*(1+np.tanh((x-x0)/l))


    def bx(x):
        return 0.


    def by(x):
        from pyphare.pharein.global_vars import sim
        L = sim.simulation_domain()[0]
        v1=-1
        v2=1.
        return v1 + (v2-v1)*(S(x,L*0.25,1) -S(x, L*0.75, 1))


    def bz(x):
        return 0.5


    def b2(x):
        return bx(x)**2 + by(x)**2 + bz(x)**2


    def T(x):
        K = 1
        return 1/density(x)*(K - b2(x)*0.5)


    def vx(x):
        return 0.


    def vy(x):
        return 0.


    def vz(x):
        return 0.


    def vthx(x):
        return T(x)


    def vthy(x):
        return T(x)


    def vthz(x):
        return T(x)


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
    config()
    simulator = Simulator(gv.sim)
    simulator.initialize()
    simulator.run()



if __name__=="__main__":
    main()
