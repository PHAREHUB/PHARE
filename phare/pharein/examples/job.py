#!/usr/bin/env python


from pharein import Simulation
from pharein import MaxwellianFluidModel
from pharein import ElectromagDiagnostics
from pharein import FluidDiagnostics
from pharein import ParticleDiagnostics



#------------------------------------
#     configure the simulation
#------------------------------------

Simulation(
    time_step_nbr=1000,                   # number of time steps (not specified if time_step and final_time provided)
    final_time=1.,                        # simulation final time (not specified if time_step and time_step_nbr given)
    boundary_types="periodic",            # boundary condition, string or tuple, length == len(cell) == len(dl)
    cells=80,                             # integer or tuple length == dimension
    dl=0.1,                               # mesh size of the root level, float or tuple
    path='test5'                          # directory where INI file and diagnostics directories will be
    # time_step = 0.005,                  # simulation time step (not specified if time_step_nbr and final_time given)
    # domain_size = 8.,                   # float or tuple, not specified if dl and cells are
    # interp_order = 1,                   # interpolation order, [default = 1] can be 1, 2, 3 or 4
    # layout = "yee",                     # grid layout, [default="yee"]
    # origin = 0.,                        # position of the origin of the domain, float or tuple (length = dimension)
    # particle_pusher = "modified_boris", # particle pusher method, [default = "modified_boris"]
    # refined_particle_nbr = 2,           # number of refined particle a particle is split into [default : ]
    # diag_export_format = 'ascii',       # export format of the diagnostics [default = 'ascii']
    # refinement = {"level":[0,1],        # AMR parameters
    #                "extent_ratio":[0.4, 0.6],
    #                "refinement_iterations":[0, 3]},

) # end Simulation







# in the following we usethe MaxwellianFluidModel

import numpy as np

def n(x):
    x0 = 5.
    return 1./np.cosh(x-x0)**2

def bx(x):
    x0=5.
    return np.tanh(x-x0)


MaxwellianFluidModel(bx=bx,
                     protons={"density":n},
                     background={})



ElectromagDiagnostics(
    diag_type="E",                  # available : ("E", "B")
    write_every=10,
    compute_every=5,
    start_teration=0,
    last_iteration=990,
    path = 'ElectromagDiagnostics1'   # where output files will be written, [default: name]
)


FluidDiagnostics(
    diag_type="density",            # choose in (rho_s, flux_s)
    write_every=10,                 # write on disk every x iterations
    compute_every=5,                # compute diagnostics every x iterations ( x <= write_every)
    start_iteration=0,              # iteration at which diag is enabled
    last_iteration=990,             # iteration at which diag is turned off
    population_name="protons"       # name of the population for which the diagnostics is made
  #,path = 'FluidDiagnostics1'      # where output files will be written, [default: name]
)


FluidDiagnostics(
    diag_type="bulk_velocity",
    write_every=10,
    compute_every=5,
    start_iteration=0,
    last_iteration=990,
    population_name="background"
)

FluidDiagnostics(
    diag_type="density",
    write_every=10,
    compute_every=5,
    start_iteration=0,
    last_iteration=990,
    population_name="all"
)

FluidDiagnostics(
    diag_type="flux",
    write_every=10,
    compute_every=5,
    start_iteration=0,
    last_iteration=990,
    population_name="background"
)

ElectromagDiagnostics(
    diag_type="B",
    write_every=10,
    compute_every=5,
    start_teration=0,
    last_iteration=990
)



ParticleDiagnostics(
        compute_every=10,
        write_every=10,
        start_iteration=0,
        last_iteration=90,
        diag_type="space_box",          # choose particles within a spatial box
        extent=(2., 4.),                # extent of the box
        population_name="protons"
        )

