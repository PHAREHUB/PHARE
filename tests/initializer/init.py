#! /usr/bin/env python3


import phare.pyphare as pp
import job


def density(x):
    return 2.*x

def vx(x):
    return 2*x

def vy(x):
    return 3*x

def vz(x):
    return 4*x


def vthx(x):
    return 5*x

def vthy(x):
    return 6*x

def vthz(x):
    return 7*x


simulation = job.ph.globals.sim

pp.add("simulation/name", "simulation_test")
pp.add("simulation/dimension", simulation.dims)


pp.add("simulation/grid/layout_type", simulation.layout)
pp.add("simulation/grid/nbr_cells/x", simulation.cells[0])
pp.add("simulation/grid/meshsize/x", simulation.dl[0])
pp.add("simulation/grid/origin/x", simulation.origin[0])

if (simulation.dims>1):
    pp.add("simulation/grid/nbr_cells/y", simulation.cells[1])
    pp.add("simulation/grid/meshsize/y", simulation.dl[1])
    pp.add("simulation/grid/origin/y", simulation.origin[1])

    if (simulation.dims >2):
        pp.add("simulation/grid/nbr_cells/z", simulation.cells[2])
        pp.add("simulation/grid/meshsize/z", simulation.dl[2])
        pp.add("simulation/grid/origin/z", simulation.origin[2])



pp.add("simulation/interp_order", simulation.interp_order)
pp.add("simulation/time_step", simulation.time_step)


pp.add("simulation/algo/pusher/name", simulation.particle_pusher)

pp.add("simulation/ions/nbr_populations", 2)
pp.add("simulation/ions/pop0/name", "protons")
pp.add("simulation/ions/pop0/mass", 1.)
pp.add("simulation/ions/pop0/particle_initializer/name", "maxwellian")

pp.addScalarFunction1D("simulation/ions/pop0/particle_initializer/density", density)

pp.addScalarFunction1D("simulation/ions/pop0/particle_initializer/bulk_velocity_x", vx)
pp.addScalarFunction1D("simulation/ions/pop0/particle_initializer/bulk_velocity_y", vy)
pp.addScalarFunction1D("simulation/ions/pop0/particle_initializer/bulk_velocity_z", vz)

pp.addScalarFunction1D("simulation/ions/pop0/particle_initializer/thermal_velocity_x", vthx)
pp.addScalarFunction1D("simulation/ions/pop0/particle_initializer/thermal_velocity_y", vthy)
pp.addScalarFunction1D("simulation/ions/pop0/particle_initializer/thermal_velocity_z", vthz)

pp.add("simulation/ions/pop0/particle_initializer/nbr_part_per_cell", 100)
pp.add("simulation/ions/pop0/particle_initializer/charge", 1.)
pp.add("simulation/ions/pop0/particle_initializer/basis", "cartesian")




