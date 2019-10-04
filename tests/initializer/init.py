#! /usr/bin/env python3


import src.initializer.pyphare as pp


def f(x):
    return 2.*x

pp.add("test", f)

#pp.add("simulation/name", "simulation_test")
#pp.add("simulation/ions/nbr_populations", 2)
#pp.add("simulation/ions/pop0/name", "protons")
#pp.add("simulation/ions/pop0/mass", 1.)
#pp.add("simulation/ions/pop0/charge", 1.)
#pp.add("simulation/ions/pop0/particle_initializer/name", "maxwellian")
#pp.add("simulation/ions/pop0/particle_initializer/density", f)



