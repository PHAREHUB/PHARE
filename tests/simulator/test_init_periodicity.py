
import unittest
import numpy as np


import pyphare.pharein.global_vars as gv
from pyphare.core  import phare_utilities
from pyphare.pharein import simulation
import pyphare.pharein as ph


def density(x, y): return 1.


def S(y, y0, l):
    return 0.5*(1. + np.tanh((y-y0)/(l)))


def by(x, y): return 0.
def bx(x, y): return 0.
def bz(x, y): return 1.


def b2(x, y):
    return bx(x,y)**2 + by(x, y)**2 + bz(x, y)**2


def T(x, y):
    K = 1
    temp = 1./density(x, y)*(K - b2(x, y)*0.5)
    assert np.all(temp >0)
    return temp

def vx(x, y):
    from pyphare.pharein.global_vars import sim
    y0 = gv.sim.simulation_domain()[1]
    x_b = gv.sim.simulation_domain()[0]
    amp = 0.1
    k = 4
    x0 = x_b*0.25
    x1 = x_b*0.75
    v0 = amp*np.sin(((2*np.pi)/y0)*y*k)*np.exp(-(((x-x0)**2)/2))
    v1 = amp*np.sin(((2*np.pi)/y0)*y*k)*np.exp(-(((x-x1)**2)/2))
    return(v0 + v1)

def vy(x, y):
    from pyphare.pharein.global_vars import sim
    x0 = gv.sim.simulation_domain()[0]
    v1 = -0.5
    v2 = 0.5
    return(v1 + (v2-v1)*(S(x,x0*0.25,0.5) - S(x,x0*0.75,0.5)))


def vz(x, y):
    return 0.

def vthxyz(*xyz): return 0.2



def config():
    ph.Simulation(
        time_step_nbr=2,
        time_step=0.005,
        cells=(200,100),
        dl=(0.4, 0.4),
        refinement="tagging",
        max_nbr_levels = 3,
        tag_buffer = 8,
        nesting_buffer=1,
        refinement_boxes={},
        hyper_resistivity=0.005,
        resistivity=0.001,
        strict=True,
    )

    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthxyz, "vthy": vthxyz, "vthz": vthxyz,
        "nbr_part_per_cell":100
    }

    ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density,  **vvv},

    )

    ph.ElectronModel(closure="isothermal", Te=0.0)




class TestSimulation(unittest.TestCase):

    def setUp(self):
        pass

    def test_periodicity(self):

        config()



if __name__ == "__main__":
    unittest.main()
