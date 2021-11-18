

import numpy as np
from pyphare.cpp import cpp_lib # must be first
cpp_lib("pybindlibs.cpp_sim_2_1_4")
import pyphare.pharein as ph

seed = 133333333337
cells, dl = 100, .2
patch_sizes = [50,100]
diag_outputs="tools/bench/real/harris/outputs"

def density(x, y):
    L = ph.global_vars.sim.simulation_domain()[1]
    return 0.2 + 1./np.cosh((y-L*0.3)/0.5)**2 + 1./np.cosh((y-L*0.7)/0.5)**2

def by(x, y):
    sim = ph.global_vars.sim
    Lx = sim.simulation_domain()[0]
    Ly = sim.simulation_domain()[1]
    w1, w2 = 0.2, 1.0
    x0 = (x - 0.5 * Lx)
    y1 = (y - 0.3 * Ly)
    y2 = (y - 0.7 * Ly)
    w3 = np.exp(-(x0*x0 + y1*y1) / (w2*w2))
    w4 = np.exp(-(x0*x0 + y2*y2) / (w2*w2))
    w5 = 2.0*w1/w2
    return (w5 * x0 * w3) + ( -w5 * x0 * w4)

def S(y, y0, l): return 0.5*(1. + np.tanh((y-y0)/l))
def bx(x, y):
    sim = ph.global_vars.sim
    Lx = sim.simulation_domain()[0]
    Ly = sim.simulation_domain()[1]
    w1, w2 = 0.2, 1.0
    x0 = (x - 0.5 * Lx)
    y1 = (y - 0.3 * Ly)
    y2 = (y - 0.7 * Ly)
    w3 = np.exp(-(x0*x0 + y1*y1) / (w2*w2))
    w4 = np.exp(-(x0*x0 + y2*y2) / (w2*w2))
    w5 = 2.0*w1/w2
    v1, v2 = -1, 1.
    return v1 + (v2-v1)*(S(y,Ly*0.3,0.5) -S(y, Ly*0.7, 0.5)) + (-w5*y1*w3) + (+w5*y2*w4)

def bz(x, y): return 0.
def b2(x, y): return bx(x,y)**2 + by(x, y)**2 + bz(x, y)**2
def T(x, y):  return 1./density(x, y)*(1 - b2(x, y)*0.5)
def vxyz(x, y): return 0.
def vthxyz(x, y): return np.sqrt(T(x, y))

def config():
    ph.Simulation(# strict=True,
        smallest_patch_size=patch_sizes[0], largest_patch_size=patch_sizes[1],
        time_step_nbr=10, time_step=0.001,
        cells=[cells] * 2, dl=[dl] * 2,
        resistivity=0.001, hyper_resistivity=0.001,
        diag_options={"format": "phareh5", "options": {"dir": diag_outputs, "mode":"overwrite"}},
        refinement_boxes={},
    )
    ph.MaxwellianFluidModel( bx=bx, by=by, bz=bz,
        protons={"charge": 1, "density": density, "init":{"seed": seed},
          **{ "nbr_part_per_cell":100,
            "vbulkx": vxyz, "vbulky": vxyz, "vbulkz": vxyz,
            "vthx": vthxyz, "vthy": vthxyz, "vthz": vthxyz,
          }
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    from tests.diagnostic import all_timestamps
    timestamps = all_timestamps(ph.global_vars.sim)
    timestamps = np.asarray([timestamps[0], timestamps[-1]])
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(
            quantity=quantity,
            write_timestamps=timestamps,
            compute_timestamps=timestamps,
        )

if ph.PHARE_EXE or __name__=="__main__":
    config()

if __name__=="__main__":
    from pyphare.simulator.simulator import Simulator
    Simulator(ph.global_vars.sim).run()
