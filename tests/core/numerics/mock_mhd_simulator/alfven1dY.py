from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator
import pyphare.mock_mhd_simulator.simulation as s
import pyphare.mock_mhd_simulator.mhd_model as m
import numpy as np

def config():
    sim = s.Simulation(
        ndim = 1,
        order = 1,

        timestep = 0.01,
        final_time = 1.0,

        cells = (50,),
        dl = (0.02,),
        origin = (0.0,),

        eta = 0.0,
        nu = 0.0,
        gamma = 5.0 / 3.0,

        terms = "ideal"
    )

    kx = 2.0 * np.pi

    def density(x):
        return 1.0

    def vx(x):
        return -1e-6 * np.cos(kx * x)

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def bx(x):
        return 1e-6 * np.cos(kx * x)

    def by(x):
        return 1.0

    def bz(x):
        return 0.0

    def p(x):
        return 0.1

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim

def main():
    MHDMockSimulator(config()).run("alfven1d.h5")

if __name__ == "__main__":
    main()
