from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator
import pyphare.mock_mhd_simulator.simulation as s
import pyphare.mock_mhd_simulator.mhd_model as m
import numpy as np

def config():
    sim = s.Simulation(
        ndim = 2,
        order = 1,

        timestep = 0.014,
        final_time = 0.5,

        cells = (128,128),
        dl = (1./128.,1./128.),
        origin = (0.0,0.0),

        eta = 0.0,
        nu = 0.0,
        gamma = 5.0 / 3.0,

        terms = "ideal"
    )

    B0 = 1./(np.sqrt(4.*np.pi))

    def density(x, y):
        return 25./(36.*np.pi)

    def vx(x, y):
        return np.sin(2.*np.pi*y)  # Changed sign

    def vy(x, y):
        return -np.sin(2.*np.pi*x)  # Changed sign

    def vz(x, y):
        return 0.0

    def bx(x, y):
        return B0*np.sin(2.*np.pi*y)  # Changed sign

    def by(x, y):
        return -B0*np.sin(4.*np.pi*x)  # Changed sign

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 5./(12.*np.pi)

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim

def main():
    MHDMockSimulator(config()).run("reversed_orszag_tang.h5", dumpfrequency = 1)

if __name__ == "__main__":
    main()
