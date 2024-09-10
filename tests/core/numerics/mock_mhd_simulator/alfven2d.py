from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator
import pyphare.mock_mhd_simulator.simulation as s
import pyphare.mock_mhd_simulator.mhd_model as m
import numpy as np

def config():
    cells = (200,100)
    dl = (0.02,0.02)
    sim = s.Simulation(
        ndim = 2,
        order = 1,

        timestep = 0.003,
        final_time = 0.7,

        cells = cells,
        dl = dl,
        origin = (0.0,0.0),

        eta = 0.0,
        nu = 0.0,
        gamma = 5.0 / 3.0,

        terms = "ideal"
    )

    lx = cells[0] * dl[0]
    ly = cells[1] * dl[1]

    alpha = np.arctan(0.5)
    cosalpha = np.cos(alpha)
    sinalpha = np.sin(alpha)
    kx = 2 * np.pi / (cosalpha * lx)
    ky = 2 * np.pi / (sinalpha * ly)

    def density(x, y):
        return 1.0

    def vx(x, y):
        return (-1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2.0 * sinalpha)

    def vy(x, y):
        return (1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2.0 * cosalpha)

    def vz(x, y):
        return 0.0

    def bx(x, y):
        return (1 - 1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2.0 * sinalpha)

    def by(x, y):
        return (1 + 1e-4 * np.sin(kx * x * cosalpha + ky * y * sinalpha)) / (2.0 * cosalpha)

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 0.1

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim

def main():
    MHDMockSimulator(config()).run("alfven2d.h5")

if __name__ == "__main__":
    main()
