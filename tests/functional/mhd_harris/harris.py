import numpy as np
import pyphare.mock_mhd_simulator.mhd_model as m
import pyphare.mock_mhd_simulator.simulation as s
from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator


def config():
    cells = (250, 250)
    dl = (0.1, 0.1)

    sim = s.Simulation(
        ndim=2,
        order=1,
        timestep=0.002,
        final_time=10,
        cells=cells,
        dl=dl,
        origin=(0.0, 0.0),
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        terms="hall",
        reconstruction="WENO3",
        limiter="vanleer",
        riemann="hll",
        integrator="tvdrk3",
    )

    Lx = cells[0] * dl[0]
    Ly = cells[1] * dl[1]

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def density(x, y):
        return (
            1.0
            + 1.0 / np.cosh((y - Ly * 0.3) / 0.5) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / 0.5) ** 2
        )

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def bx(x, y):
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        v1 = -1
        v2 = 1.0
        return (
            v1
            + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
            + (-w5 * y1 * w3)
            + (+w5 * y2 * w4)
        )

    def by(x, y):
        w1 = 0.2
        w2 = 1.0
        x0 = x - 0.5 * Lx
        y1 = y - 0.3 * Ly
        y2 = y - 0.7 * Ly
        w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
        w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
        w5 = 2.0 * w1 / w2
        return (w5 * x0 * w3) + (-w5 * x0 * w4)

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 1.0 - (bx(x, y) ** 2 + by(x, y) ** 2) / 2.0

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim


def main():
    MHDMockSimulator(config()).run("hall_harris.h5", dumpfrequency=500)


if __name__ == "__main__":
    main()
