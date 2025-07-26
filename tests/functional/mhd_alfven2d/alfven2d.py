import numpy as np
import pyphare.mock_mhd_simulator.mhd_model as m
import pyphare.mock_mhd_simulator.simulation as s
from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator


def config():
    alpha = 30.0 * np.pi / 180.0
    cosalpha = np.cos(alpha)
    sinalpha = np.sin(alpha)

    cells = (100, 100)
    dl = ((1.0 / cells[0]) * 1 / cosalpha, (1.0 / cells[1]) * 1 / sinalpha)

    sim = s.Simulation(
        ndim=2,
        order=1,
        timestep=0.002,
        final_time=1.0,
        cells=cells,
        dl=dl,
        origin=(0.0, 0.0),
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        terms="ideal",
        reconstruction="linear",
        limiter="vanleer",
        riemann="rusanov",
        integrator="tvdrk2",
    )

    def density(x, y):
        return 1.0

    def vx(x, y):
        return -0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * sinalpha

    def vy(x, y):
        return 0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * cosalpha

    def vz(x, y):
        return 0.1 * np.cos(2 * np.pi * (x * cosalpha + y * sinalpha))

    def bx(x, y):
        return (
            cosalpha
            - 0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * sinalpha
        )

    def by(x, y):
        return (
            sinalpha
            + 0.1 * np.sin(2 * np.pi * (x * cosalpha + y * sinalpha)) * cosalpha
        )

    def bz(x, y):
        return 0.1 * np.cos(2 * np.pi * (x * cosalpha + y * sinalpha))

    def p(x, y):
        return 0.1

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim


def main():
    MHDMockSimulator(config()).run("alfven2d.h5", dumpfrequency=50)


if __name__ == "__main__":
    main()
