import numpy as np
import pyphare.mock_mhd_simulator.mhd_model as m
import pyphare.mock_mhd_simulator.simulation as s
from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator


def config():
    sim = s.Simulation(
        ndim=2,
        order=1,
        timestep=0.001,
        final_time=0.5,
        cells=(128, 128),
        dl=(1.0 / 128.0, 1.0 / 128.0),
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

    B0 = 1.0 / (np.sqrt(4.0 * np.pi))

    def density(x, y):
        return 25.0 / (36.0 * np.pi)

    def vx(x, y):
        return -np.sin(2.0 * np.pi * y)

    def vy(x, y):
        return np.sin(2.0 * np.pi * x)

    def vz(x, y):
        return 0.0

    def bx(x, y):
        return -B0 * np.sin(2.0 * np.pi * y)

    def by(x, y):
        return B0 * np.sin(4.0 * np.pi * x)

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 5.0 / (12.0 * np.pi)

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim


def main():
    MHDMockSimulator(config()).run("orszag_tang.h5", dumpfrequency=80)


if __name__ == "__main__":
    main()
