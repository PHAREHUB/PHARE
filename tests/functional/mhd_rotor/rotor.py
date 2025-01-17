import numpy as np
import pyphare.mock_mhd_simulator.mhd_model as m
import pyphare.mock_mhd_simulator.simulation as s
from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator


def config():
    cells = (100, 100)
    dl = (1.0 / cells[0], 1.0 / cells[1])

    sim = s.Simulation(
        ndim=2,
        order=1,
        timestep=0.0003,
        final_time=0.15,
        cells=cells,
        dl=dl,
        origin=(0.0, 0.0),
        eta=0.0,
        nu=0.0,
        gamma=1.4,
        terms="ideal",
        reconstruction="WENO3",
        limiter="minmod",
        riemann="rusanov",
        integrator="tvdrk3",
    )

    B0 = 5 / (np.sqrt(4 * np.pi))
    v0 = 2

    r0 = 0.1
    r1 = 0.115

    def r(x, y):
        return np.sqrt((x - 0.5) ** 2 + (y - 0.5) ** 2)

    def f(r):
        return (r1 - r) / (r1 - r0)

    def density(x, y):
        r_ = r(x, y)
        f_ = f(r_)

        rho_values = np.where(r_ <= r0, 10.0, np.where(r_ < r1, 1.0 + 9.0 * f_, 1.0))
        return rho_values

    def vx(x, y):
        r_ = r(x, y)
        f_ = f(r_)

        vx_values = np.where(
            r_ <= r0,
            -v0 * (y - 0.5) / r0,
            np.where(r_ < r1, -f_ * v0 * (y - 0.5) / r_, 0.0),
        )
        return vx_values

    def vy(x, y):
        r_ = r(x, y)
        f_ = f(r_)

        vy_values = np.where(
            r_ <= r0,
            v0 * (x - 0.5) / r0,
            np.where(r_ < r1, f_ * v0 * (x - 0.5) / r_, 0.0),
        )
        return vy_values

    def vz(x, y):
        return 0.0

    def bx(x, y):
        return B0

    def by(x, y):
        return 0.0

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 1.0

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim


def main():
    MHDMockSimulator(config()).run("rotor.h5", dumpfrequency=100)


if __name__ == "__main__":
    main()
