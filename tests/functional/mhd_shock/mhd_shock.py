import numpy as np
import pyphare.mock_mhd_simulator.mhd_model as m
import pyphare.mock_mhd_simulator.simulation as s
from pyphare.mock_mhd_simulator.simulator import MHDMockSimulator


def config():
    cells = (800,)
    dl = (1.0,)

    sim = s.Simulation(
        ndim=1,
        order=1,
        timestep=0.2,
        final_time=80,
        cells=cells,
        dl=dl,
        origin=(0.0,),
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        terms="ideal",
        reconstruction="constant",
        limiter="minmod",
        riemann="hll",
        integrator="tvdrk2",
    )

    def density(x):
        return np.where(x < (cells[0] * dl[0] / 2), 1, 0.125)

    def vx(x):
        return 0.0

    def vy(x):
        return 0.0

    def vz(x):
        return 0.0

    def bx(x):
        return 0.75

    def by(x):
        return np.where(x < (cells[0] * dl[0] / 2), 1, -1)

    def bz(x):
        return 0.0

    def p(x):
        return np.where(x < (cells[0] * dl[0] / 2), 1, 0.1)

    m.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    return sim


def main():
    MHDMockSimulator(config()).run("shock.h5", dumpfrequency=10)


if __name__ == "__main__":
    main()
