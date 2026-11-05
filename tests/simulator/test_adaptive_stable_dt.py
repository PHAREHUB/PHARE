#!/usr/bin/env python3

import unittest

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator import SimulatorTest

V = (0.3, -0.2, 0.1)
B = (0.5, 0.4, 0.3)
P = 0.2
GAMMA = 5.0 / 3.0
CFL_WAVE = 0.5


def fast_magnetosonic(gamma, rho, b_dir, b2, p):
    sound = np.sqrt(gamma * p / rho)
    alfven_dir = np.sqrt(b_dir * b_dir / rho)
    alfven = np.sqrt(b2 / rho)
    c02 = sound * sound
    ca2 = alfven * alfven
    cadir2 = alfven_dir * alfven_dir
    return np.sqrt(
        (c02 + ca2) * 0.5
        + np.sqrt((c02 + ca2) * (c02 + ca2) - 4.0 * c02 * cadir2) * 0.5
    )


def whistler(inv_dx, rho, b2):
    v_a = np.sqrt(b2) / np.sqrt(rho)
    x = np.pi * inv_dx / (2.0 * np.sqrt(rho))
    return v_a * (np.sqrt(1.0 + x * x) + x)


def expected_dt(hall, dl, rho0):
    b2 = sum(b * b for b in B)
    inv_dt = 0.0
    for d, dx in enumerate(dl):
        cf = fast_magnetosonic(GAMMA, rho0, B[d], b2, P)
        hall_term = whistler(1.0 / dx, rho0, b2) - np.sqrt(b2 / rho0) if hall else 0.0
        inv_dt += (abs(V[d]) + cf + hall_term) / dx
    return CFL_WAVE / inv_dt


class AdaptiveDtCase(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super().tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def uniform_simulation(self, hall, dl, max_mhd_level=1, time_step=None, **kwargs):
        ph.global_vars.sim = None
        return ph.Simulation(
            interp_order=None,
            model_options=["MHDModel"],
            reconstruction="WENOZ",
            limiter="None",
            riemann="Rusanov",
            mhd_timestepper="TVDRK3",
            hall=hall,
            eta=0.0,
            nu=0.0,
            gamma=GAMMA,
            max_mhd_level=max_mhd_level,
            time_step=time_step or {"mode": "adaptive", "cfl_wave": CFL_WAVE},
            final_time=100.0,
            cells=(20, 10),
            dl=dl,
            **kwargs,
        )

    def uniform_model(self, rho0, pressure=None):
        def const(value):
            return lambda x, y: value + 0.0 * x

        ph.MHDModel(
            density=const(rho0),
            vx=const(V[0]),
            vy=const(V[1]),
            vz=const(V[2]),
            bx=const(B[0]),
            by=const(B[1]),
            bz=const(B[2]),
            p=pressure if pressure is not None else const(P),
        )


@ddt
class AdaptiveStableDtTest(AdaptiveDtCase):
    @data(
        (False, (0.3, 0.2), 1.0),
        (True, (10.0, 5.0), 1.0),
        (True, (0.05, 0.1), 1.0),
        (True, (0.05, 0.1), 0.25),
    )
    @unpack
    def test_uniform_state_dt(self, hall, dl, rho0):
        sim = self.uniform_simulation(hall, dl)
        self.uniform_model(rho0)
        self.simulator = Simulator(sim).initialize()
        dt = self.simulator.timeStep()
        expected = expected_dt(hall, dl, rho0)
        self.assertAlmostEqual(dt / expected, 1.0, delta=1e-12)


if __name__ == "__main__":
    startMPI()
    unittest.main()
