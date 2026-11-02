#!/usr/bin/env python3
"""
Conservation of the global invariants of the hybrid model, at one regime where
the scheme is known to behave well. Fast enough to run on every pull request.

Setup, periodic box, no resistivity, no hyper-resistivity:
- protons (m=1, q=1, n=1) and alphas (m=4, q=2, n=0.1), vth=0.3, Te=0.2
  (Ti/Te = 0.45), plasma beta ~ 0.7. The energy error here is particle noise,
  which grows with the ion thermal energy relative to the wave energy:
  at Ti = Te (vth=0.45) it is 1.3 to 2 times larger, so Ti is kept below Te
- a standing Alfven wave dB=0.2 with k along the box diagonal, so energy goes
  back and forth between the magnetic field and the ions
- a density modulation of 20%, so the electron pressure term works
- several patches

Checked, on level 0 (see tests/simulator/invariants.py for the definitions):
- ion number per species and mean B: exact, to round-off
- density stays positive, everything finite
- the wave really trades energy with the ions (W_B swings)
- total energy K + W_B + U_e: error bounded, relative to the wave energy;
  in 1D the bound is tight enough that dropping U_e from the sum fails it
- total ion momentum: error bounded, relative to the random-walk level of the
  particles. In the continuous model it is conserved: the force on the ions
  integrates to int(J x B - grad Pe) = 0 in a periodic box

The energy and momentum bounds are not exact properties of the scheme. They
are calibrated over 10 seeds with a margin of about 2.5. 1D is the physics
check: with a scratch copy of the code, 1D order 1 fails on a flipped pressure
sign (energy 31x the bound), a halved pressure term (1.7x, and the U_e check),
and a mass- instead of charge-weighted electron velocity (energy 3x, momentum
5x). 2D and 3D, being noisier, check the exact invariants and are a coarse guard
against heating regressions (a 2x regression passes): in 2D only the flipped pressure sign fails the bound.
Particle loading uses the standard library random distributions, which are
not bitwise portable across compilers: the margin covers that too.
How the energy error depends on the regime (Ti/Te, beta, resolution) is
mapped by the nightly test tests/functional/conservation/conserv.py.
"""

import time
import unittest

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest
from tests.simulator.invariants import invariants

ph.NO_GUI()

Te = 0.2
vth = 0.3  # protons, alphas get vth / sqrt(4)
dB = 0.2
dn = 0.2
dl = 0.5
species = {
    "protons": {"charge": 1, "mass": 1.0, "density": 1.0},
    "alphas": {"charge": 2, "mass": 4.0, "density": 0.1},
}
masses = {name: s["mass"] for name, s in species.items()}

# ndim, interp: cells per direction, patch size, particles per cell, steps,
# and the bounds. Steps run to T=4 in 1D (about a quarter of the wave period),
# T=2 in 2D and 3D. Bounds: max |dE| / E_wave and max |dP| / P_noise, about
# 2.5 times the worst of 10 seeds, and the minimum W_B swing / E_wave
# (worst of 10 seeds: 0.89 in 1D, 0.485 in 2D, 0.589 in 3D).
# Energy, measured worst / bound: 1D 0.0125/0.03, 0.0058/0.015, 0.0037/0.01,
# 2D 0.087/0.2, 0.035/0.09, 0.017/0.04, 3D 0.17/0.4, 0.019/0.05.
# Momentum: 1D 0.10/0.25, 0.036/0.1, 0.016/0.05, 2D 0.27/0.7, 0.11/0.3,
# 0.061/0.15, 3D 0.47/1.2, 0.28/0.7.
# Wall time, Release: about 20 s for the 8 cases.
cases = {
    (1, 1): dict(
        cells=32, patch=8, ppc=100, steps=400, energy=0.03, swing=0.75, momentum=0.25
    ),
    (1, 2): dict(
        cells=32, patch=8, ppc=100, steps=400, energy=0.015, swing=0.75, momentum=0.1
    ),
    (1, 3): dict(
        cells=32, patch=8, ppc=100, steps=400, energy=0.01, swing=0.75, momentum=0.05
    ),
    (2, 1): dict(
        cells=24, patch=12, ppc=50, steps=200, energy=0.2, swing=0.3, momentum=0.7
    ),
    (2, 2): dict(
        cells=24, patch=12, ppc=50, steps=200, energy=0.09, swing=0.3, momentum=0.3
    ),
    (2, 3): dict(
        cells=24, patch=12, ppc=50, steps=200, energy=0.04, swing=0.3, momentum=0.15
    ),
    (3, 1): dict(
        cells=12, patch=6, ppc=20, steps=200, energy=0.4, swing=0.3, momentum=1.2
    ),
    (3, 3): dict(
        cells=12, patch=6, ppc=20, steps=200, energy=0.05, swing=0.3, momentum=0.7
    ),
}
dt = 0.01
nbr_measures = 20


def config(ndim, interp, cells, patch, ppc, steps, seed):
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=dt,
        time_step_nbr=steps,
        cells=[cells] * ndim,
        dl=[dl] * ndim,
        interp_order=interp,
        largest_patch_size=patch,
        smallest_patch_size=patch,
        resistivity=0.0,
        hyper_resistivity=0.0,
        strict=True,
    )

    k = 2 * np.pi / (cells * dl)

    def phase(xyz):
        return k * sum(xyz[:ndim])

    # dB perpendicular to k (box diagonal) and to B0 (x), so div B = 0
    amplitude = {1: (dB, 0.0), 2: (0.0, dB), 3: (dB / np.sqrt(2), -dB / np.sqrt(2))}
    by_amp, bz_amp = amplitude[ndim]

    def constant(value):
        return lambda *xyz: value + 0.0 * xyz[0]

    def population(s, pop_seed):
        v = vth / np.sqrt(s["mass"])
        return {
            "charge": s["charge"],
            "mass": s["mass"],
            "density": lambda *xyz: s["density"] * (1.0 + dn * np.sin(phase(xyz))),
            "vbulkx": constant(0.0),
            "vbulky": constant(0.0),
            "vbulkz": constant(0.0),
            "vthx": constant(v),
            "vthy": constant(v),
            "vthz": constant(v),
            "nbr_part_per_cell": ppc,
            "init": {"seed": pop_seed},
        }

    ph.MaxwellianFluidModel(
        bx=constant(1.0),
        by=lambda *xyz: by_amp * np.cos(phase(xyz)),
        bz=lambda *xyz: bz_amp * np.cos(phase(xyz)),
        **{
            name: population(s, 10 * seed + i)
            for i, (name, s) in enumerate(species.items())
        },
    )
    ph.ElectronModel(closure="isothermal", Te=Te)
    return sim


def run(ndim, interp, seed=1):
    """runs one case, returns the history of the invariants"""
    case = cases[(ndim, interp)]
    sim = config(
        ndim, interp, case["cells"], case["patch"], case["ppc"], case["steps"], seed
    )
    simulator = Simulator(sim).initialize()
    try:
        history = [invariants(simulator, masses, Te, ndim)]
        every = max(1, case["steps"] // nbr_measures)
        for step in range(1, case["steps"] + 1):
            simulator.advance()
            if step % every == 0:
                history.append(invariants(simulator, masses, Te, ndim))
    finally:
        simulator.reset()
        ph.global_vars.sim = None
    return history


def errors(history, ndim, cells):
    """relative errors of each invariant over the run"""
    first = history[0]
    nbr_cells = cells**ndim
    # wave energy, measured: W_B minus the energy of the mean field
    E_wave = first["W_B"] - 0.5 * (first["B_sum"] ** 2).sum() / nbr_cells
    W_B = np.array([h["W_B"] for h in history])
    return {
        "N": max(
            abs(h["N"][name] - first["N"][name]) / first["N"][name]
            for h in history
            for name in first["N"]
        ),
        "B_sum": max(np.abs(h["B_sum"] - first["B_sum"]).max() for h in history)
        / nbr_cells,
        "n_min": min(h["n_min"] for h in history),
        "finite": all(
            np.isfinite([h["E"], h["K"], h["W_B"], h["U_e"], *h["P"]]).all()
            for h in history
        ),
        "swing": (W_B.max() - W_B.min()) / E_wave,
        "energy": max(abs(h["E"] - first["E"]) for h in history) / E_wave,
        "energy_without_U_e": max(
            abs(h["K"] + h["W_B"] - first["K"] - first["W_B"]) for h in history
        )
        / E_wave,
        "momentum": max(np.abs(h["P"] - first["P"]).max() for h in history)
        / first["P_noise"],
    }


@ddt
class ConservationTest(SimulatorTest):
    @data(*cases.keys())
    @unpack
    def test_invariants(self, ndim, interp):
        case = cases[(ndim, interp)]
        start = time.time()
        err = errors(run(ndim, interp), ndim, case["cells"])
        print(
            f"{ndim}D o{interp} ({time.time() - start:.1f} s): "
            + ", ".join(f"{k} {v:.3g}" for k, v in err.items())
        )

        self.assertTrue(err["finite"])
        self.assertLess(err["N"], 1e-12)
        self.assertLess(err["B_sum"], 1e-12)
        self.assertGreater(err["n_min"], 0.0)
        self.assertGreater(err["swing"], case["swing"])
        self.assertLess(err["energy"], case["energy"])
        self.assertLess(err["momentum"], case["momentum"])
        if ndim == 1:
            # the electron term must matter here: without U_e the energy
            # error is 7 times larger or more (10 seeds)
            self.assertLess(3 * err["energy"], err["energy_without_U_e"])


if __name__ == "__main__":
    unittest.main()
