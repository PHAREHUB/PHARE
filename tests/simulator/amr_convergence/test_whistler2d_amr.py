#!/usr/bin/env python3
"""AMR space-time convergence: 2D rotated Hall-MHD whistler wave (Toth,
JCP 227, 2008, Table 2). Hall-MHD analog of test_alfven2d_amr; the ideal
CP-Alfven wave is replaced by a weakly-dispersive whistler (k*d_i = 0.19),
which is the characteristic Hall wave and therefore the sensitive probe for
the coarse-fine time interpolation of the Hall term (B/E/J refiners). See
amr_convergence_base for the protocol and compute_errors for the norm.

Reconstruction is fixed at WENOZ -- the whistler is dispersive, so Linear
would not reach order 2 and could not expose an order defect. Hall dt scaling
(unlike Alfven's dt ~ dx): the whistler grid dispersion gives the stability
bound dt <~ dx^2/pi, so at fixed sigma dt ~ dx^2 -- large N is expensive
(N=128 needs ~2800 steps at sigma=0.4).

Requires the build permutation
  2,1,4,SSPRK4_5,WENOZ,None,Rusanov,true,false,false  (in res/sim/all.txt).
"""

import os
import unittest

import numpy as np
from ddt import data, ddt

import pyphare.pharein as ph
from tests.simulator.amr_convergence.amr_convergence_base import ConvergenceTestBase

os.environ.setdefault("PHARE_SCOPE_TIMING", "0")

ph.NO_GUI()

# Rotated whistler, Toth Table 2: alpha = arctan(0.5) ~= 26.56 deg.
tan_alpha = 0.5
cosalpha = 1.0 / np.sqrt(1 + tan_alpha**2)
sinalpha = tan_alpha / np.sqrt(1 + tan_alpha**2)

B0, RHO0, P0 = 1.0, 1.0, 1.0
GAMMA = 5.0 / 3.0
DELTA = 1e-3  # wave amplitude (Toth dB/B0 = 1e-3)

c_A = B0 / np.sqrt(RHO0)
c_s = np.sqrt(GAMMA * P0 / RHO0)
v_fast = np.sqrt(c_s**2 + c_A**2)

# Toth k*d_i ~= 0.19 (PHARE d_i = 1 -> k = 0.19); one mode over the domain.
K = 0.19
L = 2 * np.pi / K  # wavelength = domain length along e1

OMEGA = 0.5 * K**2 + np.sqrt((K * c_A) ** 2 + (0.5 * K**2) ** 2)
C_W = OMEGA / K  # whistler phase speed
V_AMP = DELTA * c_A / C_W  # transverse velocity amplitude (Toth eq. 55)

# Wave frame (x-y plane tilt): e1 = k || B0, e2/e3 transverse.
e1 = np.array([cosalpha, sinalpha, 0.0])
e2 = np.array([-sinalpha, cosalpha, 0.0])
e3 = np.cross(e1, e2)  # = (0, 0, 1)
B0_vec = B0 * e1

# Lengthened rotated domain: one wavelength along e1 fits each periodic axis.
Lx = L / cosalpha
Ly = L / sinalpha

TIMESTEPPER = "SSPRK4_5"
RECONSTRUCTION = "WENOZ"
LIMITER = "None"


@ddt
class WhistlerConvergenceTest(ConvergenceTestBase):
    name = "whistler2d"
    final_time = L / C_W  # one traversal = 2*pi/omega = exact return to IC

    # N=16 preasymptotic here too (same family as the alfven test); dt ~ dx^2
    # makes N=128 the expensive end (~2800 steps at sigma=0.4).
    SPATIAL_NS = [32, 64, 128]
    SPATIAL_SIGMA = 0.4  # proven operating point (uniform converges @ 2.00)
    # Measured on kaa 2026-07-09: order=2 slope 1.99 (segments 1.99 -> 2.00).
    SPATIAL_ORDER_BAND = (1.70, 2.25)

    # sigma sweep shared with the alfven test; gate = drift of the AMR/uniform
    # error ratio (see base class). The J-refiner fixed-N sigma-sweep study
    # (2026-07-07/08) validated this span for the whistler: sigma=0.9 is
    # stable up to N=128 (N=256 blows up there; fewer steps at N=64 leaves
    # more margin). That study measured the static C-F refiner defect as ~6%
    # L1(fine) drift over a comparable span, and the fixed (co-temporal
    # time-refiner) pipeline as flat -- the gate sits between. Calibration on
    # kaa 2026-07-09: raw errors flat to 0.1% (dt ~ dx^2 keeps SSPRK4_5's
    # temporal error far below the spatial floor), ratio drifts 0.0005 /
    # 0.0006 -- any C-F defect stands out against a near-zero baseline.
    SWEEP_N = 64
    SWEEP_SIGMAS = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    MAX_AMR_SIGMA_DRIFT = 0.05

    def cfl_dt(self, N):
        """Per-unit-sigma stable step. Hall dispersive bound dt ~ dx^2/pi
        dominates at fine resolution; the fast-mode bound dx/v_fast caps it at
        coarse N. Uses the smaller grid spacing (worst case)."""
        dxmin = min(Lx / N, Ly / N)
        return min(dxmin**2 / np.pi, dxmin / v_fast)

    def _common(self, N, n):
        return dict(
            smallest_patch_size=8,
            time_step=self.final_time / n,
            final_time=self.final_time,
            cells=(N, N),
            dl=(Lx / N, Ly / N),
            interp_order=1,
            hyper_resistivity=0.0,
            resistivity=0.0,
            strict=True,
            nesting_buffer=1,
            eta=0.0,
            nu=0.0,
            gamma=GAMMA,
            reconstruction=RECONSTRUCTION,
            limiter=LIMITER,
            riemann="Rusanov",
            mhd_timestepper=TIMESTEPPER,
            hall=True,
            res=False,
            hyper_res=False,
            model_options=["MHDModel"],
        )

    def amr_simulation(self, order, N, n):
        extra = {}
        tag = f"o{order}"
        if order != 0:
            extra["refinement_order"] = order
        base = f"phare_outputs/{self.name}_amr_convergence/{tag}_N{N}_n{n}"
        return self.simulation(
            refinement="boxes",
            refinement_boxes={"L0": {"B0": self.fine_box(N)}},
            max_mhd_level=2,
            diag_options={
                "format": "phareh5",
                "options": {"dir": base, "mode": "overwrite"},
            },
            **self._common(N, n),
            **extra,
        )

    def uniform_simulation(self, N, n):
        base = f"phare_outputs/{self.name}_amr_convergence/uniform_N{N}_n{n}"
        return self.simulation(
            refinement="tagging",
            max_mhd_level=1,
            max_nbr_levels=1,
            diag_options={
                "format": "phareh5",
                "options": {"dir": base, "mode": "overwrite"},
            },
            **self._common(N, n),
        )

    def add_model_and_diags(self):
        def phase(x, y):
            return K * (x * e1[0] + y * e1[1])  # t=0

        def dB(x, y):
            p = phase(x, y)
            # right-hand circularly polarized (Toth)
            return DELTA * (np.cos(p)[:, None] * e2 - np.sin(p)[:, None] * e3)

        def dV(x, y):
            p = phase(x, y)
            return V_AMP * (-np.cos(p)[:, None] * e2 + np.sin(p)[:, None] * e3)

        def density(x, y):
            return RHO0

        def p(x, y):
            return P0

        def vx(x, y):
            return dV(x, y)[:, 0]

        def vy(x, y):
            return dV(x, y)[:, 1]

        def vz(x, y):
            return dV(x, y)[:, 2]

        def bx(x, y):
            return B0_vec[0] + dB(x, y)[:, 0]

        def by(x, y):
            return B0_vec[1] + dB(x, y)[:, 1]

        def bz(x, y):
            return B0_vec[2] + dB(x, y)[:, 2]

        ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

        # dump the conserved set natively so the comparison variables are the
        # dumped cell/face averages directly.
        timestamps = [0.0, self.final_time]
        ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)
        for quantity in ["rho", "rhoV", "Etot"]:
            ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    @data(2)  # Linear2 runtime kernel
    def test_spatial_convergence(self, refinement_order):
        self.check_spatial_order(
            refinement_order, self.SPATIAL_NS, self.SPATIAL_SIGMA,
            self.SPATIAL_ORDER_BAND,
        )

    # legacy default refinement = order-1 refinement (order-0 polynomial;
    # linear only on non-shared edges/faces). The whistler pins its cost even
    # more sharply than the alfven case (1.58 there): measured on kaa
    # 2026-07-09 slope 1.28 (segments 1.36 -> 1.19 still declining), L1 ~ 7x
    # L0 at N=128 -- the dispersive wave amplifies the C-F pollution.
    # Documented defect, expected to fail the 2nd-order band until the legacy
    # op is replaced (an unexpected success flags that it was).
    @unittest.expectedFailure
    def test_spatial_convergence_legacy(self):
        self.check_spatial_order(
            0, self.SPATIAL_NS, self.SPATIAL_SIGMA, self.SPATIAL_ORDER_BAND
        )

    def test_temporal_sigma_sweep(self):
        self.check_sigma_sweep(order=2, N=self.SWEEP_N, sigmas=self.SWEEP_SIGMAS)


if __name__ == "__main__":
    unittest.main()
