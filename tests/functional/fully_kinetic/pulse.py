#!/usr/bin/env python3
"""
Light pulse across a 1D periodic box, fully kinetic, PHARE units: the one-level baseline of the
P4 gate T8 (reflection and transmission of a pulse crossing a static refined box). The refined
box is a placeholder (--refinement-box, cells): it only places the pulse and the windows.

Launch. The FK initializer starts from E^0 = 0 (fk_level_initializer.hpp), so the pulse is B_z
alone, b exp(-(x - x0)^2 / 2 sigma^2) sin(k0 (x - x0)), sigma = 8 cells, central wavelength 32
cells: it splits into a right-going and a left-going half of amplitude b/2 each. The sine
carrier makes the integral of B_z zero, so the particles start with the canonical
v_y = -(q/m) A_y, A_y' = B_z, and keep no memory of the initial field (waves.py, light).

Decomposition. With B_z averaged onto the E_y nodes (dual -> primal) and to E's time (the E
dumped at t is E^{n-1/2}: B dumped at t - dt and t), E+- = (E_y +- c B_z) / 2 is the right-going
and the left-going part of a vacuum light wave. The averages leave (1 - cos(k dx/2) cos(w dt/2))
/ 2 = 3e-3 of a pulse in the other part (measured 3.6e-3), where that pulse is: never in the
windows below.

Measurement (the one T8 uses with a refined box, on level 0):
- incident: ||E+|| over x_R(t_inc) +- 4 sigma, at t_inc = 8 sigma steps, x_R(t) = x0 + c t
- reflected / incident: ||E-|| over the reflection window [box_lo - 16 sigma, box_lo] at
  t_end = 36 sigma steps (both box edges have reflected and those reflections have left the
  box), over the incident norm
- transmitted amplitude error: ||E+|| over x_R(t_end) +- 4 sigma, over the incident norm, - 1
- check of the decomposition: ||E-|| of the left-going half at t_inc over the incident norm (1)
||f|| = sqrt(sum f^2 dx) over the physical nodes of the window. Geometry (default box
[176, 208] of 384 cells): x0 = box_lo - 6 sigma; windows, in cells, incident [128, 192],
reflected [48, 176], transmitted [240, 304]; the left-going half stays 8 sigma away from them,
periodic wrap included. The dumps are single precision: ratios below ~1e-7 are not resolved.

Densities.
- vacuum: n = 1e-6, 4 particles per cell and species: (w_pe / c k0)^2 = 1e-6, and particle
  noise far below the reflections T8 has to resolve (3.4e-3 expected at 64 fine cells per
  wavelength, from the 1D model of the P4 review).
- plasma: the species of waves.py (m_i/m_e = 25, Te = 0.04, Ti = 0.01) at n = 0.01, 64
  particles per cell and species: w_pe = 5, (w_pe / c k0)^2 = 0.010, group velocity 0.995 c,
  7e-4 of the pulse energy below the cutoff. At waves.py's n = 1, w_pe = 50 equals c k0 = 49:
  the pulse would sit at the cutoff and not propagate. As waves.py, v_th,e = 0.1 c is on the
  0.3 c guard of the strict setup checks (0.33 with the canonical v_y): strict=False.

Measured, seed fixed, Debug build, C++ at 1714f3d49 (one level, so the "reflection" is the
floor of the measurement):
- vacuum, orders 1 and 2: reflected / incident 2.5e-5, transmitted error 7e-8 (float32 dumps),
  left half 1.000; 3.4 s, 3.8 s
- plasma, orders 1 and 2: reflected / incident 4.6e-2 (thermal noise, E- flat at 2e-3 to 4e-3
  over the box), transmitted error +1.4e-3, left half 0.998; 37.7 s, 43.2 s
The plasma floor is above every reflection T8 expects: T8 gates on the vacuum case.
"""

import argparse
import time
from dataclasses import dataclass

import numpy as np
import waves
from waves import c_over_va, dl, mass_e, time_step

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

SIGMA = 8  # cells, every length below is in sigma
WAVELENGTH = 32  # cells
CELLS = 48 * SIGMA
BOX = (22 * SIGMA, 26 * SIGMA)  # default placeholder refined box
AMPLITUDE = 0.05
STEPS_INCIDENT = 8 * SIGMA  # c dt = dx / 2: the pulse has moved 4 sigma
STEPS_END = 36 * SIGMA
DENSITIES = {"vacuum": 1e-6, "plasma": 1e-2}
PPC = {"vacuum": 4, "plasma": 64}  # per species
TEMPERATURES = {"protons": 0.01, "electrons": 0.04}  # waves.py
diag_root = "phare_outputs/fully_kinetic_pulse"


@dataclass(frozen=True)
class Geometry:
    box: tuple  # (lo, hi), cells

    @property
    def x0(self):
        return (self.box[0] - 6 * SIGMA) * dl

    def x_right(self, steps):
        return self.x0 + c_over_va * steps * time_step

    def windows(self):
        """physical intervals of the incident, reflected and transmitted windows, and of the
        left-going half at t_inc"""
        sigma = SIGMA * dl
        lo = self.box[0] * dl
        x_inc, x_end = self.x_right(STEPS_INCIDENT), self.x_right(STEPS_END)
        x_left = self.x0 - c_over_va * STEPS_INCIDENT * time_step
        return {
            "incident": (x_inc - 4 * sigma, x_inc + 4 * sigma),
            "reflected": (lo - 16 * sigma, lo),
            "transmitted": (x_end - 4 * sigma, x_end + 4 * sigma),
            "left half": (x_left - 4 * sigma, x_left + 4 * sigma),
        }


def pulse(geometry):
    k0, sigma = 2 * np.pi / (WAVELENGTH * dl), SIGMA * dl

    def bz(x):
        u = x - geometry.x0
        return AMPLITUDE * np.exp(-(u**2) / (2 * sigma**2)) * np.sin(k0 * u)

    # A_y = int_0^x B_z, tabulated (the integral over the box is zero)
    xs = np.linspace(0, CELLS * dl, 64 * CELLS + 1)
    b = bz(xs)
    ay = np.concatenate([[0], np.cumsum((b[1:] + b[:-1]) / 2 * np.diff(xs))])

    def canonical_vy(q, m):
        def vy(x):
            return -q / m * np.interp(x, xs, ay)

        return vy

    return bz, canonical_vy


def run(case, interp_order, geometry):
    bz, canonical_vy = pulse(geometry)
    density, ppc = DENSITIES[case], PPC[case]
    populations = {
        "protons": waves.maxwellian(
            1, 1.0, density, TEMPERATURES["protons"], ppc=ppc, seed=1
        ),
        "electrons": waves.maxwellian(
            -1, mass_e, density, TEMPERATURES["electrons"], ppc=ppc, seed=2
        ),
    }
    for pop in populations.values():
        pop["vbulky"] = canonical_vy(pop["charge"], pop["mass"])

    diag_dir = f"{diag_root}/{case}_{interp_order}"
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=STEPS_END,
        cells=[CELLS],
        dl=[dl],
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=False,
    )
    ph.MaxwellianFluidModel(
        bx=waves.uniform(0.0), by=waves.uniform(0.0), bz=bz, **populations
    )
    steps = [STEPS_INCIDENT, STEPS_END]
    ph.ElectromagDiagnostics(
        quantity="E", write_timestamps=[s * time_step for s in steps]
    )
    ph.ElectromagDiagnostics(
        quantity="B",
        write_timestamps=sorted(s * time_step for n in steps for s in (n - 1, n)),
    )
    start = time.perf_counter()
    Simulator(sim).run().reset()
    return diag_dir, time.perf_counter() - start


def level0(diag_dir, quantity, time_):
    """positions and values of a component on the physical nodes of level 0, sorted"""
    hier = hierarchy_from(
        h5_filename=f"{diag_dir}/EM_{quantity[0]}.h5", times=f"{time_:.10f}"
    )
    length = CELLS * dl
    xs, vs = [], []
    for patch in hier.level(0, time_).patches:
        pdata = patch.patch_datas[quantity]
        (x,) = pdata.meshgrid()
        inside = (x > -1e-6) & (x < length - 1e-6)
        xs.append(x[inside])
        vs.append(pdata.dataset[:][inside])
    x, v = np.concatenate(xs), np.concatenate(vs).astype(np.float64)
    order = np.argsort(x)
    return x[order], v[order]


def characteristics(diag_dir, steps):
    """x and E+-, at the time and on the nodes of E_y"""
    x, ey = level0(diag_dir, "Ey", steps * time_step)
    _, bz_before = level0(diag_dir, "Bz", (steps - 1) * time_step)
    _, bz_after = level0(diag_dir, "Bz", steps * time_step)
    bz = (bz_before + bz_after) / 2  # B at t - dt/2, on the dual nodes x_{i+1/2}
    bz = (bz + np.roll(bz, 1)) / 2  # onto the primal nodes x_i (periodic)
    return x, (ey + c_over_va * bz) / 2, (ey - c_over_va * bz) / 2


def norm(x, f, window):
    """sqrt(sum f^2 dx) over the nodes inside a window, periodic"""
    length = CELLS * dl
    lo, hi = window
    inside = ((x - lo) % length) <= (hi - lo)
    return np.sqrt(np.sum(f[inside] ** 2) * dl)


def measure(diag_dir, geometry):
    w = geometry.windows()
    x, plus, minus = characteristics(diag_dir, STEPS_INCIDENT)
    incident = norm(x, plus, w["incident"])
    left_half = norm(x, minus, w["left half"]) / incident
    x, plus, minus = characteristics(diag_dir, STEPS_END)
    return {
        "reflected / incident": norm(x, minus, w["reflected"]) / incident,
        "transmitted error": norm(x, plus, w["transmitted"]) / incident - 1,
        "left half / incident": left_half,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--order", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument("--case", choices=[*DENSITIES, "both"], default="both")
    parser.add_argument(
        "--refinement-box",
        default=f"{BOX[0]},{BOX[1]}",
        help="lo,hi in cells; placeholder: places the pulse and the windows, no refinement yet",
    )
    args = parser.parse_args()
    geometry = Geometry(tuple(int(c) for c in args.refinement_box.split(",")))

    for case in DENSITIES if args.case == "both" else [args.case]:
        diag_dir, seconds = run(case, args.order, geometry)
        results = measure(diag_dir, geometry)
        print(
            f"{case} order {args.order}: "
            + ", ".join(f"{name} {value:+.3e}" for name, value in results.items())
            + f", run {seconds:.1f} s"
        )


if __name__ == "__main__":
    main()
