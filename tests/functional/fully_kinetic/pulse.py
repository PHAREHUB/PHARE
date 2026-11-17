#!/usr/bin/env python3
"""
Light pulse across a 1D periodic box, fully kinetic, PHARE units: the P4 gate T8, reflection of
a pulse crossing a static 2:1 refined box (--refinement-box, coarse cells, default [176, 208]),
and its one-level baseline (--one-level: the box only places the pulse and the windows).

Launch. The FK initializer starts from E^0 = 0 (fk_level_initializer.hpp), so the pulse is B_z
alone, b exp(-(x - x0)^2 / 2 sigma^2) sin(k0 (x - x0)), sigma = 8 coarse cells, central
wavelength --wavelength fine cells (32 by default, 16 coarse cells): it splits into a
right-going and a left-going half of amplitude b/2 each. The sine
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

T8 (--check, vacuum, orders 1 and 3, lambda = 32 and 16 fine cells, one filter pass and Marder
s = 0.5 as every refined run, inert in vacuum): the same run with and without the box, both
with the same options; R = ||E-_box - E-_one level|| over the reflection window / incident norm
of the one-level run (the floor and the decomposition leak cancel). Both box edges reflect: the
two reflections arrive 2L = 64 cells apart, wider than the pulse, so they add in norm without
interfering. Model: the 1D Yee model with the P4 coupling run on this exact pulse, geometry and
measure (docs/superpowers/specs/p4/mr_light_1d.py, t8_reflection), with PHARE's fine B ghost, the
conservative dual child B_p + sigma (B_{p+1} - B_{p-1}) / 8 (core/data/grid/gridlayoutdefs.hpp,
directionalProlongation): R = 9.69e-4, 4.06e-3, 2.49e-2, 1.66e-1 at lambda = 64, 32, 16, 8 fine
cells; R(16) / R(32) = 6.12, steeper than (k dx_f)^2 for this broad pulse. Measured (Debug, C++ at
fe21e3bab, order 1): 1.17e-3, 4.89e-3, 2.94e-2, 1.87e-1, 1.12 to 1.20 times the model (the 18%
left is not explained), ratio 6.01. Asserted: |R / R_model - 1| <= 0.3 at lambda = 32 and 16, and
R(16) / R(32) in [5, 7.5]. In vacuum R does not depend on the order (identical to 4 digits at orders
1 to 3, fe21e3bab): order 2 is left to --order 2 by hand. What the ratio rejects, in the model: fine ghosts held at the newest
coarse value instead of interpolated in time (a broken D1) 1.85 (and R x 12 to 39); the first
substep's fine B ghost left from before the reflux 3.94; a one-sided 3/4, 1/4 dual child 3.60.
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

SIGMA = 8  # coarse cells, every length below is in sigma
WAVELENGTHS = (32, 16)  # fine cells, T8
RATIO = 2  # refinement
CELLS = 48 * SIGMA
BOX = (22 * SIGMA, 26 * SIGMA)  # default placeholder refined box
AMPLITUDE = 0.05
STEPS_INCIDENT = 8 * SIGMA  # c dt = dx / 2: the pulse has moved 4 sigma
STEPS_END = 36 * SIGMA
DENSITIES = {"vacuum": 1e-6, "plasma": 1e-2}
PPC = {"vacuum": 4, "plasma": 64}  # per species
TEMPERATURES = {"protons": 0.01, "electrons": 0.04}  # waves.py
FILTER_PASSES, CLEANING = 1, 0.5  # every refined run (P4 D7), the T8 runs
MAX_MODEL_ERROR = 0.3
RATIO_RANGE = (5.0, 7.5)
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


def pulse(geometry, wavelength):
    """B_z of the pulse and the canonical v_y, wavelength in fine cells"""
    k0, sigma = 2 * np.pi / (wavelength / RATIO * dl), SIGMA * dl

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


def run(
    case,
    interp_order,
    geometry,
    wavelength=WAVELENGTHS[0],
    refined=False,
    passes=0,
    cleaning=0.0,
):
    bz, canonical_vy = pulse(geometry, wavelength)
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

    level = "box" if refined else "one"
    diag_dir = f"{diag_root}/{case}_{interp_order}_l{wavelength}_N{passes}_s{cleaning:g}_{level}"
    lo, hi = geometry.box
    refinement = (
        {"refinement_boxes": {"L0": {"B0": [(lo,), (hi - 1,)]}}} if refined else {}
    )
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=STEPS_END,
        cells=[CELLS],
        dl=[dl],
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=passes,
        divergence_cleaning=cleaning,
        **refinement,
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


def reflection(box_dir, one_dir, geometry):
    """||E-_box - E-_one level|| over the reflection window at t_end over the incident norm of
    the one-level run, and the transmitted amplitude error of the box run against it"""
    w = geometry.windows()
    x, plus, _ = characteristics(one_dir, STEPS_INCIDENT)
    incident = norm(x, plus, w["incident"])
    x, plus_one, minus_one = characteristics(one_dir, STEPS_END)
    _, plus_box, minus_box = characteristics(box_dir, STEPS_END)
    return (
        norm(x, minus_box - minus_one, w["reflected"]) / incident,
        norm(x, plus_box, w["transmitted"]) / norm(x, plus_one, w["transmitted"]) - 1,
    )


# the T8 model (docstring), per wavelength in fine cells
MODEL_R = {64: 9.69e-4, 32: 4.06e-3, 16: 2.49e-2, 8: 1.66e-1}


def model_reflection(wavelength):
    return MODEL_R.get(wavelength, float("nan"))


CHECK_ORDERS = (1, 3)


def check(geometry):
    failures = []
    for order in CHECK_ORDERS:
        R = {}
        for wavelength in WAVELENGTHS:
            options = dict(passes=FILTER_PASSES, cleaning=CLEANING)
            one, s1 = run("vacuum", order, geometry, wavelength, **options)
            box, s2 = run(
                "vacuum", order, geometry, wavelength, refined=True, **options
            )
            R[wavelength], transmitted = reflection(box, one, geometry)
            model = model_reflection(wavelength)
            print(
                f"order {order} lambda {wavelength} fine cells: R {R[wavelength]:.3e}, model "
                f"{model:.3e}, R / model - 1 {R[wavelength] / model - 1:+.2f}, transmitted error "
                f"{transmitted:+.1e}, runs {s1:.1f} + {s2:.1f} s"
            )
            if abs(R[wavelength] / model - 1) > MAX_MODEL_ERROR:
                failures.append(f"order {order} lambda {wavelength}: R off the model")
        ratio = R[16] / R[32]
        print(f"order {order}: R(16) / R(32) = {ratio:.2f}, in {RATIO_RANGE}")
        if not RATIO_RANGE[0] <= ratio <= RATIO_RANGE[1]:
            failures.append(f"order {order}: R(16) / R(32) = {ratio:.2f}")
    if failures:
        raise RuntimeError("; ".join(failures))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--order", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument("--case", choices=[*DENSITIES, "both"], default="both")
    parser.add_argument(
        "--refinement-box",
        default=f"{BOX[0]},{BOX[1]}",
        help="lo,hi in coarse cells (hi excluded): the refined box, and the pulse and windows",
    )
    parser.add_argument(
        "--wavelength", type=int, default=WAVELENGTHS[0], help="fine cells"
    )
    parser.add_argument(
        "--one-level", action="store_true", help="no refinement: the baseline"
    )
    parser.add_argument("--filter-passes", type=int, default=0)
    parser.add_argument("--cleaning", type=float, default=0.0)
    parser.add_argument("--check", action="store_true", help="T8, the registered test")
    args = parser.parse_args()
    geometry = Geometry(tuple(int(c) for c in args.refinement_box.split(",")))
    if args.check:
        return check(geometry)

    options = dict(passes=args.filter_passes, cleaning=args.cleaning)
    for case in DENSITIES if args.case == "both" else [args.case]:
        one, seconds = run(case, args.order, geometry, args.wavelength, **options)
        results = measure(one, geometry)
        print(
            f"{case} order {args.order} one level: "
            + ", ".join(f"{name} {value:+.3e}" for name, value in results.items())
            + f", run {seconds:.1f} s"
        )
        if args.one_level:
            continue
        box, seconds = run(
            case, args.order, geometry, args.wavelength, refined=True, **options
        )
        R, transmitted = reflection(box, one, geometry)
        print(
            f"{case} order {args.order} refined box: R {R:.3e} (model "
            f"{model_reflection(args.wavelength):.3e}), transmitted error {transmitted:+.1e},"
            f" run {seconds:.1f} s"
        )


if __name__ == "__main__":
    main()
