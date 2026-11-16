#!/usr/bin/env python3
"""
FullyKineticModel physics checks in 2D, periodic, PHARE units, with the conventions and grid
corrections of waves.py (1D). The grid corrections are per direction: the particle shape
factor is the product over directions, S(k)^2 = prod_d sinc(k_d dx / 2)^(2(p+1)), and the Yee
derivative along d acts as K_d = (2/dx) sin(k_d dx / 2).

- langmuir: oblique, k = (1, 1) 2 pi / L, electron velocity v0 k/|k| sin(k.x),
  w^2 = (k_x / K_x) S^2 (w_pe^2 (1 + 3 |k|^2 lambda_De^2) + w_pi^2): with k_x = k_y, the
  Esirkepov current of a displacement along d carries k_d / K_d in each direction alike
- light: oblique, both polarizations with the same dispersion,
  w^2 = S^2 w_p^2 + c^2 sum_d K_d^2 (+ warm transverse term as in 1D), canonical start
  v = -(q/m) A:
  "in plane": B_z = b cos(k.x), in-plane E and current, A = (b / |k|^2) (-k_y, k_x) sin(k.x)
  "out of plane": B = b (k_y, -k_x) / |k| cos(k.x), E_z and J_z, A_z = (b / |k|) sin(k.x):
  the current deposited out of the plane and both in-plane B components
- weibel: two electron beams +-v0 along z (out of the plane), k in the plane so k.v0 = 0: no
  two-stream. Warm in the plane (v_th = 0.72: dx / lambda_D = 2.8), cold along z. Kinetic
  dispersion of the transverse mode E_z, for a purely growing w = i gamma:
    -gamma^2 - c^2 K^2 - S^2 w_p^2 + S^2 sum_s w_ps^2 A_s (1 - sqrt(pi) y_s erfcx(y_s)) = 0,
    y_s = gamma / (sqrt(2) |k| v_th,s), A_s = (<v_z^2>_s) / v_th,s^2,
  which reduces to the cold w^4 - (w_p^2 + c^2 k^2) w^2 - w_pe^2 k^2 v0^2 = 0. In a 16^2 box
  only modes (1, 0), (0, 1) (gamma = 5.1) and (1, +-1) (gamma = 1.7) are unstable: no faster
  mode grows from the noise. Measured on (1, 0) and (0, 1) together, from the noise.
  Not a test yet (--weibel): measured in 1D with the same 16k particles per beam, noise
  2e-4 to 1e-3 and saturation 0.09 leave a ~2.8 e-folding window where the fitted gamma moves
  by -6% to +11% with the window bounds.
"""

import sys

import numpy as np
import waves
from waves import (
    c_over_va,
    check,
    dl,
    exponential_growth,
    k_over_K,
    leapfrog,
    mass_e,
    maxwellian,
    mode_amplitude,
    plasma_frequency,
    run,
    standing_fit,
)

# 1D values (waves.py), to be recalibrated over 5 seeds in 2D
TOLERANCES = {"langmuir": 7e-3, "light": 4e-3, "weibel": 5e-2}
MAX_LIGHT_DAMPING = waves.MAX_LIGHT_DAMPING


def shape_factor2(k, interp_order):
    return np.prod([waves.shape_factor2(kd, interp_order) for kd in k])


def yee_k2(k):
    return np.sum((2 / dl * np.sin(np.asarray(k) * dl / 2)) ** 2)


def wavevector(cells, mode):
    return 2 * np.pi * np.asarray(mode) / (np.asarray(cells) * dl)


def langmuir(interp_order, seed=0, v0=0.15, Te=0.04, ppc=16):
    cells, mode = (32, 32), (1, 1)
    k = wavevector(cells, mode)
    unit = k / np.linalg.norm(k)

    def velocity(direction):
        def fn(x, y):
            return v0 * unit[direction] * np.sin(k[0] * x + k[1] * y)

        return fn

    electrons = maxwellian(-1, mass_e, 1.0, Te, ppc=ppc, seed=2 + 10 * seed)
    electrons["vbulkx"], electrons["vbulky"] = velocity(0), velocity(1)
    diag_dir, times = run(
        "langmuir_2d",
        interp_order,
        cells,
        final_time=10 * 2 * np.pi / plasma_frequency(1, 1, mass_e),
        populations={
            "protons": maxwellian(1, 1.0, 1.0, 0.01, ppc=ppc, seed=1 + 10 * seed),
            "electrons": electrons,
        },
    )
    w_pe, w_pi = plasma_frequency(1, -1, mass_e), plasma_frequency(1, 1, 1.0)
    debye = np.sqrt(Te / mass_e) / w_pe
    S2 = shape_factor2(k, interp_order) * k_over_K(k[0])
    k2 = np.sum(k**2)
    expected = leapfrog(np.sqrt(S2 * (w_pe**2 * (1 + 3 * k2 * debye**2) + w_pi**2)))

    # the longitudinal field, along k
    Ex = mode_amplitude(diag_dir, "E", "x", times, cells, mode)
    Ey = mode_amplitude(diag_dir, "E", "y", times, cells, mode)
    measured, error, damping = standing_fit(
        times, unit[0] * Ex + unit[1] * Ey, guess=expected
    )
    return measured, error, damping, expected


def light(interp_order, seed=0, b=0.05, Te=0.04, ppc=16, polarization="in plane"):
    cells, mode = (32, 32), (1, 1)
    k = wavevector(cells, mode)
    k2 = np.sum(k**2)
    in_plane = polarization == "in plane"

    def b_component(amplitude):
        def fn(x, y):
            return amplitude * np.cos(k[0] * x + k[1] * y)

        return fn

    if in_plane:
        fields = {"bz": b_component(b)}
        # A = (b / |k|^2) (-k_y, k_x) sin(k.x): B_z = dA_y/dx - dA_x/dy
        potential = {"x": -b * k[1] / k2, "y": b * k[0] / k2}
    else:
        unit = k / np.sqrt(k2)
        fields = {"bx": b_component(b * unit[1]), "by": b_component(-b * unit[0])}
        # A_z = (b / |k|) sin(k.x): (B_x, B_y) = (dA_z/dy, -dA_z/dx)
        potential = {"z": b / np.sqrt(k2)}

    populations = {
        "protons": maxwellian(1, 1.0, 1.0, 0.01, ppc=ppc, seed=3 + 10 * seed),
        "electrons": maxwellian(-1, mass_e, 1.0, Te, ppc=ppc, seed=4 + 10 * seed),
    }

    def canonical(q, m, a):
        def fn(x, y):
            return -q / m * a * np.sin(k[0] * x + k[1] * y)

        return fn

    for pop in populations.values():
        for direction, a in potential.items():
            pop[f"vbulk{direction}"] = canonical(pop["charge"], pop["mass"], a)

    name = "light_2d" if in_plane else "light_2d_out_of_plane"
    diag_dir, times = run(
        name, interp_order, cells, final_time=1.0, populations=populations, **fields
    )
    w_pe2 = plasma_frequency(1, -1, mass_e) ** 2
    w_p2 = w_pe2 + plasma_frequency(1, 1, 1.0) ** 2
    S2 = shape_factor2(k, interp_order)
    cold2 = S2 * w_p2 + c_over_va**2 * yee_k2(k)
    warm2 = S2 * w_pe2 * k2 * (Te / mass_e) / cold2
    expected = leapfrog(np.sqrt(cold2 + warm2))

    if in_plane:
        signal = mode_amplitude(diag_dir, "B", "z", times, cells, mode)
    else:
        signal = unit[1] * mode_amplitude(
            diag_dir, "B", "x", times, cells, mode
        ) - unit[0] * mode_amplitude(diag_dir, "B", "y", times, cells, mode)
    measured, error, damping = standing_fit(times, signal, guess=expected)
    return measured, error, damping, expected


def light_out_of_plane(interp_order, seed=0):
    return light(interp_order, seed, polarization="out of plane")


# Weibel: beams +-v0 along z, warm in the plane, cold along z; ions at rest
weibel_v0 = 2.4
weibel_vth_plane = 0.72
weibel_vth_z = 0.05
weibel_vth_ions = 0.1


def weibel_growth(k, interp_order):
    """gamma of the purely growing transverse mode at wavevector k (in the plane)"""
    from scipy.optimize import brentq
    from scipy.special import erfcx

    S2 = shape_factor2(k, interp_order)
    k_norm = np.linalg.norm(k)
    w_pe2 = plasma_frequency(1, -1, mass_e) ** 2
    w_pi2 = plasma_frequency(1, 1, 1.0) ** 2
    anisotropy = (weibel_v0**2 + weibel_vth_z**2) / weibel_vth_plane**2

    def response(gamma, vth):
        y = gamma / (np.sqrt(2) * k_norm * vth)
        return 1 - np.sqrt(np.pi) * y * erfcx(y)

    def dispersion(gamma):
        return (
            -(gamma**2)
            - c_over_va**2 * yee_k2(k)
            - S2 * (w_pe2 + w_pi2)
            + S2 * w_pe2 * anisotropy * response(gamma, weibel_vth_plane)
            + S2 * w_pi2 * response(gamma, weibel_vth_ions)
        )

    assert dispersion(1e-9) > 0, "stable mode"
    gamma = brentq(dispersion, 1e-9, 100.0)
    # the leapfrog grows by (1 + gamma dt / 2 ...) per step: gamma -> (2/dt) asinh(gamma dt/2)
    return 2 / waves.time_step * np.arcsinh(gamma * waves.time_step / 2)


def weibel(interp_order, seed=0, ppc_beam=64, ppc_ions=32, ndim=2):
    """ndim=1 (k along x only) has the same dispersion: with as many particles per beam, it
    measures the noise and saturation levels of the 2D run at a fraction of its cost"""
    cells = (16,) * ndim
    k = wavevector(cells, (1,) + (0,) * (ndim - 1))
    expected = weibel_growth(k, interp_order)

    def beam(sign, beam_seed):
        pop = maxwellian(-1, mass_e, 0.5, 0.0, ppc=ppc_beam, seed=beam_seed)
        pop["vbulkz"] = waves.uniform(sign * weibel_v0)
        pop["vthx"] = pop["vthy"] = waves.uniform(weibel_vth_plane)
        pop["vthz"] = waves.uniform(weibel_vth_z)
        return pop

    diag_dir, times = run(
        f"weibel_{ndim}d",
        interp_order,
        cells,
        final_time=12 / expected,
        populations={
            "ions": maxwellian(
                1, 1.0, 1.0, weibel_vth_ions**2, ppc=ppc_ions, seed=5 + 10 * seed
            ),
            "beam_up": beam(+1, 6 + 10 * seed),
            "beam_down": beam(-1, 7 + 10 * seed),
        },
    )
    # k along x: B_y, k along y: B_x; both modes have the same growth rate
    By = mode_amplitude(diag_dir, "B", "y", times, cells, (1, 0)[:ndim])
    amplitude = np.abs(By)
    if ndim == 2:
        Bx = mode_amplitude(diag_dir, "B", "x", times, cells, (0, 1))
        amplitude = np.sqrt(amplitude**2 + np.abs(Bx) ** 2)
    print(
        f"weibel: |B(1,0)| at t=0 {amplitude[0]:.2e}, max {amplitude.max():.2e}"
        f" at t={times[np.argmax(amplitude)]:.3f}, expected growth {expected:.3f}"
    )
    return exponential_growth(times, amplitude, expected, start_efolds=3), expected


def main():
    with_weibel = "--weibel" in sys.argv
    interp_orders = [int(arg) for arg in sys.argv[1:] if arg != "--weibel"] or [1, 2, 3]
    results = []
    for order in interp_orders:
        for name, case in (
            ("langmuir", langmuir),
            ("light", light),
            ("light out of plane", light_out_of_plane),
        ):
            measured, error, damping, expected = case(order)
            print(
                f"{name} 2D order {order}: fit error {error / expected:.1e},"
                f" damping {damping:.3g}"
            )
            tolerance = TOLERANCES[name.split()[0]]
            results.append(
                check(f"{name} 2D order {order}", measured, expected, tolerance)
            )
            if name.startswith("light"):
                ok = damping / expected < MAX_LIGHT_DAMPING
                print(
                    f"light 2D damping / w = {damping / expected:.2e}"
                    f" {'ok' if ok else 'FAILED'}"
                )
                results.append(ok)
        if with_weibel:
            measured, expected = weibel(order)
            results.append(
                check(
                    f"weibel 2D order {order}", measured, expected, TOLERANCES["weibel"]
                )
            )
    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
