#!/usr/bin/env python3
"""
FullyKineticModel physics checks, 1D, periodic, PHARE units (B0, n0, V_A, d_i, 1/Omega_ci).
A species of density n, charge q and mass m has the plasma frequency c |q| sqrt(n/m), with
c = c/V_A.

- langmuir: electron velocity perturbation v_x = v0 sin(kx), standing Langmuir wave,
  w^2 = w_pe^2 (1 + 3 k^2 lambda_De^2) + w_pi^2 (the ions, m = 1, are not immobile)
- light: unmagnetized plasma, B_z = b cos(kx), standing light wave,
  w^2 = w_pe^2 + w_pi^2 + c^2 k^2
- two_stream: two cold electron beams +-v0 on cold ions, mode 1 seeded, linear growth rate
  from 1 = w_pi^2 / w^2 + sum over beams w_pb^2 / (w -+ k v0)^2

The scheme is compared with: k -> (2/dx) sin(k dx / 2) where k comes from a Yee derivative
(light wave), w -> (2/dt) asin(w dt / 2) for the leapfrog, and w_p^2 -> S(k)^2 w_p^2 for the
particle shape of order p, S(k) = sinc(k dx / 2)^(p+1), which deposits and gathers (Birdsall
and Langdon, Plasma Physics via Computer Simulation, ch. 8). Measured at k dx = 0.2: without
S(k) the Langmuir error grows with the order, -0.71%, -0.92%, -1.03% for p = 1, 2, 3. The dumped E is half a step
behind its time label; frequencies and growth rates do not depend on it.
"""

import sys

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

c_over_va = 10.0
mass_e = 0.04  # m_i / m_e = 25
dl = 0.04
time_step = 0.002  # w_pe dt = 0.1, c dt / dx = 0.5
diag_root = "phare_outputs/fully_kinetic_waves"

FREQUENCY_TOLERANCE = 1e-2
GROWTH_TOLERANCE = 1e-1


def plasma_frequency(n, q, m):
    return c_over_va * abs(q) * np.sqrt(n / m)


def shape_factor2(k, interp_order):
    return np.sinc(k * dl / 2 / np.pi) ** (2 * (interp_order + 1))


def leapfrog(omega):
    return 2 / time_step * np.arcsin(omega * time_step / 2)


def uniform(value):
    def fn(x):
        return value + 0.0 * x

    return fn


def maxwellian(charge, mass, density, T, vx=None, ppc=200, seed=0):
    vth = uniform(np.sqrt(T / mass))
    return {
        "charge": charge,
        "mass": mass,
        "density": uniform(density),
        "vbulkx": vx or uniform(0.0),
        "vbulky": uniform(0.0),
        "vbulkz": uniform(0.0),
        "vthx": vth,
        "vthy": vth,
        "vthz": vth,
        "nbr_part_per_cell": ppc,
        "init": {"seed": seed},
    }


def run(name, interp_order, cells, final_time, populations, by=None, bz=None, bx=0.0):
    diag_dir = f"{diag_root}/{name}_{interp_order}"
    nsteps = round(final_time / time_step)
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=nsteps,
        cells=cells,
        dl=dl,
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        # these runs explore warned regimes on purpose (cold beams, v_th,e / c = 0.3)
        strict=False,
    )
    ph.MaxwellianFluidModel(
        bx=uniform(bx), by=by or uniform(0.0), bz=bz or uniform(0.0), **populations
    )
    timestamps = time_step * np.arange(0, nsteps + 1, 2)
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    Simulator(sim).run().reset()
    return diag_dir, timestamps


def mode_amplitude(diag_dir, quantity, component, times, cells, mode):
    """complex Fourier coefficient of a field component at k = 2 pi mode / L, per time"""
    L = cells * dl
    amplitude = []
    for time in times:
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/EM_{quantity}.h5", times=f"{time:.10f}"
        )
        values = {}
        for patch in hier.level(0, time).patches:
            pdata = patch.patch_datas[f"{quantity}{component}"]
            for x, v in zip(pdata.x, pdata.dataset[:]):
                if -1e-6 < x < L - 1e-6:
                    values[round(2 * x / dl)] = (x, v)
        x, v = np.array([values[k] for k in sorted(values)]).T
        amplitude.append(np.mean(v * np.exp(-2j * np.pi * mode * x / L)))
    return np.array(amplitude)


def standing_frequency(times, amplitude, guess):
    """frequency of a standing oscillation of a complex mode amplitude: least squares fit of
    a cos(w t + p) + c to its projection on its dominant phase. Zero crossings are biased
    once a damped wave reaches the noise level."""
    from scipy.optimize import curve_fit

    phase = np.angle(amplitude[np.argmax(np.abs(amplitude))])
    signal = (amplitude * np.exp(-1j * phase)).real

    def model(t, a, w, p, c):
        return a * np.cos(w * t + p) + c

    best = None
    for p0 in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        try:
            params, _ = curve_fit(
                model, times, signal, p0=[np.abs(signal).max(), guess, p0, 0]
            )
        except RuntimeError:
            continue
        residual = np.sum((model(times, *params) - signal) ** 2)
        if best is None or residual < best[0]:
            best = (residual, abs(params[1]))
    assert best is not None, "no frequency fit converged"
    return best[1]


def check(name, measured, expected, tolerance):
    error = abs(measured - expected) / abs(expected)
    status = "ok" if error < tolerance else "FAILED"
    print(
        f"{name}: measured {measured:.4f} expected {expected:.4f} error {error:.2e} {status}"
    )
    return error < tolerance


def langmuir(interp_order):
    cells, mode, Te, v0 = 32, 1, 0.04, 0.3
    k = 2 * np.pi * mode / (cells * dl)

    def vx(x):
        return v0 * np.sin(k * x)

    diag_dir, times = run(
        "langmuir",
        interp_order,
        cells,
        final_time=20 * 2 * np.pi / plasma_frequency(1, 1, mass_e),
        populations={
            "protons": maxwellian(1, 1.0, 1.0, 0.01, ppc=400, seed=1),
            "electrons": maxwellian(-1, mass_e, 1.0, Te, vx=vx, ppc=400, seed=2),
        },
    )
    w_pe, w_pi = plasma_frequency(1, -1, mass_e), plasma_frequency(1, 1, 1.0)
    debye = np.sqrt(Te / mass_e) / w_pe
    S2 = shape_factor2(k, interp_order)
    expected = leapfrog(np.sqrt(S2 * (w_pe**2 * (1 + 3 * (k * debye) ** 2) + w_pi**2)))

    Ex = mode_amplitude(diag_dir, "E", "x", times, cells, mode)
    measured = standing_frequency(times, Ex, guess=expected)
    return check(
        f"langmuir order {interp_order}", measured, expected, FREQUENCY_TOLERANCE
    )


def light(interp_order):
    cells, mode, b = 64, 2, 0.05
    k = 2 * np.pi * mode / (cells * dl)

    def bz(x):
        return b * np.cos(k * x)

    diag_dir, times = run(
        "light",
        interp_order,
        cells,
        final_time=1.0,
        populations={
            "protons": maxwellian(1, 1.0, 1.0, 0.01, ppc=100, seed=3),
            "electrons": maxwellian(-1, mass_e, 1.0, 0.04, ppc=100, seed=4),
        },
        bz=bz,
    )
    k_grid = 2 / dl * np.sin(k * dl / 2)
    w_p2 = plasma_frequency(1, -1, mass_e) ** 2 + plasma_frequency(1, 1, 1.0) ** 2
    S2 = shape_factor2(k, interp_order)
    expected = leapfrog(np.sqrt(S2 * w_p2 + (c_over_va * k_grid) ** 2))

    Bz = mode_amplitude(diag_dir, "B", "z", times, cells, mode)
    measured = standing_frequency(times, Bz, guess=expected)
    return check(f"light order {interp_order}", measured, expected, FREQUENCY_TOLERANCE)


def two_stream_growth(k, v0, interp_order):
    """largest imaginary part of w, roots of the cold two-beam + ions dispersion relation"""
    S2 = shape_factor2(k, interp_order)
    w_pb2 = S2 * plasma_frequency(0.5, -1, mass_e) ** 2
    w_pi2 = S2 * plasma_frequency(1, 1, 1.0) ** 2
    a, b = k * v0, w_pb2
    # w^2 (w-a)^2 (w+a)^2 = w_pi^2 (w-a)^2 (w+a)^2 + b w^2 ((w+a)^2 + (w-a)^2)
    w = np.polynomial.Polynomial([0, 1])
    lhs = w**2 * (w - a) ** 2 * (w + a) ** 2
    rhs = w_pi2 * (w - a) ** 2 * (w + a) ** 2 + b * w**2 * ((w + a) ** 2 + (w - a) ** 2)
    return max(r.imag for r in (lhs - rhs).roots())


def two_stream(interp_order):
    cells, mode, v0 = 32, 1, 6.0
    k = 2 * np.pi * mode / (cells * dl)

    def beam(sign):
        def vx(x):
            return sign * v0 * (1 + 1e-3 * np.sin(k * x))

        return vx

    expected = two_stream_growth(k, v0, interp_order)
    diag_dir, times = run(
        "two_stream",
        interp_order,
        cells,
        final_time=12 / expected,
        populations={
            "protons": maxwellian(1, 1.0, 1.0, 1e-4, ppc=200, seed=5),
            "beam_right": maxwellian(
                -1, mass_e, 0.5, 1e-4, vx=beam(+1), ppc=200, seed=6
            ),
            "beam_left": maxwellian(
                -1, mass_e, 0.5, 1e-4, vx=beam(-1), ppc=200, seed=7
            ),
        },
    )
    amplitude = np.abs(mode_amplitude(diag_dir, "E", "x", times, cells, mode))

    # linear phase: between 30 times the initial level and a tenth of the maximum
    peak = np.argmax(amplitude)
    window = (amplitude[:peak] > 30 * amplitude[0]) & (
        amplitude[:peak] < 0.1 * amplitude[peak]
    )
    assert window.sum() >= 10, f"linear phase too short ({window.sum()} dumps)"
    measured = np.polyfit(times[:peak][window], np.log(amplitude[:peak][window]), 1)[0]
    return check(
        f"two-stream order {interp_order}", measured, expected, GROWTH_TOLERANCE
    )


def main():
    interp_orders = [int(arg) for arg in sys.argv[1:]] or [1, 2, 3]
    results = [
        case(order) for order in interp_orders for case in (langmuir, light, two_stream)
    ]
    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
