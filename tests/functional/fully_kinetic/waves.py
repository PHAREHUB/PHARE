#!/usr/bin/env python3
"""
FullyKineticModel physics checks, 1D, periodic, PHARE units (B0, n0, V_A, d_i, 1/Omega_ci).
A species of density n, charge q and mass m has the plasma frequency c |q| sqrt(n/m), with
c = c/V_A.

- langmuir: electron velocity perturbation v_x = v0 sin(kx), standing Langmuir wave,
  w^2 = w_pe^2 (1 + 3 k^2 lambda_De^2) + w_pi^2 (the ions, m = 1, are not immobile)
- light: unmagnetized plasma, B_z = b cos(kx), standing light wave,
  w^2 = w_pe^2 + w_pi^2 + c^2 k^2. Particles start with v_y = -(q/m) A_y, A_y = (b/k) sin(kx):
  at rest they would keep the memory of the initial field (conserved canonical momentum), a
  static B_z of w_p^2/w^2 ~ half of b that thermal motion mixes away, which a damped fit
  reads as damping (measured: 0.064 at rest, 0.015 with this start)
- two_stream: two cold electron beams +-v0 on cold ions, mode 1 seeded, linear growth rate
  from 1 = w_pi^2 / w^2 + sum over beams w_pb^2 / (w -+ k v0)^2

The scheme is compared with:
- k -> (2/dx) sin(k dx / 2) where k comes from a Yee derivative (light wave),
- w -> (2/dt) asin(w dt / 2) for the leapfrog,
- w_p^2 -> S(k)^2 w_p^2 for the particle shape of order p, S(k) = sinc(k dx / 2)^(p+1), which
  deposits and gathers (Birdsall and Langdon, Plasma Physics via Computer Simulation, ch. 8).
  Measured at k dx = 0.2: without S(k) the Langmuir error grows with the order, -0.71%, -0.92%,
  -1.03% for p = 1, 2, 3,
- w_p^2 -> w_p^2 k / K for longitudinal waves: the discrete Gauss law divergence acts as
  K = (2/dx) sin(k dx / 2) while the particle density responds to k (+0.08% on w),
- the warm transverse correction w_pe^2 -> w_pe^2 (1 + k^2 v_th,x^2 / w^2) (+0.12% on w).
Neglected: higher thermal terms (~1e-4) and the nonlinear Langmuir shift, which grows as v0^2
(measured -0.49% +- 0.14% at v0 = 0.3, -0.18% +- 0.14% at v0 = 0.15). The dumped E is half a step
behind its time label; frequencies and growth rates do not depend on it.

Measured numbers: Debug builds, recorded with the code of these commits: S(k) order errors
86fdc49b2; 5-seed tolerances, nonlinear Langmuir shift, the +0.08% and +0.12% corrections and
the x155 seed projection ea8a5f97e; light start (0.064 / 0.015), light damping and window
curvature 524f91e6d.
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

# |mean| + 3 standard deviations of a single run, measured over 5 random seeds at order 1:
# langmuir -0.18% +- 0.14%, light -0.04% +- 0.09%, two-stream growth +1.1% +- 1.1%
TOLERANCES = {"langmuir": 7e-3, "light": 4e-3, "two-stream": 5e-2}
# light wave damping, 5 seeds: gamma / w = 0.9e-4 +- 1.3e-4, consistent with no damping
# (it was 1.1e-3 with particles starting at rest, see the light case above)
MAX_LIGHT_DAMPING = 1e-3
# two-stream: growth rates fitted on the two halves of the linear window differ by up to 10%
# (5 seeds); more means the window reached the nonlinear phase
MAX_WINDOW_CURVATURE = 0.15


def plasma_frequency(n, q, m):
    return c_over_va * abs(q) * np.sqrt(n / m)


def shape_factor2(k, interp_order):
    """S(k)^2 of the particle shape (deposit and gather), product over the directions of k
    (a number in 1D, one component per direction otherwise)"""
    return np.prod(
        np.sinc(np.atleast_1d(k) * dl / 2 / np.pi) ** (2 * (interp_order + 1))
    )


def yee_k(k):
    """K_d = (2/dx) sin(k_d dx / 2), the wavenumber a Yee derivative sees, per direction"""
    return 2 / dl * np.sin(np.atleast_1d(k) * dl / 2)


def yee_k2(k):
    return np.sum(yee_k(k) ** 2)


def filter_transfer(k, passes):
    """transfer of the binomial current filter, prod_d cos^2(k_d dx / 2) per pass: it
    multiplies the plasma response (the current is filtered, the gather is not)"""
    return np.prod(np.cos(np.atleast_1d(k) * dl / 2) ** (2 * passes))


def wavevector(cells, mode):
    """k = 2 pi mode / L, per direction"""
    return 2 * np.pi * np.atleast_1d(mode) / (np.atleast_1d(cells) * dl)


def leapfrog(omega):
    return 2 / time_step * np.arcsin(omega * time_step / 2)


def uniform(value):
    def fn(*xyz):
        return value + 0.0 * xyz[0]

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


def run(
    name,
    interp_order,
    cells,
    final_time,
    populations,
    by=None,
    bz=None,
    bx=0.0,
    current_filter_passes=0,
    divergence_cleaning=0.0,
):
    diag_dir = f"{diag_root}/{name}_{interp_order}"
    nsteps = round(final_time / time_step)
    cells = list(np.atleast_1d(cells))
    ndim = len(cells)
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=nsteps,
        cells=cells,
        dl=[dl] * ndim,
        interp_order=interp_order,
        model_options=["FullyKineticModel"],
        c_over_va=c_over_va,
        current_filter_passes=current_filter_passes,
        divergence_cleaning=divergence_cleaning,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        # these runs explore warned regimes on purpose (cold beams, v_th,e / c = 0.3)
        strict=False,
    )
    ph.MaxwellianFluidModel(
        bx=bx if callable(bx) else uniform(bx),
        by=by or uniform(0.0),
        bz=bz or uniform(0.0),
        **populations,
    )
    timestamps = time_step * np.arange(0, nsteps + 1, 2)
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    Simulator(sim).run().reset()
    return diag_dir, timestamps


def mode_amplitude(diag_dir, quantity, component, times, cells, mode):
    """complex Fourier coefficient of a field component at k = 2 pi mode / L, per time.
    cells and mode are numbers in 1D, sequences of one value per direction otherwise."""
    L = np.atleast_1d(cells) * dl
    k = 2 * np.pi * np.atleast_1d(mode) / L
    amplitude = []
    for time in times:
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/EM_{quantity}.h5", times=f"{time:.10f}"
        )
        values = {}
        for patch in hier.level(0, time).patches:
            pdata = patch.patch_datas[f"{quantity}{component}"]
            mesh = pdata.meshgrid()
            inside = np.all(
                [(x > -1e-6) & (x < length - 1e-6) for x, length in zip(mesh, L)],
                axis=0,
            )
            positions = np.stack([x[inside] for x in mesh], axis=1)
            keys = np.round(2 * positions / dl).astype(int)
            for key, position, v in zip(keys, positions, pdata.dataset[:][inside]):
                values[tuple(key)] = (position, v)
        position = np.array([values[key][0] for key in sorted(values)])
        v = np.array([values[key][1] for key in sorted(values)])
        amplitude.append(np.mean(v * np.exp(-1j * position @ k)))
    return np.array(amplitude)


def standing_fit(times, amplitude, guess):
    """fit a exp(-g t) cos(w t + p) + c to the projection of a complex mode amplitude on its
    dominant phase: returns w, its standard error, and g. Zero crossings are biased once a
    damped wave reaches the noise level, and a constant amplitude fit mismatches its decay."""
    from scipy.optimize import curve_fit

    phase = np.angle(amplitude[np.argmax(np.abs(amplitude))])
    signal = (amplitude * np.exp(-1j * phase)).real

    def model(t, a, w, p, c, g):
        return a * np.exp(-g * t) * np.cos(w * t + p) + c

    best = None
    for p0 in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        try:
            params, cov = curve_fit(
                model, times, signal, p0=[np.abs(signal).max(), guess, p0, 0, 0]
            )
        except RuntimeError:
            continue
        residual = np.sum((model(times, *params) - signal) ** 2)
        if best is None or residual < best[0]:
            best = (residual, abs(params[1]), np.sqrt(cov[1, 1]), params[4])
    assert best is not None, "no frequency fit converged"
    return best[1:]


def check(name, measured, expected, tolerance):
    error = abs(measured - expected) / abs(expected)
    status = "ok" if error < tolerance else "FAILED"
    print(
        f"{name}: measured {measured:.4f} expected {expected:.4f} error {error:.2e} {status}"
    )
    return error < tolerance


def k_over_K(k):
    """longitudinal fields: in the discrete Gauss law the grid divergence acts as
    (2/dx) sin(k dx/2) = K, while the particle density responds to k: w_p^2 -> w_p^2 k / K"""
    return 1 / np.sinc(k * dl / 2 / np.pi)


def langmuir(interp_order, seed=0, v0=0.15, Te=0.04, ppc=400, divergence_cleaning=0.0):
    """divergence_cleaning (Marder, s): must not change the resolved physics. Measured, order 1,
    seed 0, Debug, C++ at c757493ea: -1.33e-3 without, -1.53e-3 with s = 0.5 (fit error
    2.6e-5)"""
    cells, mode = 32, 1
    k = 2 * np.pi * mode / (cells * dl)

    def vx(x):
        return v0 * np.sin(k * x)

    diag_dir, times = run(
        f"langmuir_cleaning{divergence_cleaning}"
        if divergence_cleaning
        else "langmuir",
        interp_order,
        cells,
        final_time=20 * 2 * np.pi / plasma_frequency(1, 1, mass_e),
        populations={
            "protons": maxwellian(1, 1.0, 1.0, 0.01, ppc=ppc, seed=1 + 10 * seed),
            "electrons": maxwellian(
                -1, mass_e, 1.0, Te, vx=vx, ppc=ppc, seed=2 + 10 * seed
            ),
        },
        divergence_cleaning=divergence_cleaning,
    )
    w_pe, w_pi = plasma_frequency(1, -1, mass_e), plasma_frequency(1, 1, 1.0)
    debye = np.sqrt(Te / mass_e) / w_pe
    S2 = shape_factor2(k, interp_order) * k_over_K(k)
    expected = leapfrog(np.sqrt(S2 * (w_pe**2 * (1 + 3 * (k * debye) ** 2) + w_pi**2)))

    Ex = mode_amplitude(diag_dir, "E", "x", times, cells, mode)
    measured, error, damping = standing_fit(times, Ex, guess=expected)
    return measured, error, damping, expected


def light(interp_order, seed=0, b=0.05, Te=0.04, ppc=100, density=1.0):
    cells, mode = 64, 2
    k = 2 * np.pi * mode / (cells * dl)

    def bz(x):
        return b * np.cos(k * x)

    populations = {
        "protons": maxwellian(1, 1.0, density, 0.01, ppc=ppc, seed=3 + 10 * seed),
        "electrons": maxwellian(-1, mass_e, density, Te, ppc=ppc, seed=4 + 10 * seed),
    }

    def canonical_vy(q, m):
        def vy(x):
            return -q / m * b / k * np.sin(k * x)

        return vy

    for pop in populations.values():
        pop["vbulky"] = canonical_vy(pop["charge"], pop["mass"])

    diag_dir, times = run(
        "light", interp_order, cells, final_time=1.0, populations=populations, bz=bz
    )
    w_pe2 = plasma_frequency(density, -1, mass_e) ** 2
    w_p2 = w_pe2 + plasma_frequency(density, 1, 1.0) ** 2
    S2 = shape_factor2(k, interp_order)
    cold2 = S2 * w_p2 + c_over_va**2 * yee_k2(k)
    # warm transverse wave: w_pe^2 -> w_pe^2 (1 + k^2 v_th,x^2 / w^2)
    warm2 = S2 * w_pe2 * k**2 * (Te / mass_e) / cold2
    expected = leapfrog(np.sqrt(cold2 + warm2))

    Bz = mode_amplitude(diag_dir, "B", "z", times, cells, mode)
    measured, error, damping = standing_fit(times, Bz, guess=expected)
    return measured, error, damping, expected


# heavy neutralizing ions: with m_i = 1, k v0 / w_pb = 0.83 sits where an ion-coupled root and
# the electron two-stream root exchange the largest growth (gamma jumps from 14.1 to 11.8
# between 0.8 and 0.9), the expected value is ill-conditioned. With m_i = 100, w_pi = 1.
two_stream_ion_mass = 100.0


def two_stream_growth(k, v0, interp_order):
    """largest imaginary part of w, roots of the cold two-beam + ions dispersion relation"""
    S2 = shape_factor2(k, interp_order) * k_over_K(k)
    w_pb2 = S2 * plasma_frequency(0.5, -1, mass_e) ** 2
    w_pi2 = S2 * plasma_frequency(1, 1, two_stream_ion_mass) ** 2
    a, b = k * v0, w_pb2
    # w^2 (w-a)^2 (w+a)^2 = w_pi^2 (w-a)^2 (w+a)^2 + b w^2 ((w+a)^2 + (w-a)^2)
    w = np.polynomial.Polynomial([0, 1])
    lhs = w**2 * (w - a) ** 2 * (w + a) ** 2
    rhs = w_pi2 * (w - a) ** 2 * (w + a) ** 2 + b * w**2 * ((w + a) ** 2 + (w - a) ** 2)
    return max(r.imag for r in (lhs - rhs).roots())


def two_stream(interp_order, seed=0):
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
            "ions": maxwellian(
                1, two_stream_ion_mass, 1.0, 1e-4, ppc=200, seed=5 + 10 * seed
            ),
            "beam_right": maxwellian(
                -1, mass_e, 0.5, 1e-4, vx=beam(+1), ppc=200, seed=6 + 10 * seed
            ),
            "beam_left": maxwellian(
                -1, mass_e, 0.5, 1e-4, vx=beam(-1), ppc=200, seed=7 + 10 * seed
            ),
        },
    )
    amplitude = np.abs(mode_amplitude(diag_dir, "E", "x", times, cells, mode))
    # the velocity-only seed still projects onto the eigenmode in the first e-foldings:
    # amplitude x155 in the first 2, the decaying root is gone after
    return exponential_growth(times, amplitude, expected, start_efolds=3), expected


def exponential_growth(times, amplitude, expected, start_efolds):
    """growth rate fitted on the linear phase of |amplitude|: from start_efolds / expected
    until the amplitude first exceeds 30% of saturation. Saturation is the first approach to
    the maximum: later bounces can exceed it."""
    peak = np.argmax(amplitude)
    saturation = amplitude[peak]
    assert peak < len(amplitude) - 10, (
        "saturation not reached before the end of the run"
    )
    assert amplitude[peak:].min() < 0.8 * saturation, "no saturation, still growing"
    end = np.argmax(amplitude > 0.3 * saturation)
    window = slice(np.argmax(times > start_efolds / expected), end)
    t, log_a = times[window], np.log(amplitude[window])
    assert len(t) >= 8, f"linear phase too short ({len(t)} dumps)"
    measured = np.polyfit(t, log_a, 1)[0]
    half = len(t) // 2
    lower = np.polyfit(t[:half], log_a[:half], 1)[0]
    upper = np.polyfit(t[half:], log_a[half:], 1)[0]
    print(f"growth on the lower/upper half of the window: {lower:.3f} / {upper:.3f}")
    assert abs(upper - lower) / measured < MAX_WINDOW_CURVATURE, (
        "window is not exponential"
    )
    return measured


def main():
    interp_orders = [int(arg) for arg in sys.argv[1:]] or [1, 2, 3]
    results = []
    for order in interp_orders:
        for name, case in (("langmuir", langmuir), ("light", light)):
            measured, error, damping, expected = case(order)
            print(
                f"{name} order {order}: fit error {error / expected:.1e}, damping {damping:.3g}"
            )
            results.append(
                check(f"{name} order {order}", measured, expected, TOLERANCES[name])
            )
            if name == "light":
                ok = damping / expected < MAX_LIGHT_DAMPING
                print(
                    f"light damping / w = {damping / expected:.2e} {'ok' if ok else 'FAILED'}"
                )
                results.append(ok)
        measured, expected = two_stream(order)
        results.append(
            check(
                f"two-stream order {order}",
                measured,
                expected,
                TOLERANCES["two-stream"],
            )
        )
    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
