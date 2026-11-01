#!/usr/bin/env python3
"""
Parallel hybrid dispersion relation: Alfven/ion-cyclotron (AIC) and whistler branches.

Cold ions, cold massless electrons, units di and 1/Omega_ci. With b = By + i Bz and
u = vy + i vz, the linear system is

    du/dt = db/dx,    db/dt = du/dx - i d2b/dx2

so a mode b ~ exp(i(qx - wt)) obeys w^2 + q^2 w - q^2 = 0, i.e.

    w_aic(q) = q^2/2 (sqrt(1 + 4/q^2) - 1)   (w > 0)
    w_wh(q)  = q^2/2 (sqrt(1 + 4/q^2) + 1)   (w < 0 root, |w| given here)

with eigenvector u = -(q/w) b. Each run loads pure AIC modes at q = +m k0 and pure
whistler modes at q = -m k0, so every spatial Fourier bin holds one frequency, read
from the slope of its phase in time. Two runs cover k in [0.06, 0.5] and [1, 4].

The scheme should not give exactly w(k), and measured w are compared to w_scheme:
- on the Yee grid d/dx acts on a mode as i k_grid, k_grid = (2/dx) sin(k dx/2);
- the order-p particle shape (B-spline of variance (p+1) dx^2/12) weakens the Hall
  term as if di = 1 - (p+1) dx^2/12. Empirical: measured for p = 1 at dx = 0.2,
  0.4, 0.8 and for p = 1, 2, 3 at dx = 0.4; it accounts for a Hall splitting
  w_wh - w_aic up to 11% too small and leaves errors on w of 0.1% or less.
With di, w_aic/wh = k (sqrt(1 + k^2 di^2/4) -/+ k di/2), and the eigenvector is unchanged.
Every mode is resolved by at least 16 cells per wavelength.
"""

import os
import sys

import matplotlib as mpl
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator, startMPI

mpl.use("Agg")


# dB/B0 of each wave. Must stay linear: at 3e-2 the AIC waves get a +2.8% frequency shift
# and decay near Omega_ci; noise is fought with particles per cell instead.
AMPLITUDE = 1e-2
INTERP_ORDERS = (1, 2, 3)

# |w_measured - w_scheme| / w_scheme, and the allowed final/initial mode amplitude ratio
TOLERANCE = 1e-2
AMPLITUDE_RANGE = (0.3, 3.0)

RUNS = {
    "low_k": dict(
        cells=128,
        dl=0.8,
        modes=(1, 2, 4, 8),
        ppc=100,
        time_step=0.02,
        final_time=100.0,
        dump_every=10,
    ),
    "high_k": dict(
        cells=64,
        dl=0.1,
        modes=(1, 2, 3, 4),
        ppc=400,
        time_step=0.001,
        final_time=10.0,
        dump_every=10,
    ),
}


def omega_aic(k, di=1.0):
    return k * (np.sqrt(1 + (k * di) ** 2 / 4) - k * di / 2)


def omega_whistler(k, di=1.0):
    return k * (np.sqrt(1 + (k * di) ** 2 / 4) + k * di / 2)


def k_grid(k, dl):
    return 2 / dl * np.sin(k * dl / 2)


def di_scheme(dl, interp_order):
    return 1 - (interp_order + 1) * dl**2 / 12


def wave_modes(lx, dl, interp_order, modes):
    """(m, q, w_scheme) of every loaded mode: AIC at +q (w > 0), whistler at -q (w < 0)."""
    k0 = 2 * np.pi / lx
    di = di_scheme(dl, interp_order)
    aic = [(m, m * k0, omega_aic(k_grid(m * k0, dl), di)) for m in modes]
    whistler = [(m, -m * k0, -omega_whistler(k_grid(m * k0, dl), di)) for m in modes]
    return aic, whistler


def diag_dir(name, interp_order):
    return f"dispersion_{name}_o{interp_order}"


def config(
    name,
    interp_order,
    cells,
    dl,
    modes,
    ppc,
    time_step,
    final_time,
    dump_every,
    seed=None,
):
    aic, whistler = wave_modes(cells * dl, dl, interp_order, modes)
    phases = np.random.default_rng(0).uniform(0, 2 * np.pi, 2 * len(modes))
    waves = [(q, w, phi) for (_, q, w), phi in zip(aic + whistler, phases)]

    def b(x):
        return sum(AMPLITUDE * np.exp(1j * (q * x + phi)) for q, w, phi in waves)

    def u(x):
        return sum(
            -np.sign(q)
            * k_grid(abs(q), dl)
            / w
            * AMPLITUDE
            * np.exp(1j * (q * x + phi))
            for q, w, phi in waves
        )

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        interp_order=interp_order,
        resistivity=0.0,
        hyper_resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir(name, interp_order), "mode": "overwrite"},
        },
    )

    init = {"seed": seed} if seed is not None else {}
    ph.MaxwellianFluidModel(
        bx=lambda x: 1.0 + 0 * x,
        by=lambda x: b(x).real,
        bz=lambda x: b(x).imag,
        protons={
            "charge": 1,
            "nbr_part_per_cell": ppc,
            "density": lambda x: 1.0 + 0 * x,
            "vbulkx": lambda x: 0.0 * x,
            "vbulky": lambda x: u(x).real,
            "vbulkz": lambda x: u(x).imag,
            "vthx": lambda x: 1e-3 + 0 * x,
            "vthy": lambda x: 1e-3 + 0 * x,
            "vthz": lambda x: 1e-3 + 0 * x,
            "init": init,
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    dump_dt = time_step * dump_every
    nbr_dumps = round(final_time / dump_dt)
    ph.ElectromagDiagnostics(
        quantity="B", write_timestamps=dump_dt * np.arange(nbr_dumps + 1)
    )
    return sim


def read_b(diag_dir):
    """Times, cell-centred x and complex b = By + i Bz on the domain (level 0)."""
    import h5py

    times, bs = [], []
    with h5py.File(os.path.join(diag_dir, "EM_B.h5"), "r") as h5:
        for time_str in sorted(h5["t"].keys(), key=float):
            patches = h5["t"][time_str]["pl0"]
            xs, by, bz = [], [], []
            for patch in patches.values():
                lower, ncells = patch.attrs["lower"][0], patch.attrs["nbrCells"][0]
                by_p = patch["EM_B_y"][:]
                nghost = (
                    by_p.size - ncells
                ) // 2  # By, Bz are dual: one value per cell
                sl = slice(nghost, nghost + ncells)
                xs.append(lower + np.arange(ncells) + 0.5)
                by.append(by_p[sl])
                bz.append(patch["EM_B_z"][:][sl])
            x = np.concatenate(xs)
            order = np.argsort(x)
            times.append(float(time_str))
            bs.append((np.concatenate(by) + 1j * np.concatenate(bz))[order])
            x_sorted = x[order]
    return np.asarray(times), x_sorted, np.asarray(bs)


def measured_omega(times, x, b, q):
    """Frequency of the Fourier bin q from the slope of its unwrapped phase."""
    bq = (b * np.exp(-1j * q * x)[None, :]).mean(axis=1)
    phase = np.unwrap(np.angle(bq))
    slope = np.polyfit(times, phase, 1)[0]
    return -slope, np.abs(bq)


def analyse(name, interp_order, cells, dl, modes, **_):
    """One row per mode: (run, branch, m, k, w_theory, w_measured, amp_end/amp_start, w_scheme)."""
    times, xi, b = read_b(diag_dir(name, interp_order))
    x = xi * dl
    aic, whistler = wave_modes(cells * dl, dl, interp_order, modes)
    rows = []
    for branch, omega, waves in (
        ("aic", omega_aic, aic),
        ("whistler", omega_whistler, whistler),
    ):
        for m, q, w_scheme in waves:
            w_num, amp = measured_omega(times, x, b, q)
            k = abs(q)
            rows.append(
                (
                    name,
                    branch,
                    m,
                    k,
                    omega(k),
                    abs(w_num),
                    amp[-1] / amp[0],
                    abs(w_scheme),
                )
            )
    return rows


def plot(rows, interp_order, filename):
    import matplotlib.pyplot as plt

    branches = {
        "whistler": ("Whistler (R)", omega_whistler, "tab:red"),
        "aic": ("Alfvén/ion-cyclotron (L)", omega_aic, "tab:green"),
    }
    markers = {"low_k": "o", "high_k": "s"}

    k = np.logspace(np.log10(0.04), np.log10(12), 300)
    fig, (ax, ax_err) = plt.subplots(
        2, 1, figsize=(8, 7.5), sharex=True, height_ratios=(3, 1.3)
    )

    for label, omega, color in branches.values():
        ax.loglog(k, omega(k), color=color, lw=1.5, label=f"{label}, cold theory")
    ax.axhline(1.0, color="gray", ls=":", lw=1)
    ax.text(0.045, 1.08, r"$\Omega_{ci}$", color="gray")

    for run, marker in markers.items():
        if not any(r[0] == run for r in rows):
            continue
        dl = RUNS[run]["dl"]
        for branch, (_, _, color) in branches.items():
            pts = [r for r in rows if r[0] == run and r[1] == branch]
            k_pts = [r[3] for r in pts]
            ax.scatter(
                k_pts,
                [r[5] for r in pts],
                marker=marker,
                s=40,
                zorder=3,
                facecolor=color,
                edgecolor="black",
                lw=0.6,
            )
            ax_err.scatter(
                k_pts,
                [100 * (r[5] - r[7]) / r[7] for r in pts],
                marker=marker,
                s=40,
                zorder=3,
                facecolor=color,
                edgecolor="black",
                lw=0.6,
            )
            for r in pts:
                ax.annotate(
                    f"m={r[2]}",
                    (r[3], r[5]),
                    textcoords="offset points",
                    xytext=(5, 7 if branch == "whistler" else -12),
                    fontsize=7,
                )
        # one legend entry per run, in black, so marker shape reads as "which run"
        ax.scatter(
            [],
            [],
            marker=marker,
            facecolor="white",
            edgecolor="black",
            label=f"PHARE, {RUNS[run]['cells']} cells, $\\Delta x$ = {dl} $d_i$",
        )

    ax.set_ylabel(r"$\omega\,/\,\Omega_{ci}$")
    ax.set_title(
        f"PHARE hybrid PIC — parallel dispersion relation, interpolation order {interp_order}\n"
        r"cold ions ($v_{th}=10^{-3}\,v_A$), $T_e=0$, $B_0\parallel x$, 1D, $\geq 16$ cells per wavelength",
        fontsize=11,
    )
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    ax_err.axhline(0.0, color="black", lw=0.8)
    ax_err.set_xscale("log")
    ax_err.set_xlabel(r"$k\,d_i$")
    ax_err.set_ylabel(r"$(\omega_{PHARE}-\omega_{scheme})/\omega_{scheme}$ [%]")
    ax_err.set_title(
        r"error vs $\omega_{scheme}=\omega(k_{grid}, d_i^{eff})$, "
        r"$k_{grid}=\frac{2}{\Delta x}\sin\frac{k\Delta x}{2}$, "
        r"$d_i^{eff}=1-\frac{(p+1)\Delta x^2}{12}$ (particle shape order $p$)",
        fontsize=9,
    )
    ax_err.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(filename, dpi=130)


def check(rows):
    """Failures, one line each, for modes off the grid dispersion relation or lost."""
    failures = []
    for run, branch, m, _, _, w_num, amp, w_scheme in rows:
        err = (w_num - w_scheme) / w_scheme
        if abs(err) > TOLERANCE:
            failures.append(
                f"{run} {branch} m={m}: w={w_num:.4f}, expected {w_scheme:.4f} ({err:+.2%})"
            )
        if not AMPLITUDE_RANGE[0] < amp < AMPLITUDE_RANGE[1]:
            failures.append(
                f"{run} {branch} m={m}: amplitude ratio {amp:.2f} out of {AMPLITUDE_RANGE}"
            )
    return failures


def main():
    startMPI()
    names = sys.argv[1:] or list(RUNS)
    failures = []
    for interp_order in INTERP_ORDERS:
        rows = []
        for name in names:
            Simulator(config(name, interp_order, **RUNS[name])).run().reset()
            ph.global_vars.sim = None
            if cpp.mpi_rank() == 0:
                rows += analyse(name, interp_order, **RUNS[name])

        if cpp.mpi_rank() == 0:
            print(f"interpolation order {interp_order}")
            print(
                f"{'run':7} {'branch':9} {'m':>2} {'k':>7} {'w_th':>9} {'w_scheme':>9} {'w_num':>9} "
                f"{'err_th':>9} {'err_grid':>9} {'amp':>6}"
            )
            for run, branch, m, k, w, w_num, amp, w_scheme in rows:
                print(
                    f"{run:7} {branch:9} {m:2d} {k:7.4f} {w:9.6f} {w_scheme:9.6f} {w_num:9.6f} "
                    f"{(w_num - w) / w:9.2e} {(w_num - w_scheme) / w_scheme:9.2e} {amp:6.2f}"
                )
            plot(rows, interp_order, f"dispersion_o{interp_order}.png")
            failures += [f"order {interp_order}, {f}" for f in check(rows)]

    if cpp.mpi_rank() == 0 and failures:
        raise AssertionError("dispersion check failed:\n  " + "\n  ".join(failures))


if __name__ == "__main__":
    main()
