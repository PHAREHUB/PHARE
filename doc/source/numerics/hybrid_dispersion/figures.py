"""
Figures of the hybrid dispersion test page, drawn at documentation build time.

Data come from the 20 calibration runs of tests/functional/dispersion/dispersion.py,
stored as small text tables next to this script (the matplotlib plot directive runs it
from this directory). Run settings below mirror RUNS in that test.
"""

import matplotlib.pyplot as plt
import numpy as np

RED, GREEN = "tab:red", "tab:green"
RUNS = {
    "low_k": dict(cells=128, dl=0.8, modes=(1, 2, 4, 8)),
    "high_k": dict(cells=64, dl=0.1, modes=(1, 2, 3, 4)),
}
TOLERANCE = 1e-2
AMPLITUDE = 1e-2


def omega_aic(k, di=1.0):
    return k * (np.sqrt(1 + (k * di) ** 2 / 4) - k * di / 2)


def omega_whistler(k, di=1.0):
    return k * (np.sqrt(1 + (k * di) ** 2 / 4) + k * di / 2)


def k_grid(k, dl):
    return 2 / dl * np.sin(k * dl / 2)


def read_table(filename):
    """Columns named by the first non-comment line (genfromtxt would take the comment)."""
    with open(filename) as f:
        lines = [line for line in f if not line.startswith("#")]
    return np.genfromtxt(lines, delimiter=",", names=True, dtype=None, encoding="utf-8")


def calibration():
    return read_table("calibration.csv")


def select(table, order, run, branch):
    rows = table[
        (table["order"] == order) & (table["run"] == run) & (table["branch"] == branch)
    ]
    return rows[np.argsort(rows["k"])]


def dispersion(order=1):
    """Measured w against theory (run 1), and error against w_scheme."""
    table = calibration()
    branches = {
        "whistler": ("Whistler (R)", omega_whistler, RED),
        "aic": ("Alfvén/ion-cyclotron (L)", omega_aic, GREEN),
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
        for branch, (_, _, color) in branches.items():
            rows = select(table, order, run, branch)
            style = dict(
                marker=marker,
                s=40,
                zorder=3,
                facecolor=color,
                edgecolor="black",
                lw=0.6,
            )
            ax.scatter(rows["k"], rows["w_run1"], **style)
            ax_err.scatter(
                rows["k"],
                100 * (rows["w_run1"] - rows["w_scheme"]) / rows["w_scheme"],
                **style,
            )
            for r in rows:
                ax.annotate(
                    f"m={r['m']}",
                    (r["k"], r["w_run1"]),
                    textcoords="offset points",
                    xytext=(5, 7 if branch == "whistler" else -12),
                    fontsize=7,
                )
        ax.scatter(
            [],
            [],
            marker=marker,
            facecolor="white",
            edgecolor="black",
            label=f"PHARE, {RUNS[run]['cells']} cells, $\\Delta x$ = {RUNS[run]['dl']} $d_i$",
        )

    ax.set_ylabel(r"$\omega\,/\,\Omega_{ci}$")
    ax.set_title(
        f"PHARE hybrid PIC — parallel dispersion relation, interpolation order {order}\n"
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


def loaded_modes():
    """(m, q, w) of the ion-scale run at order 1, as loaded by the test."""
    run = RUNS["high_k"]
    dl, k0 = run["dl"], 2 * np.pi / (run["cells"] * run["dl"])
    di = 1 - 2 * dl**2 / 12
    aic = [(m, m * k0, omega_aic(k_grid(m * k0, dl), di)) for m in run["modes"]]
    whistler = [
        (m, -m * k0, -omega_whistler(k_grid(m * k0, dl), di)) for m in run["modes"]
    ]
    return aic, whistler


def method():
    """Loaded field, its spectrum, and the measured phase of three Fourier bins."""
    run = RUNS["high_k"]
    aic, whistler = loaded_modes()
    x = (np.arange(run["cells"]) + 0.5) * run["dl"]
    phases = np.random.default_rng(0).uniform(
        0, 2 * np.pi, 2 * len(aic)
    )  # as in the test
    b = sum(
        AMPLITUDE * np.exp(1j * (q * x + phi))
        for (_, q, _), phi in zip(aic + whistler, phases)
    )

    fig, axs = plt.subplots(1, 3, figsize=(13, 3.6))
    ax = axs[0]
    ax.plot(x, b.real, color="k", lw=1.2, label=r"$B_y$")
    ax.plot(x, b.imag, color="gray", lw=1.2, ls="--", label=r"$B_z$")
    ax.set_xlabel(r"$x\ [d_i]$")
    ax.legend(
        loc="lower center", fontsize=9, ncol=2, bbox_to_anchor=(0.5, 1.0), frameon=False
    )
    ax.set_title("(a) loaded field: 8 circular waves at $t=0$", pad=22)

    ax = axs[1]
    qs = 2 * np.pi * np.fft.fftfreq(x.size, run["dl"])
    spec = np.abs(np.fft.fft(b)) / x.size
    order = np.argsort(qs)
    ax.stem(qs[order], spec[order], linefmt="0.6", markerfmt=" ", basefmt=" ")
    for modes, color in ((aic, GREEN), (whistler, RED)):
        for _, q, _ in modes:
            ax.plot(q, spec[np.argmin(abs(qs - q))], "o", color=color, mec="k")
    ax.text(
        0.97, 0.5, "AIC (L)\nat $+q$", color=GREEN, ha="right", transform=ax.transAxes
    )
    ax.text(
        0.03, 0.5, "whistler (R)\nat $-q$", color=RED, ha="left", transform=ax.transAxes
    )
    ax.set_ylim(0, 0.012)
    ax.set_xlim(-5, 5)
    ax.set_xlabel(r"$q\ [d_i^{-1}]$")
    ax.set_ylabel(r"$|\hat b(q)|$")
    ax.set_title(r"(b) Fourier transform of $b=B_y+iB_z$")

    ax = axs[2]
    measured = read_table("phases.csv")
    t = measured["t"]
    for column, (m, _, w), color in (
        ("aic_m1", aic[0], GREEN),
        ("aic_m4", aic[3], GREEN),
        ("whistler_m1", whistler[0], RED),
    ):
        ph = measured[column]
        ax.plot(t, ph, color=color, lw=2.2, alpha=0.5)
        ax.plot(t, -w * t, "k:", lw=1)
        mid = len(t) // 2 if w < 0 or m == 1 else 3 * len(t) // 4
        ax.annotate(
            f"{'AIC' if w > 0 else 'whistler'} m={m}, $\\omega$={abs(w):.2f}",
            (t[mid], ph[mid]),
            fontsize=8,
            xytext=(0, 9 if w < 0 or m == 1 else -14),
            textcoords="offset points",
            ha="right" if w < 0 else "center",
        )
    ax.set_xlabel(r"$t\ [\Omega_{ci}^{-1}]$")
    ax.set_ylabel(r"phase of $\hat b(q,t)$ [rad]")
    ax.set_title(r"(c) phase $= -\omega t$: slope gives $\omega$")
    ax.plot([], [], "k:", label="expected slope")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()


def corrections():
    """Grid derivative factor, and the measured Hall splitting deficit against the shape formula."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    ax = axs[0]
    kdx = np.linspace(1e-3, np.pi, 200)
    ax.plot(kdx, np.sin(kdx / 2) / (kdx / 2), "k")
    ax.axvspan(
        0,
        2 * np.pi / 16,
        color="tab:blue",
        alpha=0.15,
        label=r"test modes ($\geq 16$ cells / $\lambda$)",
    )
    ax.set_xlabel(r"$k\Delta x$")
    ax.set_ylabel(r"$k_{grid}/k$")
    ax.set_title(
        r"(a) grid derivative: $k_{grid}=\frac{2}{\Delta x}\sin\frac{k\Delta x}{2}$"
    )
    ax.legend(loc="lower left")
    ax.grid(alpha=0.3)

    ax = axs[1]
    split = read_table("hall_splitting.csv")
    dxs = np.linspace(0, 0.9, 100)
    colors = {1: "tab:blue", 2: "tab:orange", 3: "tab:purple"}
    for p, color in colors.items():
        ax.plot(
            dxs,
            -100 * (p + 1) * dxs**2 / 12,
            color=color,
            label=rf"$-(p+1)\Delta x^2/12$, p={p}",
        )
        for r in split[split["order"] == p]:
            yerr = 0 if np.isnan(r["err_std"]) else 100 * r["err_std"]
            ax.errorbar(
                r["dx"], 100 * r["err_mean"], yerr=yerr, fmt="o", color=color, mec="k"
            )
    ax.plot([], [], "o", color="w", mec="k", label="PHARE (measured)")
    ax.set_xlabel(r"$\Delta x\ [d_i]$")
    ax.set_ylabel(r"error on $\omega_{wh}-\omega_{aic}$ [%]")
    ax.set_title(
        r"(b) particle shape: Hall splitting $\propto d_i^{eff}=1-\frac{(p+1)\Delta x^2}{12}$"
    )
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()


def calibration_spread():
    """Error of every point, mean +- 4 sigma over the calibration runs."""
    table = calibration()
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for p, ax in zip((1, 2, 3), axs):
        for run, marker in (("low_k", "o"), ("high_k", "s")):
            for branch, color in (("aic", GREEN), ("whistler", RED)):
                rows = select(table, p, run, branch)
                ax.errorbar(
                    rows["k"],
                    100 * rows["err_mean"],
                    yerr=400 * rows["err_std"],
                    fmt=marker,
                    color=color,
                    mec="k",
                    ms=6,
                    capsize=2,
                    lw=1,
                )
        ax.axhspan(-100 * TOLERANCE, 100 * TOLERANCE, color="tab:blue", alpha=0.08)
        for sign in (1, -1):
            ax.axhline(sign * 100 * TOLERANCE, color="tab:blue", ls="--", lw=1)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xscale("log")
        ax.set_xlabel(r"$k\,d_i$")
        ax.set_title(f"interpolation order {p}")
        ax.grid(alpha=0.3, which="both")
    axs[0].set_ylabel(r"$(\omega_{PHARE}-\omega_{scheme})/\omega_{scheme}$ [%]")
    axs[0].errorbar(
        [],
        [],
        yerr=[],
        fmt="o",
        color="w",
        mec="k",
        label=r"mean $\pm 4\sigma$, 20 runs",
    )
    axs[0].plot([], [], "--", color="tab:blue", label="tolerance")
    axs[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
