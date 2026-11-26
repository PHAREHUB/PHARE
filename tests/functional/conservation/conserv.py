#!/usr/bin/env python3
"""
How well the hybrid scheme conserves energy, depending on the regime.
Nightly. The fast check at one good regime is tests/simulator/test_conservation.py.

Total energy E = K + W_B + U_e (ions, magnetic field, isothermal electrons,
see tests/simulator/invariants.py) is conserved by the continuous equations,
not by the discrete scheme: particle noise exchanges energy with the fields
and heats the plasma. How much depends on the regime. Two parts:

1. Regime map. Uniform plasma at rest, n = 1, B = x, 1D, interpolation order 1,
   Vth from 0.03 to 2 at Te = 0.1 and 0.4. The change of E and of its three
   parts, with their sign, relative to the initial ion kinetic energy K0:
   - cold ions (Ti << Te): E drifts up by several % to tens of % of K0, more
     at larger Te for the same Vth. Only the drift of E measures the scheme:
     K alone also gains the electron energy of the initial density noise
     (U_e starts about Te / (3 ppc) per cell above its smooth value), a
     relaxation, not an error. Whether Ti/Te alone sets the drift is not shown
     here (no two runs share Ti/Te). This is the known weakness of
     quasineutral hybrid codes with cold ions, a property of the scheme.
   - Ti >= Te: E is conserved to about 0.1% of K0 over the run.
   - hot ions (high beta): the ions lose energy to magnetic fluctuations that
     grow from the smooth initial field; E moves less than K.
   Checked: the bound at Ti >= Te, a ceiling on the cold-ion drift (a
   regression guard: it must not get worse), the hot-ion exchange.

2. Numerical parameters, 2D, Ti = Te, where the drift is measurable.
   Energy drift = slope of E(t), per seed; the slopes of K and W_B are shown
   too: the ions lose energy and the magnetic field gains more (noise
   fluctuations growing from the smooth initial B), so E drifts up. Scans of
   particles per cell, time step, cell size and interpolation order around one
   reference. The cell size scan keeps the box and particles per cell, and
   keeps dt / dx^2 fixed so the whistler time step stability stays the same
   fraction of its limit. Checked: drift ~ 1/ppc, independent of dt, lower at
   higher order. dt independence shows the drift is noise, not a time
   integration error. The cell size scan is shown, not checked.

Figures: conservation.png (regime map), conservation_numerics.png (scans).
The previous version of this test measured |K(t=50) - K(t=3..4)| / K only:
the ion part alone, unsigned, at one time, hiding the exchange with B.

Each run is a separate process (one simulation per process), run in parallel.
"""

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

nbr_measures = 40


# ------------------------------------------------------------------ one run


def run_one(ndim, interp, cells, dl, patch, dt, steps, ppc, vth, Te, seed):
    """uniform plasma at rest, returns time series of K, W_B, U_e"""
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator
    from tests.simulator.invariants import invariants

    ph.NO_GUI()
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

    def constant(value):
        return lambda *xyz: value + 0.0 * xyz[0]

    ph.MaxwellianFluidModel(
        bx=constant(1.0),
        by=constant(0.0),
        bz=constant(0.0),
        protons={
            "charge": 1,
            "mass": 1.0,
            "density": constant(1.0),
            **{f"vbulk{c}": constant(0.0) for c in "xyz"},
            **{f"vth{c}": constant(vth) for c in "xyz"},
            "nbr_part_per_cell": ppc,
            "init": {"seed": seed},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=Te)

    masses = {"protons": 1.0}
    # one log per process would collide: all runs are rank 0 in the same directory
    simulator = Simulator(sim, log_to_file=False).initialize()
    if sim.dry_run:  # pull-request CI only checks that the setup is valid
        return None
    every = max(1, steps // nbr_measures)
    series = {"t": [], "K": [], "W_B": [], "U_e": [], "N": [], "n_min": []}

    def measure(step):
        inv = invariants(simulator, masses, Te, ndim)
        series["t"].append(step * dt)
        for key in ("K", "W_B", "U_e"):
            series[key].append(float(inv[key]))
        series["N"].append(float(inv["N"]["protons"]))
        series["n_min"].append(float(inv["n_min"]))

    measure(0)
    for step in range(1, steps + 1):
        simulator.advance()
        if step % every == 0:
            measure(step)
    simulator.reset()
    return series


def run_all(jobs):
    """runs each job in its own process, in parallel"""
    nbr_workers = max(1, min(len(jobs), (os.cpu_count() or 2) // 2, 12))

    def launch(job):
        out = subprocess.run(
            [sys.executable, "-u", __file__, "--run", json.dumps(job)],
            capture_output=True,
            text=True,
            check=False,  # a failed run is reported below with its output
        )
        for line in out.stdout.splitlines():
            if line.startswith("RESULT "):
                return json.loads(line[len("RESULT ") :])
        raise RuntimeError(
            f"run {job} failed:\n{out.stdout[-2000:]}\n{out.stderr[-2000:]}"
        )

    with ThreadPoolExecutor(nbr_workers) as pool:
        return list(pool.map(launch, jobs))


def energies(series):
    """changes of E, K, W_B, U_e relative to K0, as arrays over time"""
    K0 = series["K"][0]
    d = {
        key: (np.asarray(series[key]) - series[key][0]) / K0
        for key in ("K", "W_B", "U_e")
    }
    d["E"] = d["K"] + d["W_B"] + d["U_e"]
    d["t"] = np.asarray(series["t"])
    d["dN"] = np.abs(np.asarray(series["N"]) / series["N"][0] - 1).max()
    d["n_min"] = min(series["n_min"])
    d["finite"] = all(np.isfinite(d[k]).all() for k in ("E", "K", "W_B", "U_e"))
    return d


# ------------------------------------------------------------- part 1: map

map_Te = (0.1, 0.4)
map_vth = (0.03, 0.1, 0.3, 0.6, 1.0, 2.0)
map_run = dict(
    ndim=1, interp=1, cells=100, dl=0.2, patch=20, dt=0.002, steps=10000, ppc=100
)

# bounds from 3 seeds, Release build:
# - Ti >= Te: worst max |dE| / K0 = 0.0015 (Vth = 2), bound 0.005
# - cold ions, final dE / K0 at Vth = 0.03: 0.066-0.114 at Te = 0.1,
#   0.099-0.266 at Te = 0.4, ceilings about 2.5 times the worst seed
# - hot ions, Vth = 2, mean over the second half: dK = -0.0034 to -0.0044,
#   dW_B = +0.0025 to +0.0039 of K0, same sign for all seeds and both Te
good_regime_bound = 0.005
# Ti < Te, final dE / K0 per (Te, vth), worst of 3 seeds: 0.114, 0.006 (Te=0.1);
# 0.266, 0.042, 0.0058 (Te=0.4). Ceilings about 2.5 times that: the cold-ion
# heating must not get worse
cold_heating_ceiling = {
    (0.1, 0.03): 0.3,
    (0.1, 0.1): 0.015,
    (0.4, 0.03): 0.7,
    (0.4, 0.1): 0.1,
    (0.4, 0.3): 0.015,
}


def regime_map():
    jobs = [dict(map_run, vth=v, Te=Te, seed=1) for Te in map_Te for v in map_vth]
    return {(j["Te"], j["vth"]): energies(r) for j, r in zip(jobs, run_all(jobs))}


def check_map(results):
    failures = []
    for (Te, vth), d in results.items():
        Ti_Te = vth**2 / Te
        if not d["finite"] or d["n_min"] <= 0:
            failures.append(
                f"Te={Te} vth={vth}: non finite energy or n_min {d['n_min']}"
            )
            continue
        if d["dN"] > 1e-12:
            failures.append(f"Te={Te} vth={vth}: ion number changed by {d['dN']:.1e}")
        final = {k: d[k][-1] for k in ("E", "K", "W_B", "U_e")}
        print(
            f"Te {Te} vth {vth} Ti/Te {Ti_Te:.3g}: dE {final['E']:+.2e} = "
            f"dK {final['K']:+.2e} + dW_B {final['W_B']:+.2e} + dU_e {final['U_e']:+.2e} (of K0)"
        )
        # every point is asserted: a cold-ion ceiling (Ti < Te) or the good bound
        # (Ti/Te >= 0.9, which is 0.8999.. in floating point for two points)
        if (Te, vth) not in cold_heating_ceiling:
            assert Ti_Te > 0.85, (Te, vth)
            worst = np.abs(d["E"]).max()
            if worst > good_regime_bound:
                failures.append(
                    f"Ti/Te={Ti_Te:.3g}: max|dE|/K0 {worst:.2e} > {good_regime_bound}"
                )
    for (Te, vth), ceiling in cold_heating_ceiling.items():
        heating = results[(Te, vth)]["E"][-1]
        if not heating < ceiling:
            failures.append(
                f"cold ions, Te={Te} vth={vth}: dE/K0 {heating:.3g} > {ceiling}"
                " (heats more than before)"
            )
    for Te in map_Te:
        # hot ions: the ions lose energy and W_B gains it, over the second half
        d = results[(Te, max(map_vth))]
        late = d["t"] >= d["t"][-1] / 2
        dK, dW = d["K"][late].mean(), d["W_B"][late].mean()
        if not (dK < 0 < dW):
            failures.append(
                f"hot ions, Te={Te}: expected dK < 0 < dW_B, got {dK:.2e}, {dW:.2e}"
            )
    return failures


def plot_map(results, filename="conservation.png"):
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    styles = {
        "E": ("k", "o-", "total E"),
        "K": ("C0", "s--", "ions K"),
        "W_B": ("C1", "^--", "magnetic W_B"),
        "U_e": ("C2", "v--", "electrons U_e"),
    }
    for ax, Te in zip(axes, map_Te):
        Ti_Te = np.array(map_vth) ** 2 / Te
        for key, (color, style, label) in styles.items():
            ax.plot(
                Ti_Te,
                [100 * results[(Te, v)][key][-1] for v in map_vth],
                style,
                color=color,
                label=label,
            )
        ax.axvline(1.0, color="grey", lw=0.8, ls=":")
        ax.axhline(0.0, color="grey", lw=0.8)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=0.1)
        ax.set_xlabel("Ti / Te")
        ax.set_title(
            f"Te = {Te}, 1D, order 1, T = {map_run['steps'] * map_run['dt']:g}"
        )
    axes[0].set_ylabel("change at the end, % of initial ion energy")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(filename, dpi=120)


# ----------------------------------------------- part 2: numerical parameters

# 2D, order 1, Ti = Te, box 12 x 12, run to T = 2. Each scan changes one
# parameter from this reference. dx changes the number of cells, not the box
scaling_run = dict(
    ndim=2,
    interp=1,
    cells=24,
    dl=0.5,
    patch=12,
    dt=0.01,
    steps=200,
    ppc=50,
    vth=0.3,
    Te=0.09,
)
box = scaling_run["cells"] * scaling_run["dl"]
final_time = scaling_run["steps"] * scaling_run["dt"]
scaling_seeds = (1, 2)
scans = {
    "particles per cell": ("ppc", (25, 50, 100, 200)),
    "time step dt": ("dt", (0.02, 0.01, 0.005, 0.0025)),
    "cell size dx": ("dl", (1.0, 0.5, 0.25)),
    "interpolation order": ("interp", (1, 2, 3)),
}
# measured over 4 seeds, energy drift in K0 per unit time: reference 2.5e-3 to
# 4.5e-3; 4x particles 3.6x less on average; dt / 2 equal to 3 digits for every
# seed; order 2 2.5x less, order 3 4.8x less


def variant(param, value):
    job = dict(scaling_run, **{param: value})
    if param == "dt":
        job["steps"] = round(final_time / value)
    if param == "dl":
        job["cells"] = round(box / value)
        # same whistler stability margin: dt ~ dx^2
        job["dt"] = scaling_run["dt"] * (value / scaling_run["dl"]) ** 2
        job["steps"] = round(final_time / job["dt"])
    return job


def drift_rates(d):
    """slopes over time of E, K and W_B, in K0 per unit time"""
    return {key: np.polyfit(d["t"], d[key], 1)[0] for key in ("E", "K", "W_B")}


def scaling():
    """drift rates (E, K, W_B) per seed for every value of every scan"""
    configs = {}
    for param, values in scans.values():
        for value in values:
            for seed in scaling_seeds:
                job = dict(variant(param, value), seed=seed)
                configs[json.dumps(job, sort_keys=True)] = job
    jobs = list(configs.values())
    rates = {}
    for job, r in zip(jobs, run_all(jobs)):
        d = energies(r)
        if not d["finite"] or d["n_min"] <= 0 or d["dN"] > 1e-12:
            raise RuntimeError(f"scaling run {job} broken: {d['n_min']=}, {d['dN']=}")
        rates[json.dumps(job, sort_keys=True)] = drift_rates(d)
    out = {}
    for name, (param, values) in scans.items():
        out[name] = {
            value: [
                rates[json.dumps(dict(variant(param, value), seed=s), sort_keys=True)]
                for s in scaling_seeds
            ]
            for value in values
        }
    return out


def seed_mean(per_seed, key):
    return float(np.mean([r[key] for r in per_seed]))


def check_scaling(scan_rates):
    for name, per_value in scan_rates.items():
        print(
            f"2D drift vs {name} (K0 per unit time, E = K + W_B + U_e): "
            + ", ".join(
                f"{v}: E {seed_mean(r, 'E'):+.2e} K {seed_mean(r, 'K'):+.2e}"
                f" W_B {seed_mean(r, 'W_B'):+.2e}"
                for v, r in per_value.items()
            )
        )
    mean = {
        name: {v: seed_mean(r, "E") for v, r in per_value.items()}
        for name, per_value in scan_rates.items()
    }
    failures = []
    ref = mean["particles per cell"][scaling_run["ppc"]]
    if not ref > 0:
        return [f"no measurable energy drift in the reference run ({ref:.2e})"]

    ppc = mean["particles per cell"]
    slope = np.polyfit(np.log(list(ppc)), np.log(list(ppc.values())), 1)[0]
    print(f"energy drift ~ ppc^{slope:.2f}")
    dts = mean["time step dt"]
    spread = max(dts.values()) / min(dts.values()) - 1
    orders = {order: abs(drift) for order, drift in mean["interpolation order"].items()}
    checks = {
        f"drift ~ 1/ppc (slope {slope:.2f}, expected -1 +- 0.4)": abs(slope + 1) < 0.4,
        "energy drifts up for every time step": min(dts.values()) > 0,
        f"time step does not matter (spread {spread:.1%} < 10%)": min(dts.values()) > 0
        and spread < 0.1,
        "order 2 drifts at least 1.5x less than order 1": orders[2] < orders[1] / 1.5,
        "order 3 drifts at least 2x less than order 1": orders[3] < orders[1] / 2,
    }
    failures += [name for name, ok in checks.items() if not ok]
    return failures


def plot_scaling(scan_rates, filename="conservation_numerics.png"):
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, NullLocator

    ink, muted = "#222222", "#888888"
    parts = {  # key: color, marker, label
        "E": ("#222222", "o", "total E"),
        "W_B": ("#d9822b", "^", "magnetic W_B"),
        "K": ("#2a6fb0", "s", "ions K"),
    }
    titles = {
        "ppc": "particles per cell",
        "dt": "time step dt",
        "dl": "cell size dx  (dt ~ dx\u00b2)",
        "interp": "interpolation order",
    }

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.2), sharey=True)
    for ax, (name, per_value) in zip(axes, scan_rates.items()):
        param = scans[name][0]
        values = np.array(list(per_value), dtype=float)
        ref = list(per_value).index(scaling_run[param])
        ax.axhline(0.0, color=muted, lw=0.8)
        for key, (color, marker, _) in parts.items():
            rates = 100 * np.array([[r[key] for r in per_value[v]] for v in per_value])
            mean = rates.mean(axis=1)
            for k in range(rates.shape[1]):
                ax.plot(
                    values,
                    rates[:, k],
                    marker,
                    color=color,
                    alpha=0.3,
                    ms=4,
                    mec="none",
                )
            ax.plot(values, mean, "-", marker=marker, color=color, lw=1.8, ms=6)
            if key == "E":
                ax.plot(values[ref], mean[ref], "o", ms=13, mfc="none", mec=ink, mew=1)
                if param == "ppc":
                    ax.plot(
                        values,
                        mean[ref] * values[ref] / values,
                        "--",
                        color=muted,
                        lw=1,
                    )
        if param != "interp":
            ax.set_xscale("log")
        ax.xaxis.set_major_locator(FixedLocator(values))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xticklabels([f"{v:g}" for v in values])
        ax.set_title(titles[param], fontsize=11, color=ink)
        ax.margins(x=0.12)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.grid(axis="y", color="#dddddd", lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(colors=ink, labelsize=9)
    axes[0].set_ylabel("drift rate\n(% of initial ion energy per unit time)", color=ink)

    handles = [
        Line2D([], [], color=c, marker=m, lw=1.8, ms=6, label=label)
        for c, m, label in parts.values()
    ] + [
        Line2D(
            [],
            [],
            ls="none",
            marker="o",
            color=muted,
            alpha=0.5,
            ms=4,
            mec="none",
            label="one seed",
        ),
        Line2D(
            [],
            [],
            ls="none",
            marker="o",
            ms=13,
            mfc="none",
            mec=ink,
            label="reference (order 1, 50 ppc, dt 0.01, dx 0.5)",
        ),
        Line2D([], [], ls="--", color=muted, lw=1, label="1 / ppc"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=9)
    fig.suptitle(
        "Energy drift vs numerical parameters  "
        f"(2D uniform plasma, Ti = Te, T = {final_time:g}; E = K + W_B + U_e)",
        fontsize=12,
        color=ink,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(filename, dpi=130)


# --------------------------------------------------------------------- main


def dry_run():
    """builds every configuration without running it (PHARE_DRY_RUN=1)"""
    jobs = [dict(map_run, vth=v, Te=Te, seed=1) for Te in map_Te for v in map_vth]
    jobs += [
        dict(variant(p, v), seed=1) for p, values in scans.values() for v in values
    ]
    for job in jobs:
        run_one(**job)
    print(f"dry run: {len(jobs)} configurations built")


def main():
    if os.environ.get("PHARE_DRY_RUN", "0") == "1":
        return dry_run()
    start = time.time()
    results = regime_map()
    failures = check_map(results)
    plot_map(results)
    scan_rates = scaling()
    plot_scaling(scan_rates)
    failures += check_scaling(scan_rates)
    print(f"total {time.time() - start:.0f} s")
    if failures:
        print("FAILED:\n  " + "\n  ".join(failures))
        sys.exit(1)
    print("OK")


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--run":
        print("RESULT " + json.dumps(run_one(**json.loads(sys.argv[2]))))
    else:
        main()
