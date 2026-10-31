import json
import os

import matplotlib as mpl
import numpy as np

import pyphare.pharein as ph
from pyphare.pharein import global_vars
from pyphare.pharein.diagnostics import FluidDiagnostics
from pyphare.pharesee.hierarchy import fromfunc, hierarchy_from
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from tests.diagnostic import all_timestamps

mpl.use("Agg")


ncell = 100
dl = 0.2
L = ncell * dl
ts = 0.01
masses = (2, 3)
charges = (1, 2)


def densityMain_1d(x):
    return 1.0


def densityBeam_1d(x):
    u = x / L - 0.5
    return np.exp(-(u**2))


def bx_1d(x):
    return 1.0


def by_1d(x):
    return 0.0


def bz_1d(x):
    return 0.0


def v0_1d(x):
    return 0.0


def vth_1d(x):
    return np.sqrt(1.0)


def config_1d():
    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=60,
        time_step=ts,
        time_step_nbr=1,
        cells=ncell,
        dl=dl,
        diag_options={
            "format": "phareh5",
            "options": {"dir": "nCheck_1d", "mode": "overwrite"},
        },
    )

    v_pop = {
        "vbulkx": v0_1d,
        "vbulky": v0_1d,
        "vbulkz": v0_1d,
        "vthx": vth_1d,
        "vthy": vth_1d,
        "vthz": vth_1d,
    }

    ph.MaxwellianFluidModel(
        bx=bx_1d,
        by=by_1d,
        bz=bz_1d,
        main={
            "mass": masses[0],
            "charge": charges[0],
            "density": densityMain_1d,
            "nbr_part_per_cell": 1000,
            **v_pop,
        },
        beam={
            "mass": masses[1],
            "charge": charges[1],
            "density": densityBeam_1d,
            "nbr_part_per_cell": 1000,
            **v_pop,
        },
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    timestamps = all_timestamps(global_vars.sim)

    for quantity in ["charge_density", "mass_density"]:
        FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    poplist = ["main", "beam"]
    for pop in poplist:
        for quantity in ["density", "charge_density"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                population_name=pop,
            )

    return sim


def densityMain_2d(x, y):
    assert len(x) == len(y)
    return 1.0 * np.ones_like(x)


def densityBeam_2d(x, y):
    assert len(x) == len(y)
    u = x / L - 0.5
    v = y / L - 0.5
    return np.exp(-(u**2) - v**2)


def bx_2d(x, y):
    return 1.0


def by_2d(x, y):
    return 0.0


def bz_2d(x, y):
    return 0.0


def v0_2d(x, y):
    return 0.0


def vth_2d(x, y):
    return np.sqrt(1.0)


def config_2d():
    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=60,
        time_step=ts,
        time_step_nbr=1,
        cells=(ncell, ncell),
        dl=(dl, dl),
        diag_options={
            "format": "phareh5",
            "options": {"dir": "nCheck_2d", "mode": "overwrite"},
        },
    )

    v_pop = {
        "vbulkx": v0_2d,
        "vbulky": v0_2d,
        "vbulkz": v0_2d,
        "vthx": vth_2d,
        "vthy": vth_2d,
        "vthz": vth_2d,
    }

    ph.MaxwellianFluidModel(
        bx=bx_2d,
        by=by_2d,
        bz=bz_2d,
        main={
            "mass": masses[0],
            "charge": charges[0],
            "density": densityMain_2d,
            "nbr_part_per_cell": 1000,
            **v_pop,
        },
        beam={
            "mass": masses[1],
            "charge": charges[1],
            "density": densityBeam_2d,
            "nbr_part_per_cell": 1000,
            **v_pop,
        },
    )

    ph.ElectronModel(closure="isothermal", Te=0.0)

    timestamps = all_timestamps(global_vars.sim)

    for quantity in ["charge_density", "mass_density"]:
        FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

    poplist = ["main", "beam"]
    for pop in poplist:
        for quantity in ["density", "charge_density"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                population_name=pop,
            )

    return sim


# Each check compares the density PHARE computes from the particles with the density
# given as input. Particles are drawn at random (no seed), so the spread of the difference
# is particle noise and changes from run to run. The largest spread accepted per patch,
# for each (dimension, quantity), is read from THRESHOLDS_FILE. Those values were measured
# over many runs by res/calibration/density_check/calibrate.py: rerun it to update them
# when this test's setup changes (cells, particles per cell, masses, charges, densities).
THRESHOLDS_FILE = "density_check_thresholds.json"


def load_thresholds():
    """{"1d": {"mass": t, "charge": t}, "2d": {...}} from THRESHOLDS_FILE.

    CMake copies the file next to this test in the build directory; when the test
    runs from the source tree, it is read from res/calibration/density_check/.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(here, THRESHOLDS_FILE),
        os.path.join(
            here,
            *[os.pardir] * 3,
            "res",
            "calibration",
            "density_check",
            THRESHOLDS_FILE,
        ),
    ]
    for path in candidates:
        if os.path.exists(path):
            with open(path) as f:
                return {
                    dim: {qty: v["threshold"] for qty, v in qties.items()}
                    for dim, qties in json.load(f)["thresholds"].items()
                }
    raise FileNotFoundError(f"{THRESHOLDS_FILE} not found in {candidates}")


def computed_and_expected(ndim, quantity, time=0.0):
    """(computed, expected) hierarchies of the total ion mass or charge density."""
    r = Run(os.path.join(os.curdir, f"nCheck_{ndim}d"))
    densities = {
        1: (densityMain_1d, densityBeam_1d),
        2: (densityMain_2d, densityBeam_2d),
    }[ndim]
    if quantity == "mass":
        h = r.GetMassDensity(time)
        func = {
            1: fromfunc.ions_mass_density_func1d,
            2: fromfunc.ions_mass_density_func2d,
        }[ndim]
        return h, hierarchy_from(hier=h, from_func=func, masses=masses, densities=densities)
    h = r.GetNi(time)
    func = {
        1: fromfunc.ions_charge_density_func1d,
        2: fromfunc.ions_charge_density_func2d,
    }[ndim]
    return h, hierarchy_from(hier=h, from_func=func, charges=charges, densities=densities)


def patch_sigmas(h, H, time=0.0):
    """Standard deviation of (computed - expected) on each patch."""
    sigmas = []
    for lvl_h, lvl_H in zip(h.levels(time).values(), H.levels(time).values()):
        for patch_h, patch_H in zip(lvl_h.patches, lvl_H.patches):
            dset_h = patch_h.patch_datas["value"][patch_h.box]
            dset_H = patch_H.patch_datas["value"][patch_H.box]
            sigmas.append(float(np.std(dset_h - dset_H)))
    return sigmas


def main():
    import matplotlib.pyplot as plt

    Simulator(config_1d()).run().reset()
    ph.global_vars.sim = None
    Simulator(config_2d()).run().reset()

    thresholds = load_thresholds()
    hiers = {}
    for ndim in (1, 2):
        for quantity in ("mass", "charge"):
            h, H = computed_and_expected(ndim, quantity)
            hiers[(ndim, quantity)] = (h, H)
            threshold = thresholds[f"{ndim}d"][quantity]
            for sigma in patch_sigmas(h, H):
                print(
                    f"dim = {ndim}, {quantity}, sigma(user v - actual v) = {sigma}, threshold = {threshold}"
                )
                assert sigma < threshold

    _fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(6, 8))
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for ax, quantity in ((ax1, "mass"), (ax2, "charge")):
        h, H = hiers[(1, quantity)]
        h.plot(ax=ax, ls="-", lw=2.0, color=cycle[0])
        H.plot(ax=ax, ls="-", lw=2.0, color=cycle[1])
        ax.set_title(f"{quantity} density : 1d")

    cmap = mpl.colormaps["viridis"]
    h1, H1 = hiers[(2, "mass")]
    h2, H2 = hiers[(2, "charge")]
    h1.plot(ax=ax3, vmin=3.75, vmax=5, cmap=cmap, title="computed mass density : 2d")
    H1.plot(ax=ax4, vmin=3.75, vmax=5, cmap=cmap, title="expected mass density : 2d")
    h2.plot(ax=ax5, vmin=2.2, vmax=3, cmap=cmap, title="computed charge density : 2d")
    H2.plot(ax=ax6, vmin=2.2, vmax=3, cmap=cmap, title="expected charge density : 2d")

    plt.tight_layout()
    plt.savefig("nCheck.pdf", dpi=300)


if __name__ == "__main__":
    main()
