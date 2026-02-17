import os
import numpy as np
import matplotlib as mpl

import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.pharein import global_vars
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharein.diagnostics import FluidDiagnostics
from pyphare.pharesee.hierarchy import fromfunc

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
        boundary_types="periodic",
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
        boundary_types=("periodic", "periodic"),
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


def main():
    import matplotlib.pyplot as plt

    Simulator(config_1d()).run().reset()
    ph.global_vars.sim = None
    Simulator(config_2d()).run().reset()

    def assert_close_enough(h, H):
        for lvl_h, lvl_H in zip(h.levels(time).values(), H.levels(time).values()):
            for patch_h, patch_H in zip(lvl_h.patches, lvl_H.patches):
                pd_h = patch_h.patch_datas["value"]
                pd_H = patch_H.patch_datas["value"]

                dset_h = pd_h[patch_h.box]
                dset_H = pd_H[patch_H.box]

                std = np.std(dset_h - dset_H)
                print("dim = {}, sigma(user v - actual v) = {}".format(pd_H.ndim, std))
                assert std < 0.062  # empirical value obtained from print just above

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(6, 8))

    # 1d stuffs
    run_path = os.path.join(os.curdir, "nCheck_1d")
    time = 0.0
    r = Run(run_path)

    h1 = r.GetMassDensity(time)
    h2 = r.GetNi(time)

    H1 = hierarchy_from(
        hier=h1,
        func=fromfunc.ions_mass_density_func1d,
        masses=masses,
        densities=(densityMain_1d, densityBeam_1d),
    )
    H2 = hierarchy_from(
        hier=h2,
        func=fromfunc.ions_charge_density_func1d,
        charges=charges,
        densities=(densityMain_1d, densityBeam_1d),
    )

    assert_close_enough(h1, H1)
    assert_close_enough(h2, H2)

    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    h1.plot(ax=ax1, ls="-", lw=2.0, color=cycle[0])
    H1.plot(ax=ax1, ls="-", lw=2.0, color=cycle[1])

    h2.plot(ax=ax2, ls="-", lw=2.0, color=cycle[0])
    H2.plot(ax=ax2, ls="-", lw=2.0, color=cycle[1])

    ax1.set_title("mass density : 1d")
    ax2.set_title("charge density : 1d")

    # 2d stuffs
    run_path = os.path.join(os.curdir, "nCheck_2d")
    time = 0.0
    r = Run(run_path)

    h1 = r.GetMassDensity(time)
    h2 = r.GetNi(time)

    H1 = hierarchy_from(
        hier=h1,
        func=fromfunc.ions_mass_density_func2d,
        masses=masses,
        densities=(densityMain_2d, densityBeam_2d),
    )
    H2 = hierarchy_from(
        hier=h2,
        func=fromfunc.ions_charge_density_func2d,
        charges=charges,
        densities=(densityMain_2d, densityBeam_2d),
    )

    assert_close_enough(h1, H1)
    assert_close_enough(h2, H2)

    cmap = mpl.colormaps["viridis"]

    h1.plot(ax=ax3, vmin=3.75, vmax=5, cmap=cmap, title="computed mass density : 2d")
    H1.plot(ax=ax4, vmin=3.75, vmax=5, cmap=cmap, title="expected mass density : 2d")
    h2.plot(ax=ax5, vmin=2.2, vmax=3, cmap=cmap, title="computed charge density : 2d")
    H2.plot(ax=ax6, vmin=2.2, vmax=3, cmap=cmap, title="expected charge density : 2d")

    plt.tight_layout()
    plt.savefig("nCheck.pdf", dpi=300)


if __name__ == "__main__":
    main()
