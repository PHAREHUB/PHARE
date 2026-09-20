import sys
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt

from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare


atol = 1e-12

outputpath = Path("phare_outputs/compare_diags/")
outputs = str(outputpath)


def safe(s):
    return "".join(c if c.isalnum() else "_" for c in s)


def plot_patch_diff(qty, reason, ref_data, cmp_data):
    layout = ref_data.layout
    diff = cmp_data.dataset[:] - ref_data.dataset[:]

    absmax = np.max(np.abs(diff))
    vmin, vmax = -absmax, absmax
    if absmax == 0:
        vmin, vmax = -1e-16, 1e-16

    fig, ax = plt.subplots()
    im = ax.pcolormesh(
        layout.yeeCoordsFor(qty, "x", withGhosts=True),
        layout.yeeCoordsFor(qty, "y", withGhosts=True),
        diff.T,
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(f"{qty} diff\n{reason}\nmax|diff|={absmax:.3e}")
    plt.colorbar(im, ax=ax)
    fig.savefig(outputs + f"/{safe(qty)}_{safe(reason)}_diff.png", dpi=200)
    plt.close(fig)


def plot_patch_side_by_side(qty, reason, ref_data, cmp_data):
    layout = ref_data.layout
    x = layout.yeeCoordsFor(qty, "x", withGhosts=True)
    y = layout.yeeCoordsFor(qty, "y", withGhosts=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    for ax, data, label in (
        (axes[0], ref_data, "ref"),
        (axes[1], cmp_data, "cmp"),
    ):
        vals = data.dataset[:]
        im = ax.pcolormesh(x, y, vals.T, cmap="Spectral_r")
        ax.set_title(label)
        plt.colorbar(im, ax=ax)
    fig.suptitle(f"{qty}\n{reason}")
    fig.savefig(outputs + f"/{safe(qty)}_{safe(reason)}_values.png", dpi=200)
    plt.close(fig)


def check_and_plot(qty_label, h0, h1, plot=True):
    eqr = hierarchy_compare(h0, h1, atol=atol)
    print(f"\n\t{qty_label}", eqr)

    try:
        if plot:
            for reason, cmp_data, ref_data in eqr.failed:
                short_reason = reason.split("\n")[0]
                plot_patch_diff(cmp_data.field_name, short_reason, ref_data, cmp_data)
                plot_patch_side_by_side(
                    cmp_data.field_name, short_reason, ref_data, cmp_data
                )
    except Exception:
        ...
    return bool(eqr)


def main():
    run0 = Run(sys.argv[1])
    run1 = Run(sys.argv[2])

    times = run0.all_times()["B"]
    print("times", times)
    for time in times:
        print("\ntime :", time)

        # for pop in run0.all_pops():
        #     for kind in ["domain", "levelGhost"]:
        #         p0 = run0.GetParticles(time, pop, type=kind)
        #         p1 = run1.GetParticles(time, pop, type=kind)
        #         check_and_plot(f"particles {kind} {pop}", p0, p1, plot=False)

        check_and_plot(
            "B", run0.GetB(time, all_primal=False), run1.GetB(time, all_primal=False)
        )
        check_and_plot(
            "E", run0.GetE(time, all_primal=False), run1.GetE(time, all_primal=False)
        )

        for pop in run0.all_pops():
            check_and_plot(
                f"rho {pop}",
                run0.GetN(time, pop_name=pop),
                run1.GetN(time, pop_name=pop),
            )
            check_and_plot(
                f"flux {pop}",
                run0.GetFlux(time, pop_name=pop),
                run1.GetFlux(time, pop_name=pop),
            )

        check_and_plot("Ni", run0.GetNi(time), run1.GetNi(time))
        check_and_plot(
            "V", run0.GetVi(time, all_primal=False), run1.GetVi(time, all_primal=False)
        )


if __name__ == "__main__":
    outputpath.mkdir(parents=True, exist_ok=True)
    main()
