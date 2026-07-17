import sys
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt

import pyphare.core.box as boxm
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.pharesee.hierarchy.hierarchy_utils import (
    hierarchy_compare,
    single_patch_for_LO,
)


atol = 1e-12

outputpath = Path("phare_outputs/comapare_diags/")
outputs = str(outputpath)


def make_fig(hier, fig_name, ilvl, collections):
    hier.plot_2d_patches(0, collections=collections).savefig(fig_name + ".png")


def plot_file(time, qty):
    return outputs + f"/{qty}_{time}.png"


def plot_ghost_box_diff(hier, time, qty, name):
    patches = hier.level(time, 0)
    all_data = [p.patch_datas[qty].dataset[:] for p in patches]
    absmax = max(np.max(np.abs(d)) for d in all_data) + 1e-6
    vmin, vmax = -absmax, absmax

    fig, ax = plt.subplots()
    im = None
    for patch in patches:
        print("patch.box", patch.box)
        pdat = patch.patch_datas[qty]
        layout = pdat.layout
        im = ax.pcolormesh(
            layout.yeeCoordsFor(qty, "x", withGhosts=True),
            layout.yeeCoordsFor(qty, "y", withGhosts=True),
            pdat.dataset[:].T,
            cmap="RdBu_r",
            vmin=vmin,
            vmax=vmax,
        )
    if im is not None:
        plt.colorbar(im, ax=ax)
    fig.savefig(plot_file(time, name + "_gb"), dpi=200)


def main():
    run0 = Run(sys.argv[1])
    run1 = Run(sys.argv[2])

    times = run0.all_times()["B"]
    print("times", times)
    for time in times:
        print("\ntime :", time)

        b0 = run0.GetB(time, all_primal=False)
        b1 = run1.GetB(time, all_primal=False)
        eqr = hierarchy_compare(b0, b1, atol=atol)
        print("\n\t E", eqr)

        if not eqr:
            bdiff = VectorField(b0) - VectorField(b1)
            for c in ["x", "y", "z"]:
                bdiff.plot(
                    filename=plot_file(time, f"B{c}"),
                    qty=f"{c}",
                    plot_patches=True,
                    vmin=-1e-12,
                    vmax=1e-12,
                )

        e0 = run0.GetE(time, all_primal=False)
        e1 = run1.GetE(time, all_primal=False)
        eqr = hierarchy_compare(e0, e1, atol=atol)
        print("\n\t E", eqr)
        if not eqr:
            ediff = VectorField(e0) - VectorField(e1)
            for c in ["x", "y", "z"]:
                ediff.plot(
                    filename=plot_file(time, f"E{c}"),
                    qty=f"{c}",
                    plot_patches=True,
                    vmin=-1e-12,
                    vmax=1e-12,
                )

        N0 = run0.GetNi(time)
        N1 = run1.GetNi(time)
        eqr = hierarchy_compare(N0, N1, atol=atol)
        print(f"\n\n Ni", eqr)
        if not eqr:
            Ndiff = ScalarField(N0) - ScalarField(N1)
            Ndiff.plot(
                filename=plot_file(time, "Ni"),
                plot_patches=True,
                vmin=-1e-12,
                vmax=1e-12,
            )

        V0 = run0.GetVi(time, all_primal=False)
        V1 = run1.GetVi(time, all_primal=False)
        eqr = hierarchy_compare(V0, V1, atol=atol)
        print(f"\n\n V ", eqr)
        if not eqr:
            Vdiff = VectorField(V0) - VectorField(V1)
            for c in ["x", "y", "z"]:
                Vdiff.plot(
                    filename=plot_file(time, f"V{c}"),
                    qty=f"{c}",
                    plot_patches=True,
                    vmin=-1e-12,
                    vmax=1e-12,
                )
                plot_ghost_box_diff(Vdiff, time, c, f"V{c}")

        pops = run0.all_pops()
        for pop in pops:
            N0 = run0.GetN(time, pop_name=pop)
            N1 = run1.GetN(time, pop_name=pop)
            eqr = hierarchy_compare(N0, N1, atol=atol)
            print(f"\n\n rho {pop}", eqr)
            if not eqr:
                Ndiff = ScalarField(N0) - ScalarField(N1)
                Ndiff.plot(
                    filename=plot_file(time, "N"),
                    plot_patches=True,
                    vmin=0,
                    vmax=1e-12,
                )
                plot_ghost_box_diff(Ndiff, time, "value", "N")

        for pop in pops:
            F0 = run0.GetFlux(time, pop_name=pop)
            F1 = run1.GetFlux(time, pop_name=pop)
            eqr = hierarchy_compare(F0, F1, atol=atol)
            print(f"\n\n flux {pop}", eqr)
            if not eqr:
                Fdiff = VectorField(F0) - VectorField(F1)
                for c in ["x", "y", "z"]:
                    Fdiff.plot(
                        filename=plot_file(time, f"F{c}"),
                        qty=c,
                        plot_patches=True,
                        vmin=0,
                        vmax=1e-12,
                    )
                    plot_ghost_box_diff(Fdiff, time, c, f"F{c}")


if __name__ == "__main__":
    outputpath.mkdir(parents=True, exist_ok=True)
    main()
