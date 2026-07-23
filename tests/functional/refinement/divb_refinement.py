#!/usr/bin/env python3
# divB e2e for the higher-order field-refinement operators (refinement_order 0/2).
#
# Idea: in the discrete Yee scheme Faraday preserves divB exactly, so on a fine level
# created by B prolongation, max|divB| is set *purely* by the B refinement at the
# coarse-fine boundary. We init B from a Harris double current sheet rotated so that
#   Bx = Bx(y)  (double tanh),  By = 0,  Bz = 0
# => divB = dBx/dx + dBy/dy = 0 to MACHINE PRECISION (Bx is constant in x, By identically
# 0). Unlike a point-sampled 2D curl (which carries a discrete common-mode interior divB
# ~3e-3 identical across all orders, so the test discriminates nothing), this init's divB is
# machine-zero analytically; after NSTEPS of evolution + diagnostic reconstruction the floor
# is a small truncation residual (~2e-7) that is itself ~equal on every level, so the fine
# level only inherits it and any operator shared-face spike (O(B/dx)) is unmistakable.
# The y-domain is tall enough that the tanh tails are flat at the periodic y boundary
# (sheets at 0.3Ly/0.7Ly, L=0.5 => sheet-to-boundary dist/L = 24); otherwise the B mismatch
# across the periodic seam would itself create real boundary divB.
# A divB-preserving refinement keeps fine-level max|divB| at the (machine-zero) coarse floor
# for every order; a misclassified / overwritten shared face spikes it to O(B/dx). The B
# tangential (dual) child correction is sign-antisymmetric about the coarse face value =>
# zero-mean per coarse face => divB preserved. This test guards that property end-to-end.
#
# Two modes, both 2D (divB is identically 0 in 1D) and both on the SAME Harris sheet:
#   boxes   - deterministic fine level via refinement_boxes (static C-F boundary straddling
#             a sheet, box edges sitting in the flat field away from the sheet)
#   tagging - tagger-created level that regrids over time; the sharp Bx(y) sheet triggers the
#             2D default tagger (max finite-diff ratio over Bx,By,Bz) and exercises the regrid
#             path with a div-free field
#
# Run: mpirun -n 12 python -u divb_refinement.py [boxes|tagging] [orders...]
#   orders    : any of 0 2            (default 0 2)

import sys
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

ph.NO_GUI()
startMPI()

# Harris double sheet. Bx = Bx(y) only => divB = dBx/dx (=0, Bx const in x) + dBy/dy (=0,
# By identically 0) = 0 to machine precision on every level. y must be tall enough that the
# tanh tails are flat at the periodic seam: sheets at 0.3Ly/0.7Ly, L=0.5 => dist/L = 24.
cells = (20, 80)
dl = (0.5, 0.5)  # domain 10 x 40
Lx, Ly = cells[0] * dl[0], cells[1] * dl[1]
L = 0.5  # sheet half-width
V1, V2 = -1.0, 1.0  # asymptotic Bx outside / between the sheets
K = 0.7  # total-pressure constant (P = (K - B^2/2)/n), guarantees T>0

time_step = 0.002
NSTEPS = {"boxes": 2, "tagging": 20}  # tagging needs steps for regrids to fire

# Interior fine box straddling the LOWER sheet (y0 = 0.3*Ly = 12 => cell 24): C-F boundaries
# in y sit in the flat field a few cells off the sheet, the sheet itself is fully refined.
FINE_BOX = [[4, 16], [15, 31]]

# A divB spike from a misclassified/overwritten shared face is O(B/dx) ~ 0.1. The refinement
# is divB-preserving iff the fine level does not AMPLIFY divB beyond the floor it inherits.
# With the Harris Bx(y) init divB is machine-zero analytically (Bx const in x, By identically
# 0); the evolved+reconstructed coarse floor is ~2e-7 for both modes, and the fine level sits
# at ~2x that (a refinement-boundary truncation factor), NOT at the O(B/dx)~0.1 a misclassified
# shared face would give. The contract: fine <= REL_TOL * max(coarse-floor, ABS_TOL) AND fine
# within REL_TOL of the order-0 baseline (isolates the operator from the init floor). ABS_TOL
# keeps a non-amplifying operator from tripping on the small coarse-floor noise.
ABS_TOL = 1e-4  # floor below which divB is "already zero" (truncation, not a bug)
REL_TOL = 5.0   # fine must not exceed REL_TOL x (coarse floor) nor REL_TOL x (order-0 fine)


# Shared Harris double-sheet init (both modes). Bx is a function of y ONLY (double tanh,
# sheets at 0.3Ly and 0.7Ly), By = Bz_tangential = 0 => divB = 0 to machine precision.
def _S(y, y0, l):
    return 0.5 * (1.0 + np.tanh((y - y0) / l))


def bx_harris(x, y):
    return V1 + (V2 - V1) * (_S(y, 0.3 * Ly, L) - _S(y, 0.7 * Ly, L)) + 0.0 * x


def by_harris(x, y):
    return 0.0 * x + 0.0 * y


def config(mode, order, diag_dir):
    nsteps = NSTEPS[mode]
    check_time = nsteps * time_step
    common = dict(
        time_step=time_step,
        time_step_nbr=nsteps,
        cells=cells,
        dl=dl,
        boundary_types=["periodic", "periodic"],
        refinement_order=order,  # <-- path under test
        strict=True,
        nesting_buffer=0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )
    if mode == "boxes":
        sim = ph.Simulation(
            refinement="boxes",
            refinement_boxes={"L0": {"B0": FINE_BOX}},
            smallest_patch_size=10,
            largest_patch_size=16,
            **common,
        )
    else:
        sim = ph.Simulation(
            refinement="tagging",
            max_nbr_levels=2,
            tag_buffer=1,
            **common,
        )

    def bz(x, y):
        return 0.0 * x

    # Harris pressure balance: total pressure P + B^2/2 = K is uniform, T = P/n > 0.
    def density(x, y):
        return (
            0.4
            + 1.0 / np.cosh((y - 0.3 * Ly) / L) ** 2
            + 1.0 / np.cosh((y - 0.7 * Ly) / L) ** 2
        )

    def temperature(x, y):
        b2 = bx_harris(x, y) ** 2 + by_harris(x, y) ** 2 + bz(x, y) ** 2
        t = (K - 0.5 * b2) / density(x, y)
        assert np.all(t > 0)
        return t

    def thermal(x, y):
        return np.sqrt(temperature(x, y))

    def zero(x, y):
        return 0.0 * x

    ph.MaxwellianFluidModel(
        bx=bx_harris,
        by=by_harris,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": thermal,
            "vthy": thermal,
            "vthz": thermal,
            "nbr_part_per_cell": 30,
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.ElectromagDiagnostics(quantity="B", write_timestamps=np.array([check_time]))
    return sim, check_time


def max_divb_per_level(diag_dir, check_time):
    # IMPORTANT: measure the DOMAIN INTERIOR only. _compute_divB builds divB from the full
    # B datasets (ghosts included) and drops ghosts_nbr, so dataset[:] still spans the ghost
    # band. The fine-level coarse-fine ghost fill is a known divB hot-spot (order-independent,
    # pre-existing) that is NOT a physical interior violation -> strip a ghost margin.
    run = Run(diag_dir)
    B = run.GetB(check_time, all_primal=False)
    blvls = B.levels(check_time)
    ng = int(blvls[min(blvls)].patches[0].patch_datas["Bx"].ghosts_nbr[0])
    divb = run.GetDivB(check_time)
    lvls = divb.levels(check_time)
    out = {}
    for lvl, level in lvls.items():
        m = 0.0
        for patch in level.patches:
            arr = np.abs(patch.patch_datas["value"].dataset[:])
            if all(s > 2 * ng for s in arr.shape):
                arr = arr[tuple(slice(ng, -ng) for _ in arr.shape)]
            if arr.size:
                m = max(m, float(np.nanmax(arr)))
        out[lvl] = m
    return out


def run_variant(mode, order):
    diag_dir = f"divb_{mode}_o{order}"
    sim, check_time = config(mode, order, diag_dir)
    label = f"order={order}"
    if cpp.mpi_rank() == 0:
        print(f"=== divB {mode} {label} ===", flush=True)
    Simulator(sim).run().reset()
    ph.global_vars.sim = None
    if cpp.mpi_rank() == 0:
        per = max_divb_per_level(diag_dir, check_time)
        finest = max(per.keys())
        coarse = per[min(per.keys())]
        print(
            f"  {label}: levels={sorted(per)}  "
            f"max|divB|_fine={per[finest]:.3e}  max|divB|_coarse={coarse:.3e}",
            flush=True,
        )
        return per
    return None


if __name__ == "__main__":
    argv = sys.argv[1:]
    mode = argv[0] if argv and argv[0] in ("boxes", "tagging") else "boxes"
    orders = [int(a) for a in argv if a.isdigit()] or [0, 2]

    res = {o: run_variant(mode, o) for o in orders}
    if cpp.mpi_rank() != 0:
        sys.exit(0)

    print(f"=== divB {mode} summary ===", flush=True)
    # baseline = legacy order-0 fine-level divB (operator-vs-legacy isolation). Falls back to
    # the first order's fine divB if order 0 was not requested.
    base_per = res.get(0, res[orders[0]])
    base = base_per[max(base_per.keys())]
    ok = True
    for order in orders:
        per = res[order]
        label = f"order={order}"
        finest = max(per.keys())
        m = per[finest]
        if finest == min(per.keys()):
            ok = False
            print(f"  {label}: NO FINE LEVEL FORMED  [FAIL]", flush=True)
            continue
        coarse = per[min(per.keys())]
        ceiling = REL_TOL * max(coarse, ABS_TOL)  # not amplified beyond the inherited floor
        floor_ok = m <= ceiling
        rel_ok = m <= REL_TOL * base if base else True  # operator vs legacy baseline
        status = "OK" if (floor_ok and rel_ok) else "FAIL"
        if not (floor_ok and rel_ok):
            ok = False
        print(
            f"  {label}: fine={m:.3e}  coarse={coarse:.3e}  base={base:.3e}  "
            f"fine/coarse={m / coarse if coarse else float('inf'):.2f}  "
            f"fine/base={m / base if base else float('inf'):.2f}  [{status}]",
            flush=True,
        )
    print(f"DIVB_{mode.upper()}_OK" if ok else f"DIVB_{mode.upper()}_FAIL", flush=True)
    sys.exit(0 if ok else 1)
