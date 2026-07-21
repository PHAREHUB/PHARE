#!/usr/bin/env python3
#
# mhd_harris_with_boundaries
#
# Discriminating test for the outer physical electric-field boundary condition.
#
# This is the *bottom half* of the mhd_harris double current sheet: the domain
# spans only the lower half of the original box (so it contains a single current
# sheet), with *reflective* (perfectly conducting wall) boundaries applied at
# ylower / yupper (x stays periodic). At a reflective wall the factory registers
# an AntiSymmetric condition on E, whose sole trigger is the solver's
# fillElectricGhosts() call. That antisymmetric E kills the tangential electric
# field on the wall so that constrained transport keeps the wall-normal field
# (By at the y boundaries) at zero.
#
# If the electric-field ghosts are NOT filled (i.e. the E BC is never applied),
# the tangential E on the wall is left uncorrected and Faraday makes By drift
# away from zero at the y boundaries. The test asserts that the wall-normal By
# stays below a tight threshold, so it passes only when the E BC is applied.
#
# The run uses two levels with x open: tagging keeps a refined level on the
# current sheet, which spans the full x extent and thus touches both open x
# physical boundaries from t=0 through every regrid. This exercises the
# init/regrid physical ghost fill on refined levels (B regrid fallback and the
# moment init refiners' boundary fill): if those ghosts are left at the NaN
# sentinel, the first fine-level flux poisons the state, which the explicit
# no-NaN assertion below catches.

import os
import numpy as np
from pathlib import Path

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest

os.environ["PHARE_SCOPE_TIMING"] = "1"  # turn on scope timing

ph.NO_GUI()

# bottom half of the mhd_harris domain: same dx, half the y extent
cells = (160, 40)
time_step = 0.005
# long enough for the refined level to go through hundreds of regrids while
# touching the open x boundaries, short enough for CI
final_time = 5.0
timestamps = np.arange(0, final_time + time_step, final_time / 5)
diag_dir = "phare_outputs/mhd_harris_with_boundaries"

hall = True
res = False
hyper_res = True

# max |By| tolerated at the wall-normal boundary rows. The reflective E BC keeps
# this at round-off; without it By drifts orders of magnitude above.
BN_WALL_TOL = 1e-6


def config():
    L = 0.5

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=(0.375, 0.375),
        refinement="tagging",
        max_mhd_level=2,
        max_nbr_levels=2,
        hyper_resistivity=0.0,
        resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        strict=True,
        nesting_buffer=1,
        hyper_mode="spatial",
        eta=0.0,
        nu=0.02,
        gamma=5.0 / 3.0,
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        mhd_timestepper="TVDRK3",
        hall=hall,
        res=res,
        hyper_res=hyper_res,
        model_options=["MHDModel"],
        boundary_types=("physical", "physical"),
        boundary_conditions={
            # x open (free outflow) so tagging puts a refined level on the x physical
            # boundary, exercising the regrid ghost-fill path; y reflective as before.
            "xlower": {"type": "open"},
            "xupper": {"type": "open"},
            "ylower": {"type": "reflective"},
            "yupper": {"type": "reflective"},
        },
    )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def density(x, y):
        # the simulated domain is the *bottom half* of the mhd_harris box: the
        # initial condition is evaluated with the original full height so only the
        # lower current sheet (at 0.25*Ly_full) falls inside; the upper sheet
        # (0.75*Ly_full) lies above yupper.
        Ly = 2.0 * sim.simulation_domain()[1]
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.25) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.75) / L) ** 2
        )

    def vx(x, y):
        return 0.0

    def vy(x, y):
        return 0.0

    def vz(x, y):
        return 0.0

    def bx(x, y):
        Lx = sim.simulation_domain()[0]
        # the simulated domain is the *bottom half* of the mhd_harris box: the
        # initial condition is evaluated with the original full height so only the
        # lower current sheet (at 0.25*Ly_full) falls inside; the upper sheet
        # (0.75*Ly_full) lies above yupper.
        Ly = 2.0 * sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.25 * Ly
        y2 = y - 0.75 * Ly

        dBx1 = -2 * dB * y1 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBx2 = 2 * dB * y2 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        v1 = -1
        v2 = 1.0
        return v1 + (v2 - v1) * (S(y, Ly * 0.25, L) - S(y, Ly * 0.75, L)) + dBx1 + dBx2

    def by(x, y):
        Lx = sim.simulation_domain()[0]
        # the simulated domain is the *bottom half* of the mhd_harris box: the
        # initial condition is evaluated with the original full height so only the
        # lower current sheet (at 0.25*Ly_full) falls inside; the upper sheet
        # (0.75*Ly_full) lies above yupper.
        Ly = 2.0 * sim.simulation_domain()[1]
        sigma = 1.0
        dB = 0.1

        x0 = x - 0.5 * Lx
        y1 = y - 0.25 * Ly
        y2 = y - 0.75 * Ly

        dBy1 = 2 * dB * x0 * np.exp(-(x0**2 + y1**2) / (sigma) ** 2)
        dBy2 = -2 * dB * x0 * np.exp(-(x0**2 + y2**2) / (sigma) ** 2)

        return dBy1 + dBy2

    def bz(x, y):
        return 0.0

    def p(x, y):
        return 1.0 - (bx(x, y) ** 2 + by(x, y) ** 2) / 2.0

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p)

    ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

    for quantity in ["rho", "V", "P"]:
        ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

    return sim


def plot_file_for_qty(plot_dir, qty, time):
    return f"{plot_dir}/harris_{qty}_t{time}.png"


def plot(diag_dir, plot_dir):
    run = Run(diag_dir)
    for time in timestamps:
        run.GetDivB(time).plot(
            filename=plot_file_for_qty(plot_dir, "divb", time),
            plot_patches=True,
            vmin=-1e-11,
            vmax=1e-11,
        )
        run.GetRanks(time).plot(
            filename=plot_file_for_qty(plot_dir, "Ranks", time), plot_patches=True
        )
        run.GetMHDrho(time).plot(
            filename=plot_file_for_qty(plot_dir, "rho", time), plot_patches=True
        )
        for c in ["x", "y", "z"]:
            run.GetMHDV(time).plot(
                filename=plot_file_for_qty(plot_dir, f"v{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
            run.GetB(time).plot(
                filename=plot_file_for_qty(plot_dir, f"b{c}", time),
                plot_patches=True,
                qty=f"{c}",
            )
        run.GetMHDP(time).plot(
            filename=plot_file_for_qty(plot_dir, "p", time), plot_patches=True
        )
        if hall:
            run.GetJ(time).plot(
                filename=plot_file_for_qty(plot_dir, "jz", time),
                qty="z",
                plot_patches=True,
            )


def max_normal_B_at_y_walls(diag_dir, time):
    """Return max |By| in the boundary-adjacent rows at ylower and yupper."""
    run = Run(diag_dir)
    B = run.GetB(time)
    lvl = B.patch_levels[0][0]
    dy = 0.375

    # collect (y, |By|) pairs over the finite (non-ghost) interior of every patch
    Ys, Vs = [], []
    for p in lvl.patches:
        pd = p.patch_datas["y"]  # By component
        vals = np.asarray(pd.dataset[:])
        X, Y = np.meshgrid(pd.x, pd.y, indexing="ij")
        m = np.isfinite(vals)
        Ys.append(Y[m])
        Vs.append(np.abs(vals[m]))
    yall = np.concatenate(Ys)
    vall = np.concatenate(Vs)
    ymin, ymax = yall.min(), yall.max()

    near_lower = yall <= ymin + 0.5 * dy
    near_upper = yall >= ymax - 0.5 * dy
    lower = vall[near_lower] if near_lower.any() else np.array([0.0])
    upper = vall[near_upper] if near_upper.any() else np.array([0.0])
    return float(lower.max()), float(upper.max())


def count_nans(diag_dir):
    """Total NaN count over the raw B dump (all times, levels, patches, ghosts included).

    Reads the h5 file directly: the pharesee face-centered B reconstruction pads with
    NaN and would give false positives.
    """
    import h5py

    n = 0
    with h5py.File(f"{diag_dir}/EM_B.h5", "r") as f:

        def visit(_, obj):
            nonlocal n
            if isinstance(obj, h5py.Dataset) and obj.dtype.kind == "f":
                n += int(np.isnan(obj[...]).sum())

        f.visititems(visit)
    return n


class HarrisBoundariesTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(HarrisBoundariesTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(HarrisBoundariesTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_run(self):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()
        if cpp.mpi_rank() == 0:
            plot_dir = Path(f"{diag_dir}_plots") / str(cpp.mpi_size())
            plot_dir.mkdir(parents=True, exist_ok=True)
            plot(diag_dir, plot_dir)
            times = get_times_from_h5(f"{diag_dir}/EM_B.h5")
            # the refined level touches the open x boundaries: any init/regrid
            # physical ghost left at the NaN sentinel poisons the state through
            # the first fine-level flux (the wall check below masks non-finite
            # values, so it would not catch this on its own)
            self.assertEqual(
                count_nans(diag_dir),
                0,
                "NaN in B: refined-level physical-boundary ghosts not filled "
                "at init/regrid",
            )
            bn_lower, bn_upper = max_normal_B_at_y_walls(diag_dir, times[-1])
            print(
                f"max |By| at walls @ t={times[-1]:.3f}: "
                f"ylower={bn_lower:.3e} yupper={bn_upper:.3e} (tol={BN_WALL_TOL:.1e})"
            )
            self.assertLess(
                max(bn_lower, bn_upper),
                BN_WALL_TOL,
                "wall-normal By deviated from zero: reflective E boundary "
                "condition not applied (fillElectricGhosts missing?)",
            )
        cpp.mpi_barrier()
        return self


if __name__ == "__main__":
    startMPI()
    HarrisBoundariesTest().test_run().tearDown()
