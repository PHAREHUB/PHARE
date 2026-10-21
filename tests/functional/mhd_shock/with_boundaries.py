#!/usr/bin/env python3

from pathlib import Path

import numpy as np
from common import BX, LEFT_INIT, RIGHT_INIT, ShockTest, dx

from pyphare import cpp
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import startMPI

cf_max = max(LEFT_INIT.cf_fast, RIGHT_INIT.cf_fast)
U_INIT = 1.05 * cf_max
cfl = 0.5
time_step = cfl * dx / (U_INIT + cf_max)
nsafe_cells = 10

nsteps = 1000
final_time = nsteps * time_step
dump_freq = 100
timestamps = np.arange(0.0, final_time, dump_freq * time_step)
diag_dir = "phare_outputs/shock_super_magnetofast_inflow"


def print_case_info():
    things_to_print = {
        "left state": LEFT_INIT,
        "right state": RIGHT_INIT,
        "cf_max": cf_max,
        "dt": time_step,
    }
    if cpp.mpi_rank() == 0:
        for key, thing in things_to_print.items():
            print(f"{key} = {thing}")


def assert_inflow_holds_left_state(diag_dir):
    if cpp.mpi_rank() != 0:
        return

    diag_path = Path(diag_dir)
    if not diag_path.exists():
        # diagnostics were not produced (e.g. skipped run) -- nothing to check
        return

    run = Run(diag_dir)

    checked_any = False
    for t in timestamps:
        try:
            sf = run.GetMHDrho(t)
        except Exception:
            continue

        level0 = sf.level(0)
        if len(level0.patches) == 0:
            continue

        # finiteness check over the interior of every level-0 patch. Ghost layers are
        # excluded: at the t=0 dump the physical-boundary ghosts still hold the fill
        # sentinel, since boundary conditions are applied during the advance.
        for patch in level0.patches:
            pd = patch.patch_datas["value"]
            ng = pd.ghosts_nbr[0]
            interior = pd.dataset[ng:-ng] if ng else pd.dataset[:]
            assert np.all(np.isfinite(interior)), (
                f"non-finite density found at t={t} in patch {patch.box}"
            )

        # leftmost patch = the one touching the xlower physical boundary
        leftmost_patch = min(level0.patches, key=lambda p: p.origin[0])
        pd = leftmost_patch.patch_datas["value"]
        ng = pd.ghosts_nbr[0]
        rho_left = pd.dataset[ng]  # first interior (non-ghost) cell

        assert np.isclose(rho_left, LEFT_INIT.rho, rtol=5e-2), (
            f"inflow density at t={t} drifted from the prescribed left state: "
            f"got {rho_left}, expected ~{LEFT_INIT.rho}"
        )
        checked_any = True

    assert checked_any, "no timestamps were available to check the inflow density"


class ShockWithBoundariesTest(ShockTest):
    def test_run(self):
        return self.run_case(
            time_step,
            final_time,
            timestamps,
            diag_dir,
            x0=nsafe_cells * dx,
            u=U_INIT,
            domain_boundaries={
                "xlower": {
                    "type": "super-magnetofast-inflow",
                    "velocity": U_INIT,
                    "density": LEFT_INIT.rho,
                    "pressure": LEFT_INIT.p,
                    "B": [BX, LEFT_INIT.by, 0.0],
                },
                "xupper": {"type": "open"},
            },
        )

    def check(self, diag_dir):
        print_case_info()
        assert_inflow_holds_left_state(diag_dir)


if __name__ == "__main__":
    startMPI()
    ShockWithBoundariesTest().test_run().tearDown()
