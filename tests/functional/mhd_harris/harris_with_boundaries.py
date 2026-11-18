#!/usr/bin/env python3

import numpy as np
from common import HarrisTest, dl, full_cells

from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import startMPI

BN_WALL_TOL = 1e-6


def max_normal_B_at_y_walls(diag_dir, time):
    """Return max |By| in the boundary-adjacent rows at ylower and yupper."""
    run = Run(diag_dir)
    B = run.GetB(time)
    lvl = B.level(0, time)
    dy = dl[1]

    # collect (y, |By|) pairs over the finite (non-ghost) interior of every patch
    Ys, Vs = [], []
    for p in lvl.patches:
        pd = p.patch_datas["y"]  # By component
        vals = np.asarray(pd.dataset[:])
        _, Y = np.meshgrid(pd.x, pd.y, indexing="ij")
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


class HarrisWithBoundariesTest(HarrisTest):
    def test_run(self):
        return self.run_case(
            (full_cells[0], full_cells[1] // 2),
            final_time=5.0,
            diag_dir="phare_outputs/mhd_harris_with_boundaries",
            domain_boundaries={
                "xlower": {"type": "open"},
                "xupper": {"type": "open"},
                "ylower": {"type": "reflective"},
                "yupper": {"type": "reflective"},
            },
        )

    def check(self, diag_dir, final_time):
        super().check(diag_dir, final_time)
        times = get_times_from_h5(f"{diag_dir}/EM_B.h5")
        bn_lower, bn_upper = max_normal_B_at_y_walls(diag_dir, times[-1])
        print(
            f"max |By| at walls @ t={times[-1]:.3f}: "
            f"ylower={bn_lower:.3e} yupper={bn_upper:.3e} (tol={BN_WALL_TOL:.1e})"
        )
        self.assertLess(
            max(bn_lower, bn_upper),
            BN_WALL_TOL,
            "wall-normal By deviated from zero: reflective boundary "
            "condition not applied (fillElectricGhosts missing?)",
        )


if __name__ == "__main__":
    startMPI()
    HarrisWithBoundariesTest().test_run().tearDown()
