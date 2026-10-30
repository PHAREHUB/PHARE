"""
Every patch holding a given B face or ghost must hold the same value for it.

MHD computes the faces and edges shared by several patches redundantly, once per patch, so
any disagreement between their ghosts persists. The layout below nests level 2 one level-1
cell inside level 1: the level-2 ghosts reach past level 1 and are filled through a coarse
temporary interpolated from level 0, where the magnetic refinement must be provided with
whole coarse cells. The second case closes the y direction with physical boundaries, so the
level-0 patches on it also share ghosts beyond the domain.
"""

import unittest

import numpy as np
from ddt import data, ddt

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core import box as boxm
from pyphare.core.box import Box2D
from pyphare.pharesee.geometry import hierarchy_overlaps
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = (40, 40)
time_step = 0.001
time_step_nbr = 2
diag_dir = "phare_outputs/mhd_shared_faces"

# level 1 spans level-1 cells [20, 59], level 2 one cell less on each side
refinement_boxes = {"L0": {"B0": Box2D(10, 29)}, "L1": {"B0": Box2D(21, 58)}}


# x stays periodic, and so does the Orszag-Tang profile along it
y_boundaries = {"ylower": {"type": "reflective"}, "yupper": {"type": "open"}}


def config(**kwargs):
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=(1.0 / cells[0], 1.0 / cells[1]),
        refinement_boxes=refinement_boxes,
        smallest_patch_size=6,
        largest_patch_size=10,
        nesting_buffer=1,
        strict=True,
        hyper_mode="spatial",
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="WENOZ",
        limiter="None",
        riemann="Rusanov",
        mhd_timestepper="TVDRK3",
        hall=False,
        model_options=["MHDModel"],
        max_mhd_level=3,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        **kwargs,
    )

    # Orszag-Tang
    B0 = 1.0 / np.sqrt(4.0 * np.pi)

    def density(x, y):
        return 25.0 / (36.0 * np.pi) + 0.0 * x

    def vx(x, y):
        return -np.sin(2.0 * np.pi * y)

    def vy(x, y):
        return np.sin(2.0 * np.pi * x)

    def bx(x, y):
        return -B0 * np.sin(2.0 * np.pi * y)

    def by(x, y):
        return B0 * np.sin(4.0 * np.pi * x)

    def zero(x, y):
        return 0.0 * x

    def p(x, y):
        return 5.0 / (12.0 * np.pi) + 0.0 * x

    ph.MHDModel(density=density, vx=vx, vy=vy, vz=zero, bx=bx, by=by, bz=zero, p=p)
    ph.ElectromagDiagnostics(
        quantity="B", write_timestamps=time_step * np.arange(time_step_nbr + 1)
    )
    return sim


@ddt
class MHDSharedFacesTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    @data({}, {"domain_boundaries": y_boundaries})
    def test_shared_B_values_agree(self, kwargs):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config(**kwargs)).run().reset()
        if cpp.mpi_rank() > 0:
            return

        hier = hierarchy_from(h5_filename=diag_dir + "/EM_B.h5")
        for time in time_step * np.arange(time_step_nbr + 1):
            overlaps = hierarchy_overlaps(hier, time)
            self.assertEqual(sorted(overlaps), [0, 1, 2])
            for ilvl, level_overlaps in overlaps.items():
                for overlap in level_overlaps:
                    box, offsets = overlap["box"], overlap["offset"]
                    if "domain_boundaries" in kwargs and any(o[1] for o in offsets):
                        continue  # hierarchy_overlaps takes every direction as periodic
                    pd0, pd1 = overlap["pdatas"]
                    values = [
                        boxm.select(
                            pd.dataset,
                            boxm.amr_to_local(box, boxm.shift(pd.ghost_box, offset)),
                        )
                        for pd, offset in zip((pd0, pd1), offsets)
                    ]
                    np.testing.assert_array_equal(
                        *values, err_msg=f"level {ilvl} at t={time} on {box}"
                    )


if __name__ == "__main__":
    startMPI()
    unittest.main()
