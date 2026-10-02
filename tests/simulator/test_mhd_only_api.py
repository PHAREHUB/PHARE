"""
MHD-only (no hybrid model) python API.

interp_order and refined_particle_nbr are hybrid-only keys: a run without a hybrid model
omits them, and post-processing of such a run must not need them either.
"""

import unittest

import numpy as np
import pyphare.pharein as ph
from pyphare.cpp import simulator_id
from pyphare.core.box import Box
from pyphare.core.gridlayout import MHDGridLayoutFor
from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy
from pyphare.pharesee.hierarchy.hierarchy_utils import single_patch_for_LO
from pyphare.pharesee.hierarchy.patch import Patch
from pyphare.pharesee.hierarchy.patchdata import FieldData
from pyphare.pharesee.hierarchy.patchlevel import PatchLevel

from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = 20
dl = 0.1
reconstruction = "linear"


def mhd_simulation(**kwargs):
    ph.global_vars.sim = None
    return ph.Simulation(
        time_step=0.001,
        time_step_nbr=1,
        cells=cells,
        dl=dl,
        model_options=["MHDModel"],
        reconstruction=reconstruction.capitalize(),
        limiter="VanLeer",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
        **kwargs,
    )


def field_data(box, name):
    layout = MHDGridLayoutFor(
        box, [box.lower[0] * dl], [dl], reconstruction=reconstruction
    )
    pd = FieldData(
        layout, name, None, centering=["primal"], ghosts_nbr=layout.ghosts_nbr
    )
    pd.dataset = np.arange(pd.size[0], dtype=float) + box.lower[0]
    return pd


def mhd_hierarchy(sim):
    """two level-0 patches tiling the domain, as an MHD run's hierarchy would be"""
    boxes = [Box([0], [cells // 2 - 1]), Box([cells // 2], [cells - 1])]
    patches = [Patch({"rho": field_data(box, "rho")}, f"0#{i}") for i, box in enumerate(boxes)]
    hier = PatchHierarchy(
        [{0: PatchLevel(0, patches)}], Box([0], [cells - 1]), times=[0.0]
    )
    hier._sim = sim
    return hier


class MHDOnlyAPITest(unittest.TestCase):
    def test_no_hybrid_keys_needed(self):
        sim = mhd_simulation()
        self.assertEqual(sim.interp_order, 0)  # internal, C++-facing value
        self.assertEqual(sim.refined_particle_nbr, 0)

    def test_explicit_interp_order_0_is_invalid(self):
        with self.assertRaises(ValueError):
            mhd_simulation(interp_order=0)

    def test_refined_particle_nbr_is_invalid(self):
        for nbr in (0, 2):
            with self.assertRaises(ValueError):
                mhd_simulation(refined_particle_nbr=nbr)

    def test_single_patch_for_LO(self):
        """regression: this used to build a GridLayout with no ghosts_nbr, which an
        MHD-only run cannot derive from its (absent) interp_order"""
        hier = mhd_hierarchy(mhd_simulation())
        merged = single_patch_for_LO(hier)

        patches = merged.level(0).patches
        self.assertEqual(len(patches), 1)
        for patch in hier.level(0).patches:
            np.testing.assert_array_equal(
                patches[0].patch_datas["rho"][patch.box],
                patch.patch_datas["rho"][patch.box],
            )


class MHDOnlyThroughTestHelperTest(SimulatorTest):
    """The same contract as MHDOnlyAPITest, but for a deck built through
    SimulatorTest.simulation() rather than ph.Simulation() directly -- that is
    the path every tests/simulator deck takes."""

    def mhd_simulation(self, **kwargs):
        return self.simulation(
            time_step=0.001,
            time_step_nbr=1,
            cells=cells,
            dl=dl,
            model_options=["MHDModel"],
            reconstruction=reconstruction.capitalize(),
            limiter="VanLeer",
            riemann="Rusanov",
            mhd_timestepper="TVDRK2",
            diag_options={
                "format": "phareh5",
                "options": {"dir": "phare_outputs/mhd_only_api", "mode": "overwrite"},
            },
            **kwargs,
        )

    def test_helper_injects_no_hybrid_keys(self):
        """regression: the helper defaulted interp_order to 1 and forwarded it
        unconditionally, overriding pharein's derivation from model_options. An
        MHD-only deck then asked for the coupled hybrid+MHD pybind module, which
        no res/sim MHD-only line builds."""
        sim = self.mhd_simulation()
        self.assertEqual(sim.interp_order, 0)  # internal, C++-facing value
        self.assertEqual(sim.refined_particle_nbr, 0)
        self.assertEqual(simulator_id(sim), "1_Linear_VanLeer_Rusanov_false")


if __name__ == "__main__":
    unittest.main()
