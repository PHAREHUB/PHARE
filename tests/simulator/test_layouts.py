#!/usr/bin/env python3
"""
Tests that tile-based particle layouts produce identical initialization
and time-advance results to the reference (AoSMapped) layout.
"""

import unittest
import itertools
from ddt import data, ddt, unpack

from pyphare import cpp
from pyphare.core.box import nDBox
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare

from tests.simulator.initialize.test_init_hybrid import HybridInitializationTest
from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest


ndim_list = [1]
interp_orders = [1]
_ref_layout = "AoSMapped"
cells = 32
ppc_per_dim = [100, 66, 20]
# os.environ["PHARE_TILING_MIN_BEFORE_SPLIT"] = "1000"


def permute():
    cmp_layouts = [
        l for l in cpp.supported_particle_layouts() if _ref_layout not in str(l)
    ]

    return [
        dict(ndim=ndim, interp_order=interp, cmp_layout=cmp_layout)
        for ndim, interp, cmp_layout in itertools.product(
            ndim_list, interp_orders, cmp_layouts
        )
    ]


class ALayoutInitTest(HybridInitializationTest):
    def compare_init(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        common = dict(
            qty=qty,
            cells=cells,
            time_step_nbr=1,
            nbr_part_per_cell=ppc_per_dim[ndim - 1],
            largest_patch_size=None,
            block_merging_particles=True,
            extra_diag_options={"fine_dump_lvl_max": 10},
            **kwargs,
        )

        diag_outputs = f"init_layout_{type(self).__name__}"
        ref_hier = self.getHierarchy(
            ndim,
            interp_order,
            diag_outputs=f"{diag_outputs}_{_ref_layout}",
            **common,
        )
        cmp_hier = self.getHierarchy(
            ndim,
            interp_order,
            sim_setup_kwargs=dict(layout=cmp_layout),
            diag_outputs=f"{diag_outputs}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            eqr = hierarchy_compare(ref_hier, cmp_hier, atol=atol)
            print(eqr)
            self.assertTrue(eqr)


@ddt
class LayoutInitL0Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, qty, cmp_layout, atol=0):
        super().compare_init(
            ndim, interp_order, qty, cmp_layout, atol=atol, refinement_boxes=None
        )

    @data(*permute())
    @unpack
    def test_B_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "b", cmp_layout, atol=0)

    @data(*permute())
    @unpack
    def test_E_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "e", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_moments_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "moments", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_particles_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "particles", cmp_layout, atol=0)


@ddt
class LayoutInitL1Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, qty, cmp_layout, atol=0):
        super().compare_init(
            ndim,
            interp_order,
            qty,
            cmp_layout,
            atol=atol,
            refinement_boxes={"L0": {"B0": nDBox(ndim, 5, 14)}},
        )

    @data(*permute())
    @unpack
    def test_B_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "b", cmp_layout, atol=0)

    @data(*permute())
    @unpack
    def test_E_init_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "e", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_moments_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "moments", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_particles_init_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_init(ndim, interp_order, "particles", cmp_layout, atol=0)


class ALayoutAdvanceTest(HybridAdvanceTest):
    def compare_advance(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        common = dict(
            qty=qty,
            cells=cells,
            time_step_nbr=kwargs.pop("time_step_nbr", 1),
            model_init={"seed": 1337},
            nbr_part_per_cell=ppc_per_dim[ndim - 1],
            largest_patch_size=None,
            block_merging_particles=True,
            extra_diag_options={"fine_dump_lvl_max": 10},
            **kwargs,
        )

        diag_outputs = f"adv_layout_{type(self).__name__}"
        ref_hier = self.getHierarchy(
            ndim,
            interp_order,
            diag_outputs=f"{diag_outputs}_{_ref_layout}",
            **common,
        )
        cmp_hier = self.getHierarchy(
            ndim,
            interp_order,
            sim_setup_kwargs=dict(layout=cmp_layout),
            diag_outputs=f"{diag_outputs}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            if kwargs.get("as_run", False):
                return ref_hier, cmp_hier
            eqr = hierarchy_compare(ref_hier, cmp_hier, atol=atol)
            print(eqr)
            self.assertTrue(eqr)


@ddt
class LayoutAdvanceL0Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        return super().compare_advance(
            ndim, interp_order, qty, cmp_layout, atol=atol, refinement_boxes=None
        )

    @data(*permute())
    @unpack
    def test_B_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "b", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_E_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "e", cmp_layout, atol=1e-12)

    @data(*permute())
    @unpack
    def test_moments_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "moments", cmp_layout, atol=1e-14)

    @data(*permute())
    @unpack
    def test_particles_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "particles", cmp_layout, atol=1e-14)


@ddt
class LayoutAdvanceL1Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        super().compare_advance(
            ndim,
            interp_order,
            qty,
            cmp_layout,
            atol=atol,
            refinement_boxes={"L0": {"B0": nDBox(ndim, 5, 14)}},
        )

    @data(*permute())
    @unpack
    def test_B_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "b", cmp_layout, atol=5e-15)

    @data(*permute())
    @unpack
    def test_E_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "e", cmp_layout, atol=8e-15)

    @data(*permute())
    @unpack
    def test_moments_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "moments", cmp_layout, atol=5e-15)

    @data(*permute())
    @unpack
    def test_particles_advance_matches_reference_layout(
        self, ndim, interp_order, cmp_layout
    ):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        self._compare_advance(ndim, interp_order, "particles", cmp_layout, atol=5e-15)


@ddt
class LayoutAdvanceL2Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, qty, cmp_layout, atol=0, **kwargs):
        return super().compare_advance(
            ndim,
            interp_order,
            qty,
            cmp_layout,
            atol=atol,
            refinement_boxes={
                "L0": {"B0": nDBox(ndim, 5, 24)},
                "L1": {"B0": nDBox(ndim, 14, 40)},
            },
            time_step_nbr=10,
            as_run=True,
        )

    @data(*permute())
    @unpack
    def test_advance_matches_reference_layout(self, ndim, interp_order, cmp_layout):
        print(
            f"{self._testMethodName} ndim={ndim} interp={interp_order} layout={cmp_layout}"
        )
        atol = 5e-15
        run0, run1 = self._compare_advance(
            ndim, interp_order, "b", cmp_layout, atol=atol
        )
        times = run0.all_times()["B"]
        print("times", times)
        for time in times:
            print("\ntime :", time)

            b0 = run0.GetB(time, all_primal=False)
            b1 = run1.GetB(time, all_primal=False)
            eqr = hierarchy_compare(b0, b1, atol=atol)
            print("\n\t E", eqr)

            e0 = run0.GetE(time, all_primal=False)
            e1 = run1.GetE(time, all_primal=False)
            eqr = hierarchy_compare(e0, e1, atol=atol)
            print("\n\t E", eqr)

            N0 = run0.GetNi(time)
            N1 = run1.GetNi(time)
            eqr = hierarchy_compare(N0, N1, atol=atol)
            print(f"\n\n Ni", eqr)

            V0 = run0.GetVi(time, all_primal=False)
            V1 = run1.GetVi(time, all_primal=False)
            eqr = hierarchy_compare(V0, V1, atol=atol)
            print(f"\n\n V ", eqr)

            pops = run0.all_pops()
            for pop in pops:
                N0 = run0.GetN(time, pop_name=pop)
                N1 = run1.GetN(time, pop_name=pop)
                eqr = hierarchy_compare(N0, N1, atol=atol)
                print(f"\n\n rho {pop}", eqr)

            for pop in pops:
                F0 = run0.GetFlux(time, pop_name=pop)
                F1 = run1.GetFlux(time, pop_name=pop)
                eqr = hierarchy_compare(F0, F1, atol=atol)
                print(f"\n\n flux {pop}", eqr)


if __name__ == "__main__":
    from pyphare.simulator.simulator import startMPI

    # raise "?"
    startMPI()
    unittest.main()
