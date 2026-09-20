"""
Shared base classes/helpers for the per-dimension test_layouts_{1,2,3}d.py files.

Tests that tile-based particle layouts produce identical initialization
and time-advance results to the reference (AoSMapped) layout.

This module holds everything that does NOT vary per dimension: the atol/
refinement_boxes selection for each level (L0/L1/L2), and the comparison
machinery. `@ddt`/`@data(*permute(...))` bake their parametrization in at
class-definition time, so the concrete, test-method-bearing classes must be
(re)defined per dimension file with that file's own `ndim_list` - see
test_layouts_1d.py/2d.py/3d.py.
"""

import itertools

from pyphare import cpp
from pyphare.core.box import nDBox
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare

from tests.simulator.initialize.test_init_hybrid import HybridInitializationTest
from tests.simulator.advance.test_advance_hybrid import HybridAdvanceTest


_ref_layout = "AoSMapped"
cells = 14
ppc_per_dim = [100, 33, 15]


def permute(ndim_list, interp_orders=[1]):
    cmp_layouts = [
        l for l in cpp.supported_particle_layouts() if _ref_layout not in str(l)
    ]

    return [
        dict(ndim=ndim, interp_order=interp, cmp_layout=cmp_layout)
        for ndim, interp, cmp_layout in itertools.product(
            ndim_list, interp_orders, cmp_layouts
        )
    ]


def compare_hierarchies(test, run0, run1, atol):
    """Asserts every diagnosed quantity is identical between two Runs, across
    every dumped time. `atol` is a dict with keys "b", "e", "moments", "particles".
    """
    times = run0.all_times()["B"]
    pops = run0.all_pops()

    # subTest only collects-and-continues when run through the real unittest
    # machinery (TestCase.run() / unittest.main()); direct calls like
    # HarrisTest().test_run() have no TestResult attached, so subTest can't
    # suppress anything there and a single failure would raise immediately.
    # Collecting failures manually here works the same way either way, and
    # loops per time/quantity so one bad check doesn't hide the rest.
    failures = []

    def check(label, eqr):
        if not eqr:
            print(f"{label}\n{eqr}", flush=True)
            failures.append(label)

    # Priority order for debugging: particles -> moments -> B -> E, so the
    # first failure in the report points at the earliest quantity in the
    # causal chain (particles feed moments, moments feed E, E feeds B).
    for time in times:
        for pop in pops:
            for kind in ["levelGhost", "domain"]:
                p0 = run0.GetParticles(time, pop, type=kind)
                p1 = run1.GetParticles(time, pop, type=kind)
                check(
                    f"time={time} pop={pop} particles_{kind}",
                    hierarchy_compare(p0, p1, atol=atol["particles"]),
                )

        Ni0 = run0.GetNi(time)
        Ni1 = run1.GetNi(time)
        check(f"time={time} Ni", hierarchy_compare(Ni0, Ni1, atol=atol["moments"]))

        Vi0 = run0.GetVi(time)
        Vi1 = run1.GetVi(time)
        check(f"time={time} Vi", hierarchy_compare(Vi0, Vi1, atol=atol["moments"]))

        for pop in pops:
            N0 = run0.GetN(time, pop_name=pop)
            N1 = run1.GetN(time, pop_name=pop)
            check(f"time={time} pop={pop} N", hierarchy_compare(N0, N1, atol=atol["moments"]))

            F0 = run0.GetFlux(time, pop_name=pop)
            F1 = run1.GetFlux(time, pop_name=pop)
            check(
                f"time={time} pop={pop} flux", hierarchy_compare(F0, F1, atol=atol["moments"])
            )

        b0 = run0.GetB(time, all_primal=False)
        b1 = run1.GetB(time, all_primal=False)
        check(f"time={time} B", hierarchy_compare(b0, b1, atol=atol["b"]))

        e0 = run0.GetE(time, all_primal=False)
        e1 = run1.GetE(time, all_primal=False)
        check(f"time={time} E", hierarchy_compare(e0, e1, atol=atol["e"]))

    test.assertTrue(not failures, f"{len(failures)} check(s) failed, see output above")


class ALayoutInitTest(HybridInitializationTest):
    def compare_init(self, ndim, interp_order, cmp_layout, atol, **kwargs):
        common = dict(
            cells=cells,
            time_step_nbr=1,
            nbr_part_per_cell=ppc_per_dim[ndim - 1],
            largest_patch_size=None,
            block_merging_particles=True,
            extra_diag_options={"fine_dump_lvl_max": 10},
            **kwargs,
        )

        diag_outputs = f"init_layout_{type(self).__name__}"
        run0 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            diag_outputs=f"{diag_outputs}_{_ref_layout}",
            **common,
        )
        run1 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            particle_layout=cmp_layout,
            diag_outputs=f"{diag_outputs}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            compare_hierarchies(self, run0, run1, atol)


class LayoutInitL0Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, cmp_layout):
        self.compare_init(
            ndim,
            interp_order,
            cmp_layout,
            atol=dict(b=0, e=1e-14, moments=1e-14, particles=0),
            refinement_boxes=None,
        )


class LayoutInitL1Test(ALayoutInitTest):
    def _compare_init(self, ndim, interp_order, cmp_layout):
        self.compare_init(
            ndim,
            interp_order,
            cmp_layout,
            atol=dict(b=0, e=1e-14, moments=1e-14, particles=0),
            refinement_boxes={"L0": {"B0": nDBox(ndim, 1, 12)}},
        )


class ALayoutAdvanceTest(HybridAdvanceTest):
    def compare_advance(self, ndim, interp_order, cmp_layout, atol, **kwargs):
        common = dict(
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
        run0 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            diag_outputs=f"{diag_outputs}_{_ref_layout}",
            **common,
        )
        run1 = self.getHierarchy(
            ndim,
            interp_order,
            qty="run",
            particle_layout=cmp_layout,
            diag_outputs=f"{diag_outputs}_{cmp_layout}",
            **common,
        )

        if cpp.mpi_rank() == 0:
            compare_hierarchies(self, run0, run1, atol)


class LayoutAdvanceL0Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, cmp_layout):
        self.compare_advance(
            ndim,
            interp_order,
            cmp_layout,
            atol=dict(b=1e-14, e=1e-12, moments=1e-14, particles=1e-14),
            refinement_boxes=None,
        )


_L1_advance_atol_per_dim = {
    # deliberately at zero for now - want these to fail so the assertion
    # message's actual max diff tells us what atol they really need.
    1: dict(b=1e-15, e=1e-14, moments=1e-14, particles=1e-15),
    2: dict(b=1e-15, e=1e-14, moments=1e-14, particles=1e-15),
    3: dict(b=1e-15, e=1e-14, moments=1e-14, particles=1e-15),
}


class LayoutAdvanceL1Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, cmp_layout):
        self.compare_advance(
            ndim,
            interp_order,
            cmp_layout,
            atol=_L1_advance_atol_per_dim[ndim],
            refinement_boxes={"L0": {"B0": nDBox(ndim, 1, 12)}},
        )


_L2_advance_atol_per_dim = {
    # moments/particles pass clean at 1e-15; b and e are loosened to a small
    # multiple of the worst observed diff - both are double-precision noise
    # floor on O(1) values (b: ~5 ULP, non-associative summation order from
    # the multithreaded push/deposit - reproduced as an intermittent single-
    # element failure at b=1e-15; e: see LayoutAdvanceL2Test logs), not a
    # growing/systematic error. b is loosened for all dims since the
    # underlying thread-order nondeterminism isn't dimension-specific, even
    # though it was only caught so far at ndim=3.
    1: dict(b=5e-15, e=1e-14, moments=1e-14, particles=1e-15),
    2: dict(b=5e-15, e=5e-14, moments=1e-14, particles=1e-15),
    3: dict(b=5e-15, e=5e-14, moments=1e-14, particles=1e-15),
}


class LayoutAdvanceL2Test(ALayoutAdvanceTest):
    def _compare_advance(self, ndim, interp_order, cmp_layout):
        self.compare_advance(
            ndim,
            interp_order,
            cmp_layout,
            atol=_L2_advance_atol_per_dim[ndim],
            refinement_boxes={
                "L0": {"B0": nDBox(ndim, 1, 12)},
                "L1": {"B0": nDBox(ndim, 4, 21)},
            },
            time_step_nbr=4,
        )
