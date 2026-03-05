#
#  Common base class across hybrid and mhd tests
#   see
#      tests/simulator/advance/test_advance_mhd.py
#      tests/simulator/advance/test_advance_hybrid.py
#


import unittest
import numpy as np
from ddt import ddt

import pyphare.core.box as boxm
from pyphare.core.box import amr_to_local


from pyphare import cpp
import pyphare.core.box as boxm
from pyphare.core.box import Box
from pyphare.core.phare_utilities import assert_fp_any_all_close
from pyphare.pharesee.geometry import hierarchy_overlaps, level_ghost_boxes
from pyphare.pharesee.hierarchy.hierarchy import format_timestamp

from tests.simulator import SimulatorTest, diff_boxes


@ddt
class AdvanceTestBase(SimulatorTest):
    def base_test_overlaped_fields_are_equal(self, datahier, coarsest_time):
        """
        here overlaps are calculated between patches at the same level
        """
        success_test_nbr = 0
        for ilvl, overlaps in hierarchy_overlaps(datahier, coarsest_time).items():
            for overlap in overlaps:
                pd1, pd2 = overlap["pdatas"]
                ovrlp_box = overlap["box"]
                offsets = overlap["offset"]

                self.assertEqual(pd1.quantity, pd2.quantity)
                self.assertEqual(pd1.quantity, "field")

                # we need to transform the AMR overlap box, which is thus
                # (because AMR) common to both pd1 and pd2 into local index
                # boxes that will allow to slice the data

                # the patchData ghost box that serves as a reference box
                # to transfrom AMR to local indexes first needs to be
                # shifted by the overlap offset associated to it
                # this is because the overlap box has been calculated from
                # the intersection of possibly shifted patch data ghost boxes

                box_pd1 = amr_to_local(ovrlp_box, boxm.shift(pd1.ghost_box, offsets[0]))
                box_pd2 = amr_to_local(ovrlp_box, boxm.shift(pd2.ghost_box, offsets[1]))

                slice1 = boxm.select(pd1.dataset, box_pd1)
                slice2 = boxm.select(pd2.dataset, box_pd2)

                loc_b1 = boxm.amr_to_local(
                    ovrlp_box, boxm.shift(pd1.ghost_box, offsets[0])
                )
                loc_b2 = boxm.amr_to_local(
                    ovrlp_box, boxm.shift(pd2.ghost_box, offsets[1])
                )

                try:
                    # empirical max absolute observed 5.2e-15
                    # https://hephaistos.lpp.polytechnique.fr/teamcity/buildConfiguration/Phare_Phare_BuildGithubPrClang/78544
                    # seems correct considering ghosts are filled with schedules
                    # involving linear/spatial interpolations and so on where
                    # rounding errors may occur.... setting atol to 5.5e-15
                    assert_fp_any_all_close(slice1, slice2, atol=5.5e-15, rtol=0)
                    success_test_nbr += 1
                except AssertionError as e:
                    import matplotlib.pyplot as plt
                    from matplotlib.patches import Rectangle

                    box = ovrlp_box
                    if box.ndim == 1:
                        failed_i = np.where(np.abs(slice1 - slice2) > 5.5e-15)


                    if ovrlp_box.ndim == 2:
                        failed_i, failed_j = np.where(np.abs(slice1 - slice2) > 5.5e-15)

                        def makerec(lower, upper, dl, fc="none", ec="g", lw=1, ls="-"):
                            origin = (lower[0] * dl[0], lower[1] * dl[1])
                            sizex, sizey = [
                                (u - l) * d for u, l, d in zip(upper, lower, dl)
                            ]
                            print(f"makerec: {origin}, {sizex}, {sizey}")
                            return Rectangle(
                                origin, sizex, sizey, fc=fc, ec=ec, ls=ls, lw=lw
                            )

                        datahier.plot(
                            qty=pd1.name,
                            plot_patches=True,
                            filename=pd1.name + ".png",
                            patchcolors=["k", "blue"],
                        )
                        for level_idx in range(datahier.levelNbr()):
                            fig, ax = datahier.plot(
                                qty=pd1.name,
                                plot_patches=True,
                                title=f"{pd1.name} at level {level_idx}",
                                levels=(level_idx,),
                            )
                            for patch in datahier.level(level_idx).patches:
                                ax.text(
                                    patch.patch_datas[pd1.name].origin[0],
                                    patch.patch_datas[pd1.name].origin[1],
                                    patch.id,
                                )

                            # add the overlap box only on the level
                            # where the failing overlap is
                            if level_idx == ilvl:
                                ax.add_patch(
                                    makerec(
                                        ovrlp_box.lower,
                                        ovrlp_box.upper,
                                        pd1.layout.dl,
                                        fc="none",
                                        ec="r",
                                    )
                                )
                                print("making recs for ghost boxes")
                                ax.add_patch(
                                    makerec(
                                        pd1.ghost_box.lower,
                                        pd1.ghost_box.upper,
                                        pd1.layout.dl,
                                        fc="none",
                                        ec="b",
                                        ls="--",
                                        lw=2,
                                    )
                                )
                                ax.add_patch(
                                    makerec(
                                        pd2.ghost_box.lower,
                                        pd2.ghost_box.upper,
                                        pd2.layout.dl,
                                        fc="none",
                                        ec="b",
                                        ls="--",
                                        lw=2,
                                    )
                                )
                                for i, j in zip(failed_i, failed_j):
                                    x = i + pd2.ghost_box.lower[0] + box_pd2.lower[0]
                                    x *= pd2.layout.dl[0]
                                    y = j + pd2.ghost_box.lower[1] + box_pd2.lower[1]
                                    y *= pd2.layout.dl[1]
                                    ax.plot(x, y, marker="+", color="r")

                                    x = i + pd1.ghost_box.lower[0] + box_pd1.lower[0]
                                    x *= pd1.layout.dl[0]
                                    y = j + pd1.ghost_box.lower[1] + box_pd1.lower[1]
                                    y *= pd1.layout.dl[1]
                                    ax.plot(x, y, marker="o", color="r")
                                ax.set_title(
                                    f"max error: {np.abs(slice1 - slice2).max()}, min error: {np.abs(slice1[failed_i, failed_j] - slice2[failed_i, failed_j]).min()}"
                                )
                                fig.savefig(
                                    f"{pd1.name}_level_{level_idx}_box_lower{box.lower}_upper{box.upper}.png"
                                )
                    print("coarsest time: ", coarsest_time)
                    print("AssertionError", pd1.name, e)
                    print(f"overlap box {box} (shape {box.shape})")
                    print(f"offsets: {offsets}")
                    print(
                        f"pd1 ghost box {pd1.ghost_box} (shape {pd1.ghost_box.shape}) and box {pd1.box} (shape {pd1.box.shape})"
                    )
                    print(
                        f"pd2 ghost box {pd2.ghost_box} (shape {pd2.ghost_box.shape}) and box {pd2.box} (shape {pd2.box.shape})"
                    )
                    print("interp_order: ", pd1.layout.interp_order)
                    if box.ndim == 1:
                        print(f"failing cells: {failed_i}")
                    elif box.ndim == 2:
                        print(f"failing cells: {failed_i}, {failed_j}")
                    print(coarsest_time)
                    if self.rethrow_:
                        raise e

        return success_test_nbr

    def _test_overlaped_fields_are_equal(self, datahier, time_step_nbr, time_step):
        if cpp.mpi_rank() > 0:
            return

        successful_test_nbr = 0
        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time = time_step_idx * time_step
            successful_test_nbr += self.base_test_overlaped_fields_are_equal(
                datahier, coarsest_time
            )

        self.assertGreater(successful_test_nbr, time_step_nbr)
        self.assertEqual(successful_test_nbr % (time_step_nbr + 1), 0)

    def _test_field_coarsening_via_subcycles(self, dim, interp_order, **kwargs):
        print(
            "test_field_coarsening_via_subcycles for dim/interp : {}/{}".format(
                dim, interp_order
            )
        )

        from pyphare.pharein import global_vars

        from .utilities.field_coarsening import coarsen

        time_step_nbr = 3

        datahier = self.getHierarchy(
            dim,
            interp_order,
            qty="fields",
            cells=60,
            time_step=0.001,
            extra_diag_options={"fine_dump_lvl_max": 10},
            time_step_nbr=time_step_nbr,
            largest_patch_size=30,
            **kwargs,
        )

        qties = ["rho"]
        qties += [f"{qty}{xyz}" for qty in ["E", "V"] for xyz in ["x", "y", "z"]]
        lvl_steps = global_vars.sim.level_time_steps
        print("LEVELSTEPS === ", lvl_steps)
        assert len(lvl_steps) > 1, "this test makes no sense with only 1 level"

        finestTimeStep = lvl_steps[-1]
        secondFinestTimeStep = lvl_steps[-2]

        finest_level_step_nbr = global_vars.sim.level_step_nbr[-1]
        uniqTimes = set([0])

        for step in range(1, finest_level_step_nbr + 1):
            checkTime = format_timestamp(finestTimeStep * step)
            self.assertIn(checkTime, datahier.times())
            uniqTimes.add(checkTime)

        self.assertEqual(len(uniqTimes), len(datahier.time_hier.items()))

        syncSteps = global_vars.sim.level_step_nbr[-2]  # ignore finest subcycles

        # FIX THIS AFTER NO MORE REGRIDS
        #  SEE: https://github.com/PHAREHUB/PHARE/issues/400
        assert syncSteps % time_step_nbr == 0  # perfect division
        startStep = (
            int(syncSteps / time_step_nbr) + 1
        )  # skip first coarsest step due to issue 400

        for step in range(startStep, syncSteps + 1):
            checkTime = format_timestamp(secondFinestTimeStep * step)
            self.assertIn(checkTime, datahier.times())
            nLevels = datahier.levelNbr(checkTime)
            self.assertGreaterEqual(nLevels, 2)
            levelNbrs = datahier.levelNbrs(checkTime)
            finestLevelNbr = max(levelNbrs)
            coarsestLevelNbr = min(levelNbrs)

            for coarseLevelNbr in range(coarsestLevelNbr, finestLevelNbr):
                coarsePatches = datahier.level(coarseLevelNbr, checkTime).patches
                finePatches = datahier.level(coarseLevelNbr + 1, checkTime).patches

                for coarsePatch in coarsePatches:
                    for finePatch in finePatches:
                        lvlOverlap = boxm.refine(coarsePatch.box, 2) * finePatch.box
                        if lvlOverlap is not None:
                            for qty in qties:
                                coarse_pd = coarsePatch.patch_datas[qty]
                                fine_pd = finePatch.patch_datas[qty]
                                coarseBox = boxm.coarsen(lvlOverlap, 2)

                                coarse_pdDataset = coarse_pd.dataset[:]
                                fine_pdDataset = fine_pd.dataset[:]

                                # local index of the coarsened fine box in the coarse patch box
                                # this is important if the L0 patch is not the lower left, its
                                # box will have a non-zero offset.
                                coarseOffset = (
                                    coarseBox.lower - coarse_pd.layout.box.lower
                                )

                                # local index 0 is in the ghost layer, so to get domain data
                                # we need to shift by the nbr of ghosts
                                dataBox_lower = (
                                    coarseOffset + coarse_pd.layout.nbrGhostFor(qty)
                                )
                                dataBox = Box(
                                    dataBox_lower, dataBox_lower + coarseBox.shape - 1
                                )
                                afterCoarse = np.copy(coarse_pdDataset)

                                # change values that should be updated to make failure obvious
                                boxm.DataSelector(afterCoarse)[dataBox] = -144123

                                coarsen(
                                    qty,
                                    coarse_pd,
                                    fine_pd,
                                    coarseBox,
                                    fine_pdDataset,
                                    afterCoarse,
                                )

                                # precision 2e-15 from empirical testing...
                                atol = 2e-15
                                try:
                                    assert_fp_any_all_close(
                                        coarse_pdDataset,
                                        afterCoarse,
                                        atol=atol,
                                        rtol=0,
                                    )
                                except AssertionError as e:
                                    print("failing for {}".format(qty))
                                    print(checkTime)
                                    print(np.abs(coarse_pdDataset - afterCoarse).max())
                                    print(boxm.DataSelector(coarse_pdDataset)[dataBox])
                                    print(boxm.DataSelector(afterCoarse)[dataBox])
                                    print("coarseBox", coarseBox)
                                    print("dataBox", dataBox)
                                    print("coarseOffset", coarseOffset)
                                    print("dataBox_lower", dataBox_lower)
                                    print("finePatch.box", finePatch.box)
                                    print("coarsePatch.box", coarsePatch.box)
                                    print("overlap", lvlOverlap)
                                    idx = zip(
                                        *np.where(
                                            ~np.isclose(
                                                coarse_pdDataset,
                                                afterCoarse,
                                                atol=atol,
                                                rtol=0,
                                            )
                                        )
                                    )
                                    print("idx", list(idx))
                                    raise e

    def base_test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
        self, L0_datahier, L0L1_datahier, quantities=None
    ):
        """
        this function groups code extracted from _test_field_level_ghosts_via_subcycles_and_coarser_interpolation
        because also used in test_2d_2_core.py and test_2d_10_core.py
        """
        if quantities is None:
            quantities = ["rho", "Vx", "Vy", "Vz"]

        from pyphare.pharein import global_vars

        from .utilities.test_refine_field import refine_time_interpolate

        def assert_time_in_hier(*ts):
            for t in ts:
                self.assertIn(format_timestamp(t), L0L1_datahier.times())

        successful_test_nbr = 0
        ndim = global_vars.sim.ndim
        lvl_steps = global_vars.sim.level_time_steps
        assert len(lvl_steps) == 2, (
            "this test is only configured for L0 -> L1 refinement comparisons"
        )

        coarse_ilvl = 0
        fine_ilvl = 1
        coarsest_time_before = 0  # init
        coarsest_time_after = coarsest_time_before + lvl_steps[coarse_ilvl]
        assert_time_in_hier(coarsest_time_before, coarsest_time_after)

        fine_subcycle_times = []
        for fine_subcycle in range(global_vars.sim.level_step_nbr[fine_ilvl] + 1):
            fine_subcycle_time = coarsest_time_before + (
                lvl_steps[fine_ilvl] * fine_subcycle
            )
            assert_time_in_hier(fine_subcycle_time)
            fine_subcycle_times += [fine_subcycle_time]

        interpolated_fields = refine_time_interpolate(
            L0_datahier,
            quantities,
            coarse_ilvl,
            coarsest_time_before,
            coarsest_time_after,
            fine_subcycle_times,
        )

        error_boxes = []
        successful_test_nbr = 0
        for fine_subcycle_time in fine_subcycle_times:
            fine_level_qty_ghost_boxes = level_ghost_boxes(
                L0L1_datahier, quantities, fine_ilvl, fine_subcycle_time
            )
            for qty in quantities:
                for fine_level_ghost_box_data in fine_level_qty_ghost_boxes[qty]:
                    fine_subcycle_pd = fine_level_ghost_box_data["pdata"]

                    for fine_level_ghost_box_info in fine_level_ghost_box_data["boxes"]:
                        # trim the border level ghost nodes from the primal fields to ignore them in comparison checks
                        fine_level_ghost_boxes = fine_level_ghost_box_info - boxm.grow(
                            fine_subcycle_pd.box, fine_subcycle_pd.primal_directions()
                        )

                        if ndim == 1:
                            self.assertEqual(
                                len(fine_level_ghost_boxes), 1
                            )  # should not be possibly > 1 in 1d
                            np.testing.assert_equal(
                                fine_level_ghost_boxes[0].shape,
                                fine_level_ghost_box_info.shape
                                - fine_subcycle_pd.primal_directions(),
                            )

                        for fine_level_ghost_box in fine_level_ghost_boxes:
                            upper_dims = (
                                fine_level_ghost_box.lower > fine_subcycle_pd.box.upper
                            )
                            for refinedInterpolatedField in interpolated_fields[qty][
                                fine_subcycle_time
                            ]:
                                lvlOverlap = (
                                    refinedInterpolatedField.box * fine_level_ghost_box
                                )
                                if lvlOverlap is not None:
                                    box = lvlOverlap
                                    fine_ghostbox_data = fine_subcycle_pd[box]
                                    refinedInterpGhostBox_data = (
                                        refinedInterpolatedField[box]
                                    )
                                    assert (
                                        fine_ghostbox_data.shape
                                        == refinedInterpGhostBox_data.shape
                                    )

                                    fine_ds = fine_subcycle_pd.dataset
                                    if (
                                        ndim == 1
                                    ):  # verify selecting start/end of L1 dataset from ghost box
                                        if upper_dims[0]:
                                            assert all(
                                                fine_ghostbox_data
                                                == fine_ds[
                                                    -fine_ghostbox_data.shape[0] :
                                                ]
                                            )
                                        else:
                                            assert all(
                                                fine_ghostbox_data
                                                == fine_ds[
                                                    : fine_ghostbox_data.shape[0]
                                                ]
                                            )
                                        assert (
                                            refinedInterpGhostBox_data.shape
                                            == fine_subcycle_pd.ghosts_nbr
                                        )
                                        assert (
                                            fine_ghostbox_data.shape
                                            == fine_subcycle_pd.ghosts_nbr
                                        )

                                    try:
                                        # empirical max absolute observed < 6.5e-15
                                        atol = 6.5e-15
                                        assert_fp_any_all_close(
                                            fine_ghostbox_data,
                                            refinedInterpGhostBox_data,
                                            atol=atol,
                                            rtol=0,
                                        )
                                    except AssertionError as e:
                                        print(
                                            f"FAIL level ghost subcycle_coarsening qty {qty}",
                                            fine_ghostbox_data,
                                            refinedInterpGhostBox_data,
                                            e,
                                        )
                                        if self.rethrow_:
                                            raise e
                                        error_boxes += diff_boxes(
                                            fine_ghostbox_data,
                                            refinedInterpGhostBox_data,
                                            box,
                                            atol=atol,
                                        )
                                    successful_test_nbr += 1
        if len(error_boxes):
            return error_boxes
        return successful_test_nbr

    def _test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
        self, ndim, interp_order, refinement_boxes
    ):
        """
        This test intends to check that level ghost field values during substeps
        are indeed the result of the space and time interpolation of the coarser level values.

        This requires:
        - to dump diagnostics at substeps
        - refine spatially and temporally L0 values and compare them to L1 level ghost values

        The time interpolation needs both t and t+dt coarse values.
        However, L0 values at t+dt are not available in diags since they are written after
        L0 receives coarser values from L1, while the L0 values at t+dt used in the simulation
        to do the time interpolation are the one **before** the coarsening correction.

        To achieve the test, we thus also run a L0-only simulation witht the same exact initial
        setup, and use the t0 and t0+dt diagnostics of the L0-only run to perform the
        space and time interpolation values to compare to L1 level ghosts of the former simulation.

        This test runs two virtually identical simulations for one step.
          L0_datahier has no refined levels
          L0L1_datahier has one refined level

        The simulations are no longer comparable after the first advance, so this test cannot work beyond that.


        WARNING: this test is now skipped in nD test field advance because it is
        as it is now, working on E and B, which, since the divB correction, are not
        time refined anymore. Only n, Vi and J are time refined and the test should
        thus be changed accordingly.

        """
        print(
            "test_field_coarsening_via_subcycles for dim/interp : {}/{}".format(
                ndim, interp_order
            )
        )

        import random

        rando = random.randint(0, int(1e10))

        def _getHier(diag_dir, boxes=[]):
            return self.getHierarchy(
                ndim,
                interp_order,
                boxes,
                "moments",  # only N, Vi and J are space/time interpolated, only test moments
                cells=30,
                time_step_nbr=1,
                largest_patch_size=15,
                extra_diag_options={"fine_dump_lvl_max": 10},
                time_step=0.001,
                model_init={"seed": rando},
                diag_outputs=diag_dir,
            )

        L0_datahier = _getHier("L0_diags")
        L0L1_datahier = _getHier("L0L1_diags", refinement_boxes)

        quantities = ["rho", "Vx", "Vy", "Vz"]
        successful_test_nbr = (
            self.base_test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
                L0_datahier, L0L1_datahier, quantities
            )
        )
        self.assertGreater(
            successful_test_nbr, len(refinement_boxes["L0"]) * len(quantities)
        )



if __name__ == "__main__":
    unittest.main()
