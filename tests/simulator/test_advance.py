from pyphare.cpp import cpp_lib

cpp = cpp_lib()

import unittest

import numpy as np
import pyphare.core.box as boxm
from ddt import ddt
from pyphare.core.box import Box
from pyphare.core.phare_utilities import assert_fp_any_all_close, np_array_ify
from pyphare.pharein import ElectronModel, MaxwellianFluidModel
from pyphare.pharein.diagnostics import (
    ElectromagDiagnostics,
    FluidDiagnostics,
    ParticleDiagnostics,
)
from pyphare.pharein.simulation import Simulation
from pyphare.pharesee.geometry import hierarchy_overlaps, level_ghost_boxes
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
from pyphare.simulator.simulator import Simulator

from tests.diagnostic import all_timestamps
from tests.simulator import SimulatorTest, diff_boxes


@ddt
class AdvanceTestBase(SimulatorTest):
    def _density(*xyz):
        from pyphare.pharein.global_vars import sim

        hL = np.array(sim.simulation_domain()) / 2
        _ = lambda i: -((xyz[i] - hL[i]) ** 2)
        return 0.3 + np.exp(sum([_(i) for i in range(len(xyz))]))

    def getHierarchy(
        self,
        ndim,
        interp_order,
        refinement_boxes,
        qty,
        nbr_part_per_cell=100,
        density=_density,
        smallest_patch_size=None,
        largest_patch_size=20,
        cells=120,
        time_step=0.001,
        model_init={},
        dl=0.2,
        extra_diag_options={},
        time_step_nbr=1,
        timestamps=None,
        block_merging_particles=False,
        diag_outputs="",
    ):
        diag_outputs = self.unique_diag_dir_for_test_case(
            "phare_outputs/advance", ndim, interp_order, diag_outputs
        )

        from pyphare.pharein import global_vars

        global_vars.sim = None
        if smallest_patch_size is None:
            from pyphare.pharein.simulation import check_patch_size

            _, smallest_patch_size = check_patch_size(
                ndim, interp_order=interp_order, cells=cells
            )

        extra_diag_options["mode"] = "overwrite"
        extra_diag_options["dir"] = diag_outputs
        self.register_diag_dir_for_cleanup(diag_outputs)
        sim = Simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=largest_patch_size,
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            boundary_types=["periodic"] * ndim,
            cells=np_array_ify(cells, ndim),
            dl=np_array_ify(dl, ndim),
            interp_order=interp_order,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5", "options": extra_diag_options},
            strict=True,
        )

        def bx(*xyz):
            return 1.0

        def by(*xyz):
            L = sim.simulation_domain()
            _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
            return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

        def bz(*xyz):
            L = sim.simulation_domain()
            _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
            return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

        def vx(*xyz):
            L = sim.simulation_domain()
            _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
            return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

        def vy(*xyz):
            L = sim.simulation_domain()
            _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
            return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

        def vz(*xyz):
            L = sim.simulation_domain()
            _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
            return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

        def vth(*xyz):
            return 0.01 + np.zeros_like(xyz[0])

        def vthx(*xyz):
            return vth(*xyz)

        def vthy(*xyz):
            return vth(*xyz)

        def vthz(*xyz):
            return vth(*xyz)

        MaxwellianFluidModel(
            bx=bx,
            by=by,
            bz=bz,
            protons={
                "charge": 1,
                "density": density,
                "vbulkx": vx,
                "vbulky": vy,
                "vbulkz": vz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
                "nbr_part_per_cell": nbr_part_per_cell,
                "init": model_init,
            },
        )

        ElectronModel(closure="isothermal", Te=0.12)

        if timestamps is None:
            timestamps = all_timestamps(global_vars.sim)

        for quantity in ["E", "B"]:
            ElectromagDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps,
            )

        for quantity in ["density", "bulkVelocity"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps,
            )

        poplist = ["protons"]
        for pop in poplist:
            for quantity in ["density", "flux"]:
                FluidDiagnostics(
                    quantity=quantity,
                    write_timestamps=timestamps,
                    compute_timestamps=timestamps,
                    population_name=pop,
                )

            for quantity in ["domain", "levelGhost", "patchGhost"]:
                ParticleDiagnostics(
                    quantity=quantity,
                    compute_timestamps=timestamps,
                    write_timestamps=timestamps,
                    population_name=pop,
                )

        Simulator(global_vars.sim).run()

        eb_hier = None
        if qty in ["e", "eb", "fields"]:
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_E.h5", hier=eb_hier
            )
        if qty in ["b", "eb", "fields"]:
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_B.h5", hier=eb_hier
            )
        if qty in ["e", "b", "eb"]:
            return eb_hier

        is_particle_type = qty == "particles" or qty == "particles_patch_ghost"

        if is_particle_type:
            particle_hier = None

        if qty == "particles":
            particle_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_pop_protons_domain.h5"
            )
            particle_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_pop_protons_levelGhost.h5",
                hier=particle_hier,
            )

        if is_particle_type:
            particle_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_pop_protons_patchGhost.h5",
                hier=particle_hier,
            )

        if not block_merging_particles and qty == "particles":
            merge_particles(particle_hier)

        if is_particle_type:
            return particle_hier

        if qty == "fields":
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_E.h5", hier=eb_hier
            )
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_B.h5", hier=eb_hier
            )
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_density.h5", hier=eb_hier
            )
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_bulkVelocity.h5", hier=mom_hier
            )
            return mom_hier

        if qty == "moments" or qty == "fields":
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_density.h5", hier=eb_hier
            )
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_bulkVelocity.h5", hier=mom_hier
            )
            return mom_hier

    def base_test_overlaped_fields_are_equal(self, datahier, coarsest_time):
        checks = 0
        for ilvl, overlaps in hierarchy_overlaps(datahier, coarsest_time).items():
            for overlap in overlaps:
                pd1, pd2 = overlap["pdatas"]
                box = overlap["box"]
                offsets = overlap["offset"]

                self.assertEqual(pd1.quantity, pd2.quantity)

                if pd1.quantity == "field":
                    # we need to transform the AMR overlap box, which is thus
                    # (because AMR) common to both pd1 and pd2 into local index
                    # boxes that will allow to slice the data

                    # the patchData ghost box that serves as a reference box
                    # to transfrom AMR to local indexes first needs to be
                    # shifted by the overlap offset associated to it
                    # this is because the overlap box has been calculated from
                    # the intersection of possibly shifted patch data ghost boxes

                    loc_b1 = boxm.amr_to_local(
                        box, boxm.shift(pd1.ghost_box, offsets[0])
                    )
                    loc_b2 = boxm.amr_to_local(
                        box, boxm.shift(pd2.ghost_box, offsets[1])
                    )

                    data1 = pd1.dataset
                    data2 = pd2.dataset

                    if box.ndim == 1:
                        slice1 = data1[loc_b1.lower[0] : loc_b1.upper[0] + 1]
                        slice2 = data2[loc_b2.lower[0] : loc_b2.upper[0] + 1]

                    if box.ndim == 2:
                        slice1 = data1[
                            loc_b1.lower[0] : loc_b1.upper[0] + 1,
                            loc_b1.lower[1] : loc_b1.upper[1] + 1,
                        ]
                        slice2 = data2[
                            loc_b2.lower[0] : loc_b2.upper[0] + 1,
                            loc_b2.lower[1] : loc_b2.upper[1] + 1,
                        ]

                    try:
                        # empirical max absolute observed 5.2e-15
                        # https://hephaistos.lpp.polytechnique.fr/teamcity/buildConfiguration/Phare_Phare_BuildGithubPrClang/78544
                        # seems correct considering ghosts are filled with schedules
                        # involving linear/spatial interpolations and so on where
                        # rounding errors may occur.... setting atol to 5.5e-15
                        assert_fp_any_all_close(slice1, slice2, atol=5.5e-15, rtol=0)
                        checks += 1
                    except AssertionError as e:
                        print("AssertionError", pd1.name, e)
                        print(pd1.box, pd2.box)
                        print(pd1.x.mean())
                        print(pd1.y.mean())
                        print(pd2.x.mean())
                        print(pd2.y.mean())
                        print(loc_b1)
                        print(loc_b2)
                        print(coarsest_time)
                        print(slice1)
                        print(slice2)
                        print(data1[:])
                        if self.rethrow_:
                            raise e
                        return diff_boxes(slice1, slice2, box)

        return checks

    def _test_overlaped_fields_are_equal(self, datahier, time_step_nbr, time_step):
        if cpp.mpi_rank() > 0:
            return

        checks = 0
        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time = time_step_idx * time_step
            checks += self.base_test_overlaped_fields_are_equal(datahier, coarsest_time)

        self.assertGreater(checks, time_step_nbr)
        self.assertEqual(checks % (time_step_nbr + 1), 0)

    def _test_overlapped_particledatas_have_identical_particles(
        self, ndim, interp_order, refinement_boxes, ppc=100, **kwargs
    ):
        print(
            "test_overlapped_particledatas_have_identical_particles, interporder : {}".format(
                interp_order
            )
        )
        from copy import copy

        time_step_nbr = 3
        time_step = 0.001

        datahier = self.getHierarchy(
            ndim,
            interp_order,
            refinement_boxes,
            "particles",
            time_step=time_step,
            time_step_nbr=time_step_nbr,
            nbr_part_per_cell=ppc,
            **kwargs,
        )

        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time = time_step_idx * time_step

            overlaps = hierarchy_overlaps(datahier, coarsest_time)

            for ilvl, lvl in datahier.patch_levels.items():
                print("testing level {}".format(ilvl))
                for overlap in overlaps[ilvl]:
                    pd1, pd2 = overlap["pdatas"]
                    box = overlap["box"]
                    offsets = overlap["offset"]

                    self.assertEqual(pd1.quantity, pd2.quantity)

                    if "particles" in pd1.quantity:
                        # the following uses 'offset', we need to remember that offset
                        # is the quantity by which a patch has been moved to detect
                        # overlap with the other one.
                        # so shift by +offset when evaluating patch data in overlap box
                        # index space, and by -offset when we want to shift box indexes
                        # to the associated patch index space.

                        # overlap box must be shifted by -offset to select data in the patches
                        part1 = copy(
                            pd1.dataset.select(boxm.shift(box, -np.asarray(offsets[0])))
                        )
                        part2 = copy(
                            pd2.dataset.select(boxm.shift(box, -np.asarray(offsets[1])))
                        )

                        # periodic icell overlaps need shifting to be the same
                        part1.iCells = part1.iCells + offsets[0]
                        part2.iCells = part2.iCells + offsets[1]
                        self.assertEqual(part1, part2)

    def _test_L0_particle_number_conservation(self, ndim, interp_order, ppc=100):
        cells = 120
        time_step_nbr = 10
        time_step = 0.001

        n_particles = ppc * (cells**ndim)

        datahier = self.getHierarchy(
            ndim,
            interp_order,
            None,
            "particles",
            time_step=time_step,
            time_step_nbr=time_step_nbr,
            nbr_part_per_cell=ppc,
            cells=cells,
        )

        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time = time_step_idx * time_step
            n_particles_at_t = 0
            for patch in datahier.level(0, coarsest_time).patches:
                n_particles_at_t += (
                    patch.patch_datas["protons_particles"].dataset[patch.box].size()
                )
            self.assertEqual(n_particles, n_particles_at_t)

    def _test_field_coarsening_via_subcycles(
        self, dim, interp_order, refinement_boxes, **kwargs
    ):
        print(
            "test_field_coarsening_via_subcycles for dim/interp : {}/{}".format(
                dim, interp_order
            )
        )

        from pyphare.pharein import global_vars

        from tests.amr.data.field.coarsening.test_coarsen_field import coarsen

        time_step_nbr = 3

        datahier = self.getHierarchy(
            dim,
            interp_order,
            refinement_boxes,
            "fields",
            cells=60,
            time_step=0.001,
            extra_diag_options={"fine_dump_lvl_max": 10},
            time_step_nbr=time_step_nbr,
            largest_patch_size=30,
            **kwargs,
        )

        qties = ["rho"]
        qties += [f"{qty}{xyz}" for qty in ["E", "B", "V"] for xyz in ["x", "y", "z"]]
        lvl_steps = global_vars.sim.level_time_steps
        print("LEVELSTEPS === ", lvl_steps)
        assert len(lvl_steps) > 1, "this test makes no sense with only 1 level"

        finestTimeStep = lvl_steps[-1]
        secondFinestTimeStep = lvl_steps[-2]

        finest_level_step_nbr = global_vars.sim.level_step_nbr[-1]
        uniqTimes = set([0])

        for step in range(1, finest_level_step_nbr + 1):
            checkTime = datahier.format_timestamp(finestTimeStep * step)
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
            checkTime = datahier.format_timestamp(secondFinestTimeStep * step)
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

                                coarseOffset = (
                                    coarseBox.lower - coarse_pd.layout.box.lower
                                )
                                dataBox_lower = (
                                    coarseOffset + coarse_pd.layout.nbrGhostFor(qty)
                                )
                                dataBox = Box(
                                    dataBox_lower, dataBox_lower + coarseBox.shape - 1
                                )
                                afterCoarse = np.copy(coarse_pdDataset)

                                # change values that should be updated to make failure obvious
                                assert dim < 3  # update
                                if dim == 1:
                                    afterCoarse[
                                        dataBox.lower[0] : dataBox.upper[0] + 1
                                    ] = -144123
                                if dim == 2:
                                    afterCoarse[
                                        dataBox.lower[0] : dataBox.upper[0] + 1,
                                        dataBox.lower[1] : dataBox.upper[1] + 1,
                                    ] = -144123

                                coarsen(
                                    qty,
                                    coarse_pd,
                                    fine_pd,
                                    coarseBox,
                                    fine_pdDataset,
                                    afterCoarse,
                                )

                                # https://hephaistos.lpp.polytechnique.fr/teamcity/buildConfiguration/Phare_Phare_BuildGithubPrClang/78507?hideProblemsFromDependencies=false&hideTestsFromDependencies=false&expandPull+Request+Details=true&expandBuildProblemsSection=true&expandBuildChangesSection=true&showLog=78507_7510_4947&logView=flowAware
                                # raising the bar at 2e-15 since clang failed for mpi test at 1.07.... e-15
                                try:
                                    assert_fp_any_all_close(
                                        coarse_pdDataset,
                                        afterCoarse,
                                        atol=2e-15,
                                        rtol=0,
                                    )
                                except AssertionError as e:
                                    print("failing for {}".format(qty))
                                    print(np.abs(coarse_pdDataset - afterCoarse).max())
                                    print(coarse_pdDataset)
                                    print(afterCoarse)
                                    raise e

    def base_test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
        self, L0_datahier, L0L1_datahier, quantities=None
    ):
        if quantities is None:
            quantities = [f"{EM}{xyz}" for EM in ["E", "B"] for xyz in ["x", "y", "z"]]

        from pyphare.pharein import global_vars

        from tests.amr.data.field.refine.test_refine_field import (
            refine_time_interpolate,
        )

        def assert_time_in_hier(*ts):
            for t in ts:
                self.assertIn(L0L1_datahier.format_timestamp(t), L0L1_datahier.times())

        checks = 0
        ndim = global_vars.sim.ndim
        lvl_steps = global_vars.sim.level_time_steps
        assert (
            len(lvl_steps) == 2
        ), "this test is only configured for L0 -> L1 refinement comparisons"

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
        checks = 0
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
                                        assert_fp_any_all_close(
                                            fine_ghostbox_data,
                                            refinedInterpGhostBox_data,
                                            atol=1e-15,
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
                                        atol=1e-15,
                                    )
                                    checks += 1
        if len(error_boxes):
            return error_boxes
        return checks

    def _test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
        self, ndim, interp_order, refinement_boxes
    ):
        """
        This test runs two virtually identical simulations for one step.
          L0_datahier has no refined levels
          L0L1_datahier has one refined level

        This is done to compare L0 values that haven't received the coarsened values of L1 because there is no L1,
          to the level field ghost of L1 of L0L1_datahier

        The simulations are no longer comparable after the first advance, so this test cannot work beyond that.
        """
        print(
            "test_field_coarsening_via_subcycles for dim/interp : {}/{}".format(
                ndim, interp_order
            )
        )

        import random

        rando = random.randint(0, 1e10)

        def _getHier(diag_dir, boxes=[]):
            return self.getHierarchy(
                ndim,
                interp_order,
                boxes,
                "moments",
                cells=30,
                time_step_nbr=1,
                largest_patch_size=15,
                extra_diag_options={"fine_dump_lvl_max": 10},
                time_step=0.001,
                model_init={"seed": rando},
                diag_outputs=diag_dir,
            )

        L0_datahier = _getHier(f"L0_diags")
        L0L1_datahier = _getHier(f"L0L1_diags", refinement_boxes)

        quantities = [f"{EM}{xyz}" for EM in ["E", "B"] for xyz in ["x", "y", "z"]]
        checks = (
            self.base_test_field_level_ghosts_via_subcycles_and_coarser_interpolation(
                L0_datahier, L0L1_datahier, quantities
            )
        )
        self.assertGreater(checks, len(refinement_boxes["L0"]) * len(quantities))

    def base_test_domain_particles_on_refined_level(self, datahier, new_time=None):
        """
        !! test assumes only domain particle patch_datas are present !!
        """
        times = datahier.times() if new_time is None else [new_time]
        checks = 0
        for coarsest_time in times:
            for patch in datahier.level(1, coarsest_time).patches:
                for pd_key, pd in patch.patch_datas.items():
                    if pd_key.endswith("_domain"):
                        self.assertGreater(pd[pd.box].size(), 0)
                        self.assertEqual(pd[pd.box].size(), pd[pd.ghost_box].size())
                        checks += 1
        self.assertGreater(checks, 0)

    def _test_domain_particles_on_refined_level(
        self, ndim, interp_order, refinement_boxes, **kwargs
    ):
        import pyphare.pharein as ph

        time_step_nbr = 5
        time_step = 0.001

        out = "domain_particles"

        self.base_test_domain_particles_on_refined_level(
            self.getHierarchy(
                ndim,
                interp_order,
                refinement_boxes,
                qty="particles",
                time_step=time_step,
                time_step_nbr=time_step_nbr,
                block_merging_particles=True,
                **kwargs,
            )
        )


if __name__ == "__main__":
    unittest.main()
