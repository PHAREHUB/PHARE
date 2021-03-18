
from pybindlibs import cpp

from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein.diagnostics import ParticleDiagnostics, FluidDiagnostics, ElectromagDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.pharein.simulation import Simulation
from pyphare.pharesee.geometry import level_ghost_boxes, hierarchy_overlaps
from pyphare.core.gridlayout import GridLayout, yee_element_is_primal
from pyphare.pharesee.particles import aggregate as aggregate_particles
import pyphare.core.box as boxm
from pyphare.core.box import Box, Box1D
import numpy as np
import unittest
from ddt import ddt, data, unpack


# AdvanceTest.test_field_level_ghosts_via_subcycles_and_coarser_interpolation_1

@ddt
class AdvanceTest(unittest.TestCase):

    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]

    def getHierarchy(self, interp_order, refinement_boxes, qty, nbr_part_per_cell=100,
                     diag_outputs="phare_outputs",
                     smallest_patch_size=5, largest_patch_size=5,
                     cells= 120, time_step=0.001,
                     dl=0.3, extra_diag_options={}, time_step_nbr=1, timestamps=None):

        from pyphare.pharein import global_vars
        global_vars.sim = None
        startMPI()
        extra_diag_options["mode"] = "overwrite"
        extra_diag_options["dir"] = diag_outputs
        Simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=largest_patch_size,
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            boundary_types="periodic",
            cells=cells,
            dl=dl,
            interp_order=interp_order,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5",
                          "options": extra_diag_options}
        )

        def density(x):
            return 1.

        def S(x,x0,l):
            return 0.5*(1+np.tanh((x-x0)/l))

        def bx(x):
            return 1.

        def by(x):
            L = global_vars.sim.simulation_domain()[0]
            v1=-1
            v2=1.
            return v1 + (v2-v1)*(S(x,L*0.25,1) -S(x, L*0.75, 1))

        def bz(x):
            return 0.5

        def b2(x):
            return bx(x)**2 + by(x)**2 + bz(x)**2

        def T(x):
            K = 1
            return 1/density(x)*(K - b2(x)*0.5)

        def vx(x):
            return 0.

        def vy(x):
            return 0.

        def vz(x):
            return 0.

        def vthx(x):
            return T(x)

        def vthy(x):
            return T(x)

        def vthz(x):
            return T(x)


        MaxwellianFluidModel(bx=bx, by=by, bz=bz,
                             protons={"charge": 1,
                                      "density": density,
                                      "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
                                      "vthx": vthx, "vthy": vthy, "vthz": vthz,
                                      "nbr_part_per_cell": nbr_part_per_cell,
                                      "init": {"seed": 1337}})

        ElectronModel(closure="isothermal", Te=0.12)

        if timestamps is None:
            timestamps = np.arange(0, global_vars.sim.final_time + global_vars.sim.time_step, global_vars.sim.time_step)

        for quantity in ["E", "B"]:
            ElectromagDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps
            )

        for quantity in ["density", "bulkVelocity"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps
            )

        for pop in ["protons"]:
            for quantity in ["density", "flux"]:
                FluidDiagnostics(quantity=quantity,
                                 write_timestamps=timestamps,
                                 compute_timestamps=timestamps,
                                 population_name=pop)

        Simulator(global_vars.sim).initialize().run().reset()

        eb_hier = None
        if qty in ["e", "eb"]:
            eb_hier = hierarchy_from(h5_filename=diag_outputs+"/EM_E.h5", hier=eb_hier)
        if qty in ["b", "eb"]:
            eb_hier = hierarchy_from(h5_filename=diag_outputs+"/EM_B.h5", hier=eb_hier)
        if qty in ["e", "b", "eb"]:
            return eb_hier

        if qty == "moments":
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_density.h5")
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_bulkVelocity.h5", hier=mom_hier)
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_density.h5", hier=mom_hier)
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_flux.h5", hier=mom_hier)
            return mom_hier




    def _test_field_coarsening_via_subcycles(self, dim, interp_order, refinement_boxes):
        print("test_field_coarsening_via_subcycles for dim/interp : {}/{}".format(dim, interp_order))

        from tests.amr.data.field.coarsening.test_coarsen_field import coarsen
        from pyphare.pharein import global_vars

        time_step_nbr=3

        diag_outputs=f"phare_outputs_subcycle_coarsening_{self.ddt_test_id()}"
        datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", cells=30,
                                      diag_outputs=diag_outputs, time_step=0.001,
                                      extra_diag_options={"fine_dump_lvl_max": 10},
                                      time_step_nbr=time_step_nbr, smallest_patch_size=5,
                                      largest_patch_size=30)

        lvlSteps = global_vars.sim.level_time_steps
        print("LEVELSTEPS === ", lvlSteps)
        assert len(lvlSteps) > 1  # this test makes no sense with only 1 level

        finestTimeStep = lvlSteps[-1]
        secondFinestTimeStep = lvlSteps[-2]

        finest_level_step_nbr = global_vars.sim.level_step_nbr[-1]
        uniqTimes = set([0])

        for step in range(1, finest_level_step_nbr + 1):
            checkTime = datahier.format_timestamp(finestTimeStep * step)
            self.assertIn(checkTime, datahier.times())
            uniqTimes.add(checkTime)

        self.assertEqual(len(uniqTimes), len(datahier.time_hier.items()))

        syncSteps = global_vars.sim.level_step_nbr[-2] # ignore finest subcycles

        # FIX THIS AFTER NO MORE REGRIDS
        #  SEE: https://github.com/PHAREHUB/PHARE/issues/400
        assert syncSteps % time_step_nbr == 0 # perfect division
        startStep = int(syncSteps / time_step_nbr) + 1 # skip first coarsest step due to issue 400

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
                            for EM in ["E", "B"]:
                                for xyz in ["x", "y", "z"]:
                                    qty = f"{EM}{xyz}"
                                    coarse_pd = coarsePatch.patch_datas[qty]
                                    fine_pd  = finePatch.patch_datas[qty]
                                    coarseBox = boxm.coarsen(lvlOverlap, 2)

                                    nGhosts = coarse_pd.layout.nbrGhostFor(qty)

                                    coarse_pdDataset = coarse_pd.dataset[:]
                                    fine_pdDataset = fine_pd.dataset[:]

                                    coarseOffset = coarseBox.lower - coarse_pd.layout.box.lower
                                    dataBox_lower = coarseOffset + nGhosts
                                    dataBox = Box(dataBox_lower, dataBox_lower + coarseBox.shape - 1)

                                    afterCoarse = np.copy(coarse_pdDataset)
                                    afterCoarse[dataBox.lower[0] : dataBox.upper[0] + 1] = 0

                                    coarsen(qty, coarse_pd.layout, fine_pd.layout, coarseBox, fine_pdDataset, afterCoarse)
                                    np.testing.assert_allclose(coarse_pdDataset, afterCoarse, atol=1e-6)


    @data(
       ({"L0": {"B0": Box1D(10, 19)}}),
       ({"L0": {"B0": Box1D(10, 14), "B1": Box1D(15, 19)}}),
       ({"L0": {"B0": Box1D(6, 23)}}),
       ({"L0": {"B0": Box1D( 2, 12), "B1": Box1D(13, 25)}}),
       ({"L0": {"B0": Box1D( 5, 20)}, "L1": {"B0": Box1D(15, 19)}}),
       ({"L0": {"B0": Box1D( 5, 20)}, "L1": {"B0": Box1D(12, 38)}, "L2": {"B0": Box1D(30, 52)} }),
    )
    def test_field_coarsening_via_subcycles(self, refinement_boxes):
        dim = refinement_boxes["L0"]["B0"].ndim
        self._test_field_coarsening_via_subcycles(dim, interp_order=1, refinement_boxes=refinement_boxes)




    def _test_field_level_ghosts_via_subcycles_and_coarser_interpolation(self, ndim, interp_order, refinement_boxes, refinement_ratio=2):
        print("test_field_coarsening_via_subcycles for dim/interp : {}/{}".format(ndim, interp_order))

        from tests.amr.data.field.refine.test_refine_field import refine
        from pyphare.pharein import global_vars

        time_step_nbr = 1 # after first advance, simulations diverge

        # "un" = unrefined simulation - subject to change!
        un_diag_outputs=f"phare_outputs_subcycle_coarsening_un_{self.ddt_test_id()}"
        un_datahier = self.getHierarchy(interp_order, [], "eb", cells=30,
                                      diag_outputs=un_diag_outputs, time_step=0.001,
                                      extra_diag_options={"fine_dump_lvl_max": 10},
                                      time_step_nbr=time_step_nbr, smallest_patch_size=5,
                                      largest_patch_size=30)

        diag_outputs=f"phare_outputs_subcycle_coarsening_{self.ddt_test_id()}"
        datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", cells=30,
                                      diag_outputs=diag_outputs, time_step=0.001,
                                      extra_diag_options={"fine_dump_lvl_max": 10},
                                      time_step_nbr=time_step_nbr, smallest_patch_size=5,
                                      largest_patch_size=30)

        lvlSteps = global_vars.sim.level_time_steps
        assert len(lvlSteps) == 2  # this test is only configured for L0 -> L1 refinement comparisons

        def assert_time_in_hier(*ts):
            for t in ts:
                print(f"checking time {t} exists")
                self.assertIn(datahier.format_timestamp(t), datahier.times())

        coarsest_time_before = 0 # init
        coarsest_time_after = coarsest_time_before + lvlSteps[0]
        assert_time_in_hier(coarsest_time_before, coarsest_time_after)

        coarseUnBeforePatches = un_datahier.level(0, coarsest_time_before).patches
        coarseUnAfterPatches = un_datahier.level(0, coarsest_time_after).patches

        assert len(coarseUnBeforePatches) == len(coarseUnAfterPatches)

        for coarsePatch_idx in range(len(coarseUnBeforePatches)):
            coarseUnBeforePatch = coarseUnBeforePatches[coarsePatch_idx]
            coarseUnAfterPatch  = coarseUnAfterPatches[coarsePatch_idx]

            assert coarseUnBeforePatch.box == coarseUnAfterPatch.box
            coarsePatchBox      = coarseUnBeforePatch.box

            for finer_subcycle in range(0, 5):

                finer_subcycle_time = coarsest_time_before + (lvlSteps[1] * finer_subcycle)
                assert_time_in_hier(finer_subcycle_time)

                alpha = (finer_subcycle_time - coarsest_time_before) / (coarsest_time_after - coarsest_time_before)
                finer_subcycle_patches = datahier.level(1, finer_subcycle_time).patches

                for finer_subcycle_patch in finer_subcycle_patches: # L0 might be multiple patches if running with MPI
                    if boxm.coarsen(finer_subcycle_patch.box, refinement_ratio) in coarsePatchBox:

                        for EM in ["E", "B"]:
                            for xyz in ["x", "y", "z"]:
                                qty = f"{EM}{xyz}"
                                assert qty in coarseUnBeforePatch.patch_datas

                                coarseUnBefore_pd = coarseUnBeforePatch.patch_datas[qty]
                                coarseUnAfter_pd  = coarseUnAfterPatch.patch_datas[qty]
                                coarseUnBefore_dataset = coarseUnBefore_pd.dataset
                                coarseUnAfter_dataset = coarseUnAfter_pd.dataset

                                is_primal = coarseUnBefore_pd.primal_directions()

                                assert coarseUnBefore_dataset.shape == coarseUnAfter_dataset.shape
                                coarseUnInterp_dataset = coarseUnAfter_dataset[:].copy()

                                assert coarsePatchBox.ndim == 1 # update for > 1d interpolation values
                                if coarsePatchBox.ndim == 1:
                                    for ix in range(coarseUnInterp_dataset.shape[0]):
                                        coarseUnInterp_dataset[ix] = (1. - alpha) * coarseUnBefore_dataset[ix] + alpha * coarseUnAfter_dataset[ix]

                                refinedUnInterpolatedField = refine(coarseUnBefore_pd, data=coarseUnInterp_dataset)
                                finer_subcycle_pd = finer_subcycle_patch.patch_datas[qty]
                                ghosts = finer_subcycle_pd.ghosts_nbr
                                finer_subcycle_dataset = finer_subcycle_pd.dataset[:]

                                finer_ghost_boxes = finer_subcycle_pd.ghost_box - finer_subcycle_patch.box
                                for finer_ghost_box in finer_ghost_boxes:

                                    is_upper = np.array(finer_ghost_box.lower > finer_subcycle_patch.box.lower, dtype=int)
                                    finer_ghost_box = boxm.shift(finer_ghost_box, is_upper and is_primal)
                                    lvlOverlap = refinedUnInterpolatedField.box * finer_ghost_box

                                    if lvlOverlap is not None:
                                        finer_ghostbox_data = finer_subcycle_pd[finer_ghost_box]
                                        if ndim == 1:
                                            if is_upper[0]:
                                                np.testing.assert_allclose(finer_ghostbox_data, finer_subcycle_pd.dataset[-finer_subcycle_pd.ghosts_nbr[0]:], atol=1e-6)
                                            else:
                                                np.testing.assert_allclose(finer_ghostbox_data, finer_subcycle_pd.dataset[:finer_subcycle_pd.ghosts_nbr[0]], atol=1e-6)

                                        refinedInterp_data_box = Box(finer_ghost_box.lower[0] + ghosts[0], finer_ghost_box.upper[0] + ghosts[0])
                                        refinedUnInterpGhostBox_data = refinedUnInterpolatedField[refinedInterp_data_box]
                                        np.testing.assert_allclose(finer_ghostbox_data, refinedUnInterpGhostBox_data, atol=1e-6)





    @data( # only supports a hierarchy with 2 levels
       ({"L0": [Box1D(5, 9)]}),
       ({"L0": [Box1D(5, 9), Box1D(20, 24)]}),
    )
    def test_field_level_ghosts_via_subcycles_and_coarser_interpolation(self, refinement_boxes):
        dim = refinement_boxes["L0"][0].ndim
        for interp in [1, 2, 3]:
            self._test_field_level_ghosts_via_subcycles_and_coarser_interpolation(dim, interp_order=interp, refinement_boxes=refinement_boxes)


    @data(
       ({"L0": {"B0": Box1D(10, 14), "B1": Box1D(15, 19)}}),
    )
    def test_hierarchy_timestamp_cadence(self, refinement_boxes):
        dim = refinement_boxes["L0"]["B0"].ndim

        time_step     = .001
        # time_step_nbr chosen to force diagnostics dumping double imprecision cadence calculations accuracy testing
        time_step_nbr = 101
        final_time    = time_step * time_step_nbr

        for trailing in [0, 1]: # 1 = skip init dumps
            for i in [2, 3]:
                timestamps = np.arange(0, final_time, time_step*i)[trailing:]

                diag_outputs=f"phare_outputs_hierarchy_timestamp_cadence_{self.ddt_test_id()}_{i}"
                hier = self.getHierarchy(interp_order=1, refinement_boxes=refinement_boxes, qty="eb", cells=30,
                                              diag_outputs=diag_outputs, time_step=time_step,
                                              time_step_nbr=time_step_nbr, smallest_patch_size=5,
                                              largest_patch_size=30, timestamps=timestamps)

                time_hier_keys = list(hier.time_hier.keys())
                self.assertEqual(len(time_hier_keys), len(timestamps))

                for i, timestamp in enumerate(time_hier_keys):
                    self.assertEqual(hier.format_timestamp(timestamps[i]), timestamp)


if __name__ == "__main__":
    unittest.main()
