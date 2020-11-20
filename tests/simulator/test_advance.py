
from pybindlibs import cpp

from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein.diagnostics import ParticleDiagnostics, FluidDiagnostics, ElectromagDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.pharein.simulation import Simulation
from pyphare.pharesee.geometry import level_ghost_boxes, hierarchy_overlaps
from pyphare.pharesee.particles import aggregate as aggregate_particles
import pyphare.core.box as boxm
from pyphare.core.box import Box, Box1D
import numpy as np
import unittest
from ddt import ddt, data, unpack


@ddt
class AdvanceTest(unittest.TestCase):

    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]

    def getHierarchy(self, interp_order, refinement_boxes, qty, nbr_part_per_cell=100,
                     diag_outputs="phare_outputs",
                     smallest_patch_size=5, largest_patch_size=5,
                     cells= 120,
                     dl=0.1, extra_diag_options={}, advances=1):

        from pyphare.pharein import global_vars
        global_vars.sim = None
        startMPI()
        extra_diag_options["mode"] = "overwrite"
        extra_diag_options["dir"] = diag_outputs
        Simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=largest_patch_size,
            time_step_nbr=30000,
            final_time=30.,
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
            return 0.

        def by(x):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()[0]
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

        for quantity in ["E", "B"]:
            ElectromagDiagnostics(
                quantity=quantity,
                write_timestamps=np.zeros(advances+1),
                compute_timestamps=np.zeros(advances+1)
            )

        for quantity in ["density", "bulkVelocity"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=np.zeros(advances+1),
                compute_timestamps=np.zeros(advances+1)
            )

        for pop in ["protons"]:
            for quantity in ["density", "flux"]:
                FluidDiagnostics(quantity=quantity,
                                 write_timestamps=np.zeros(advances+1),
                                 compute_timestamps=np.zeros(advances+1),
                                 population_name=pop)

        simulator = Simulator(global_vars.sim).initialize()
        for i in range(advances):
            simulator.advance()

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

        advances = 3
        stepDiff = .1
        nSubcycles = 10
        coarsestTimeStep = 0.001
        levelNumbers = [i for i in range(len(refinement_boxes.keys()) + 1)]

        lvlSteps = [coarsestTimeStep * (stepDiff ** (ilvl)) for ilvl in levelNumbers]
        finestTimeStep = lvlSteps[-1]
        secondFinestTimeStep = lvlSteps[-2] # this test makes no sense with only 1 level
        totalSteps = nSubcycles ** levelNumbers[-1] * advances
        uniqTimes = set([0])

        diag_outputs=f"phare_outputs_subcycle_coarsening_{self.ddt_test_id()}"
        datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", cells=30,
                                      diag_outputs=diag_outputs,
                                      extra_diag_options={"fine_dump_lvl_max": 10},
                                      advances=advances, smallest_patch_size=5, largest_patch_size=30)

        times = datahier.times()
        for step in range(1, totalSteps + 1):
            checkTime = float("{:.6f}".format(finestTimeStep * step))
            self.assertIn(checkTime, times)
            uniqTimes.add(checkTime)

        self.assertEqual(len(uniqTimes), len(datahier.time_hier.items()))

        syncSteps = nSubcycles ** levelNumbers[-2] * advances # ignore most fine subcycles

        # FIX THIS AFTER NO MORE REGRIDS
        #  SEE: https://github.com/PHAREHUB/PHARE/issues/400
        assert syncSteps % advances == 0 # perfect division
        startStep = int(syncSteps / advances) + 1 # skip coarsest step due to issue 400
        for step in range(startStep, syncSteps + 1):
            checkTime = float("{:.6f}".format(secondFinestTimeStep * step))
            self.assertIn(checkTime, datahier.times())
            nLevels = len(datahier.time_hier[checkTime].items())
            self.assertGreaterEqual(nLevels, 2)
            levelNbrs = list(datahier.time_hier[checkTime].keys())
            finestLevelNbr = max(levelNbrs)
            coarsestLevelNbr = min(levelNbrs)

            for coarseLevelNbr in range(coarsestLevelNbr, finestLevelNbr):
                coarsePatches = datahier.time_hier[checkTime][coarseLevelNbr].patches
                finePatches = datahier.time_hier[checkTime][coarseLevelNbr + 1].patches

                for coarsePatch in coarsePatches:
                    for finePatch in finePatches:
                        lvlOverlap = boxm.refine(coarsePatch.box, 2) * finePatch.box
                        if lvlOverlap is not None:
                            for EM in ["E", "B"]:
                                for xyz in ["x", "y", "z"]:
                                    qty = f"EM_{EM}_{xyz}"
                                    coarse_pd = coarsePatch.patch_datas[qty]
                                    fine_pd  = finePatch.patch_datas[qty]
                                    coarseBox = boxm.coarsen(lvlOverlap, 2)

                                    qty = f"{EM}{xyz}"
                                    nGhosts = coarse_pd.layout.nbrGhostFor(qty)

                                    coarse_pdDataset = coarse_pd.dataset[:]
                                    fine_pdDataset = fine_pd.dataset[:]

                                    coarseOffset = coarseBox.lower - coarse_pd.layout.box.lower
                                    dataBox_lower = coarseOffset + nGhosts
                                    dataBox = Box(dataBox_lower, dataBox_lower + coarseBox.shape() - 1)

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
        dim = refinement_boxes["L0"]["B0"].dim()
        self._test_field_coarsening_via_subcycles(dim, interp_order=1, refinement_boxes=refinement_boxes)


if __name__ == "__main__":
    unittest.main()
