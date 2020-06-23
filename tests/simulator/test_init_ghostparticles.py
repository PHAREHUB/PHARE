
from pybindlibs import cpp
from pyphare.simulator.simulator import create_simulator  # keep that guy first
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein import ElectronModel
from pyphare.pharein.simulation import Simulation
from pyphare.pharesee.geometry import level_ghost_boxes, hierarchy_overlaps, touch_domain_border
import pyphare.core.box as boxm
from pyphare.core.box import Box
import numpy as np
import unittest
from ddt import ddt, data, unpack


@ddt
class ParticleInitializationTest(unittest.TestCase):

    def getHierarchy(self, interp_order, refinement_boxes):

        from pyphare.pharein import globals
        globals.sim =None

        Simulation(
            smallest_patch_size=10,
            largest_patch_size=10,
            time_step_nbr=30000,
            final_time=30.,
            boundary_types="periodic",
            cells=40,
            dl=0.3,
            interp_order=interp_order,
            max_nbr_levels=len(refinement_boxes)+1,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5",
                          "options": {"dir": "phare_outputs"}}
        )

        def density(x):
            return 1.

        def by(x):
            from pyphare.pharein.globals import sim
            L = sim.simulation_domain()
            return 0.1*np.cos(2*np.pi*x/L[0])

        def bz(x):
            from pyphare.pharein.globals import sim
            L = sim.simulation_domain()
            return 0.1*np.sin(2*np.pi*x/L[0])

        def bx(x):
            return 1.

        def vx(x):
            return 0.

        def vy(x):
            from pyphare.pharein.globals import sim
            L = sim.simulation_domain()
            return 0.1*np.cos(2*np.pi*x/L[0])

        def vz(x):
            from pyphare.pharein.globals import sim
            L = sim.simulation_domain()
            return 0.1*np.sin(2*np.pi*x/L[0])

        def vthx(x):
            return 0.01

        def vthy(x):
            return 0.01

        def vthz(x):
            return 0.01

        MaxwellianFluidModel(bx=bx, by=by, bz=bz,
                             protons={"charge": 1,
                                      "density": density,
                                      "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
                                      "vthx": vthx, "vthy": vthy, "vthz": vthz,
                                      "init": {"seed": 1337}})

        ElectronModel(closure="isothermal", Te=0.12)

        dman, simulator, hier = create_simulator()
        simulator.initialize()

        datahier = hierarchy_from(simulator, hier, "particles", pop="protons")
        return dman,simulator, hier, datahier


    @data((1, {"L0": {"B0": [(10, ), (20, )]}}),
          (2, {"L0": {"B0": [(10, ), (20, )]}}),
          (3, {"L0": {"B0": [(10, ), (20, )]}}),
          (1, {"L0": {"B0": [(2, ), (12, )], "B1":[(15,),(25,)]}}),
          (2, {"L0": {"B0": [(2, ), (12, )], "B1":[(15,),(25,)]}}),
          (3, {"L0": {"B0": [(2, ), (12, )], "B1":[(15,),(25,)]}}))
    @unpack
    def test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(self,interp_order, refinement_boxes):

        dman, simulator, hier, datahier = self.getHierarchy(interp_order, refinement_boxes)
        print("interporder : {}".format(interp_order))
        overlaps = hierarchy_overlaps(datahier)
        for ilvl, lvl_overlaps in overlaps.items():
            for overlap in lvl_overlaps:

                if ilvl != 0: #only root level tested here
                    continue

                if overlap["pdatas"][0].quantity != "particles":
                    continue

                ref_pd, cmp_pd = overlap["pdatas"]

                box = overlap["box"]
                print("overlap box : {}, reference patchdata box : {}, ghostbox {},"
                " comp. patchdata box : {} ghostbox {}".format(box,ref_pd.box,ref_pd.ghost_box,cmp_pd.box, cmp_pd.ghost_box))
                offsets = overlap["offset"]

                # first let's shift the overlap box over the AMR
                # indices of the patchdata. The box has been created
                # by shifting the patchdata ghost box by 'offset' so here
                # the box is shifted by -offset to get over patchdata
                shift_refbox, shift_cmpbox = [boxm.shift(box, -off) for off in offsets]

                # the overlap box overlaps both ghost and domain cells
                # we need to extract the domain ones to later select domain
                # particles
                ovlped_refdom = ref_pd.box * shift_refbox
                ovlped_cmpdom = cmp_pd.box * shift_cmpbox

                # on lvl 0 patches are adjacent
                # therefore the overlap box must overlap the
                # patchData box. 1 cell in interporder1, 2 cells for higher
                #self.assertIsNotNone(ovlped_refdom)
                #self.assertIsNotNone(ovlped_cmpdom)
                assert(ovlped_cmpdom is not None)
                assert(ovlped_refdom is not None)

                refdomain = ref_pd.dataset.select(ovlped_refdom)
                cmpdomain = cmp_pd.dataset.select(ovlped_cmpdom)

                # now get the ghost cells of each patch data overlaped by
                # the overlap box. To do this we need to intersect the shifted
                # overlap box with the patchdata ghost box, and remove interior cells
                # note that in 1D we don't expect remove to return more than 1 box, hence [0]
                ovlped_refghost = boxm.remove(ref_pd.ghost_box * shift_refbox, ref_pd.box)[0]
                ovlped_cmpghost = boxm.remove(cmp_pd.ghost_box * shift_cmpbox, cmp_pd.box)[0]

                refghost  = ref_pd.dataset.select(ovlped_refghost)
                cmpghost  = cmp_pd.dataset.select(ovlped_cmpghost)

                print("ghost box {} has {} particles".format(ovlped_refghost, len(refghost.iCells)))
                print("ghost box {} has {} particles".format(ovlped_cmpghost, len(cmpghost.iCells)))

                # before comparing the particles we need to be sure particles of both patchdatas
                # are sorted in the same order. We do that by sorting by x position
                sort_refdomain_idx = np.argsort(refdomain.iCells + refdomain.deltas)
                sort_cmpdomain_idx = np.argsort(cmpdomain.iCells + cmpdomain.deltas)
                sort_refghost_idx = np.argsort(refghost.iCells + refghost.deltas)
                sort_cmpghost_idx = np.argsort(cmpghost.iCells + cmpghost.deltas)

                assert(sort_refdomain_idx.size != 0)
                assert(sort_cmpdomain_idx.size != 0)
                assert(sort_refdomain_idx.size != 0)
                assert(sort_cmpghost_idx.size != 0)

                np.testing.assert_allclose(refdomain.deltas[sort_refdomain_idx], cmpghost.deltas[sort_cmpghost_idx], atol=1e-12)
                np.testing.assert_allclose(cmpdomain.deltas[sort_cmpdomain_idx], refghost.deltas[sort_refghost_idx], atol=1e-12)


        del (dman, simulator, hier)
        cpp.reset()

if __name__ == "__main__":
    unittest.main()
