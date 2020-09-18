
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein.diagnostics import ParticleDiagnostics, FluidDiagnostics, ElectromagDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.pharein.simulation import Simulation
from pyphare.pharesee.geometry import level_ghost_boxes, hierarchy_overlaps, touch_domain_border
from pyphare.pharesee.particles import aggregate as aggregate_particles

import pyphare.core.box as boxm
from pyphare.core.box import Box, Box1D, Box2D, Box3D
import numpy as np
import unittest
from ddt import ddt, data, unpack



@ddt
class InitializationTest(unittest.TestCase):

    def getHierarchy(self, interp_order, refinement_boxes, qty, nbr_part_per_cell=100,
                     diag_outputs="phare_outputs",
                     density = lambda x: 0.3 + 1./np.cosh((x-6)/4.)**2,
                     beam = False,
                     smallest_patch_size=10, largest_patch_size=10,
                     cells= 120,
                     dl=0.1):

        from pyphare.pharein import global_vars
        global_vars.sim =None

        Simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=largest_patch_size,
            time_step_nbr=30000,
            final_time=30.,
            boundary_types="periodic",
            cells=cells,
            dl=dl,
            interp_order=interp_order,
            max_nbr_levels=len(refinement_boxes)+1 if refinement_boxes is not None else 1,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5",
                          "options": {"dir": diag_outputs, "mode":"overwrite"}}
        )
        def beam_density(x):
            return np.zeros_like(x)+0.3

        def by(x):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            return 0.1*np.cos(2*np.pi*x/L[0])

        def bz(x):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            return 0.1*np.sin(2*np.pi*x/L[0])

        def bx(x):
            return 1.

        def vx(x):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            return 0.1*np.cos(2*np.pi*x/L[0]) + 0.2

        def vy(x):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            return 0.1*np.cos(2*np.pi*x/L[0])

        def vz(x):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            return 0.1*np.sin(2*np.pi*x/L[0])

        def vthx(x):
            return 0.01 + np.zeros_like(x)

        def vthy(x):
            return 0.01 + np.zeros_like(x)

        def vthz(x):
            return 0.01 + np.zeros_like(x)

        if beam:
            MaxwellianFluidModel(bx=bx, by=by, bz=bz,
                                 protons={"charge": 1,
                                          "density": density,
                                          "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
                                          "vthx": vthx, "vthy": vthy, "vthz": vthz,
                                          "nbr_part_per_cell": nbr_part_per_cell,
                                          "init": {"seed": 1337}},

                                 beam={"charge": 1,
                                       "density": beam_density,
                                       "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
                                       "vthx": vthx, "vthy": vthy, "vthz": vthz,
                                       "nbr_part_per_cell": nbr_part_per_cell,
                                       "init": {"seed": 1337}})

        else:
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
                write_timestamps=np.zeros(1),
                compute_timestamps=np.zeros(1)
            )



        for quantity in ["density", "bulkVelocity"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=np.zeros(1),
                compute_timestamps=np.zeros(1)
            )

        poplist = ["protons", "beam"] if beam else ["protons"]
        for pop in poplist:
            for quantity in ["density", "flux"]:
                FluidDiagnostics(quantity=quantity,
                                 write_timestamps=np.zeros(1),
                                 compute_timestamps=np.zeros(1),
                                 population_name=pop)

            for quantity in ['domain', 'levelGhost', 'patchGhost']:
                ParticleDiagnostics(quantity=quantity,
                                    compute_timestamps=np.zeros(1),
                                    write_timestamps=np.zeros(1),
                                    population_name=pop)

        simulator = Simulator(global_vars.sim)
        simulator.initialize()

        if qty == "b":
            b_hier = hierarchy_from(h5_filename=diag_outputs+"/EM_B.h5")
            return b_hier


        if qty == "particles":
            particle_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_domain.h5")
            particle_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_patchGhost.h5", hier=particle_hier)
            particle_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_levelGhost.h5", hier=particle_hier)
            merge_particles(particle_hier)
            return particle_hier

        if qty == "moments":
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_density.h5")
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_bulkVelocity.h5", hier=mom_hier)
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_density.h5", hier=mom_hier)
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_flux.h5", hier=mom_hier)
            if beam:
                mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_beam_density.h5", hier=mom_hier)
                mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_beam_flux.h5", hier=mom_hier)
            return mom_hier





    @data(1,2,3)
    def test_B_is_as_provided_by_user(self, interp_order):
        print("test_B_is_as_provided_by_user : interp_order : {}".format(interp_order))
        hier = self.getHierarchy(interp_order, {"L0": {"B0": [(10, ), (20, )]}}, "b")

        from pyphare.pharein import global_vars
        model = global_vars.sim.model
        bx_fn = model.model_dict["bx"]
        by_fn = model.model_dict["by"]
        bz_fn = model.model_dict["bz"]

        for ilvl, level in hier.levels().items():
            print("checking level {}".format(ilvl))
            for ip, patch in enumerate(level.patches):

                xbx   = patch.patch_datas["EM_B_x"].x[:]
                bx  = patch.patch_datas["EM_B_x"].dataset[:]
                np.testing.assert_allclose(bx, bx_fn(xbx), atol=3e-5)

                xby   = patch.patch_datas["EM_B_y"].x[:]
                by  = patch.patch_datas["EM_B_y"].dataset[:]
                np.testing.assert_allclose(by, by_fn(xby), atol=3e-5)

                xbz   = patch.patch_datas["EM_B_z"].x[:]
                bz  = patch.patch_datas["EM_B_z"].dataset[:]
                np.testing.assert_allclose(bz, bz_fn(xbz), atol=3e-5)




    @data((1, {"L0": {"B0": [(10, ), (20, )]}}),
          (2, {"L0": {"B0": [(10, ), (20, )]}}),
          (3, {"L0": {"B0": [(10, ), (20, )]}}),
          (1, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}),
          (2, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}),
          (3, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}))
    @unpack
    def test_overlaped_fields_are_equal(self, interp_order, refinement_boxes):
        print("test_overlaped_fields_are_equal")
        hier = self.getHierarchy(interp_order, refinement_boxes, "b")

        overlaps = hierarchy_overlaps(hier)
        check=0
        for ilvl, lvl in hier.levels().items():

            for overlap in overlaps[ilvl]:
                pd1, pd2 = overlap["pdatas"]
                box      = overlap["box"]
                offsets  = overlap["offset"]

                self.assertEqual(pd1.quantity, pd2.quantity)

                if pd1.quantity == 'field':
                    check+=1

                    # we need to transform the AMR overlap box, which is thus
                    # (because AMR) common to both pd1 and pd2 into local index
                    # boxes that will allow to slice the data

                    # the patchData ghost box that serves as a reference box
                    # to transfrom AMR to local indexes first needs to be
                    # shifted by the overlap offset associated to it
                    # this is because the overlap box has been calculated from
                    # the intersection of possibly shifted patch data ghost boxes

                    loc_b1 = boxm.amr_to_local(box, boxm.shift(pd1.ghost_box, offsets[0]))
                    loc_b2 = boxm.amr_to_local(box, boxm.shift(pd2.ghost_box, offsets[1]))

                    data1 = pd1.dataset
                    data2 = pd2.dataset

                    slice1 = data1[loc_b1.lower[0]:loc_b1.upper[0] + 1]
                    slice2 = data2[loc_b2.lower[0]:loc_b2.upper[0] + 1]

                    self.assertTrue(np.allclose(slice1, slice2, atol=1e-12))

        self.assertTrue(check>0)




    @data(1, 2, 3)
    def test_bulkvel_is_as_provided_by_user(self, interp_order):
        print("test_density_is_as_provided_by_user : interp_order : {}".format(interp_order))
        hier = self.getHierarchy(interp_order,
                                 {"L0": {"B0": [(10, ), (20, )]}},
                                 "moments",
                                 nbr_part_per_cell=10000, beam=True)

        from pyphare.pharein import global_vars
        model = global_vars.sim.model
        # protons and beam have same bulk vel here so take
        # only proton func.
        vx_fn = model.model_dict["protons"]["vx"]
        vy_fn = model.model_dict["protons"]["vy"]
        vz_fn = model.model_dict["protons"]["vz"]
        nprot = model.model_dict["protons"]["density"]
        nbeam = model.model_dict["beam"]["density"]


        for ilvl, level in hier.levels().items():
            print("checking density on level {}".format(ilvl))
            for ip,patch in enumerate(level.patches):
                print("patch {}".format(ip))

                layout    = patch.patch_datas["protons_flux_x"].layout
                centering = layout.centering["X"][patch.patch_datas["protons_flux_x"].field_name]
                nbrGhosts = layout.nbrGhosts(interp_order, centering)

                x   = patch.patch_datas["protons_flux_x"].x[nbrGhosts:-nbrGhosts]
                fpx = patch.patch_datas["protons_flux_x"].dataset[nbrGhosts:-nbrGhosts]
                fpy = patch.patch_datas["protons_flux_y"].dataset[nbrGhosts:-nbrGhosts]
                fpz = patch.patch_datas["protons_flux_z"].dataset[nbrGhosts:-nbrGhosts]
                fbx = patch.patch_datas["protons_flux_x"].dataset[nbrGhosts:-nbrGhosts]
                fby = patch.patch_datas["protons_flux_y"].dataset[nbrGhosts:-nbrGhosts]
                fbz = patch.patch_datas["protons_flux_z"].dataset[nbrGhosts:-nbrGhosts]
                ni  = patch.patch_datas["density"].dataset[nbrGhosts:-nbrGhosts]

                vxact = (fpx + fbx)/ni
                vyact = (fpy + fby)/ni
                vzact = (fpz + fbz)/ni

                vxexp =(nprot(x) * vx_fn(x) + nbeam(x) * vx_fn(x))/(nprot(x)+nbeam(x))
                vyexp =(nprot(x) * vy_fn(x) + nbeam(x) * vy_fn(x))/(nprot(x)+nbeam(x))
                vzexp =(nprot(x) * vz_fn(x) + nbeam(x) * vz_fn(x))/(nprot(x)+nbeam(x))

                for vexp, vact in zip((vxexp, vyexp, vzexp), (vxact, vyact, vzact)):
                    std = np.std(vexp-vact)
                    print("sigma(user v - actual v) = {}".format(std))
                    self.assertTrue(std < 1e-2) # empirical value obtained from print just above





    @data(1, 2, 3)
    def test_density_is_as_provided_by_user(self, interp_order):
        print("test_density_is_as_provided_by_user : interp_order : {}".format(interp_order))
        hier = self.getHierarchy(interp_order,
                                 {"L0": {"B0": [(10, ), (20, )]}},
                                 "moments",
                                 nbr_part_per_cell=10000, beam=True)

        from pyphare.pharein import global_vars
        model = global_vars.sim.model
        proton_density_fn = model.model_dict["protons"]["density"]
        beam_density_fn = model.model_dict["beam"]["density"]

        for ilvl, level in hier.levels().items():
            print("checking density on level {}".format(ilvl))
            for ip,patch in enumerate(level.patches):
                print("patch {}".format(ip))

                ion_density     = patch.patch_datas["density"].dataset[:]
                proton_density  = patch.patch_datas["protons_density"].dataset[:]
                beam_density    = patch.patch_datas["beam_density"].dataset[:]
                x               = patch.patch_datas["density"].x

                layout    = patch.patch_datas["density"].layout
                centering = layout.centering["X"][patch.patch_datas["density"].field_name]
                nbrGhosts = layout.nbrGhosts(interp_order, centering)

                protons_expected = proton_density_fn(x[nbrGhosts:-nbrGhosts])
                beam_expected    = beam_density_fn(x[nbrGhosts:-nbrGhosts])
                ion_expected     = protons_expected + beam_expected

                ion_actual     = ion_density[nbrGhosts:-nbrGhosts]
                beam_actual    = beam_density[nbrGhosts:-nbrGhosts]
                protons_actual = proton_density[nbrGhosts:-nbrGhosts]

                names    = ("ions", "protons", "beam")
                expected = (ion_expected, protons_expected, beam_expected)
                actual   = (ion_actual, protons_actual, beam_actual)
                devs = {name:np.std(expected-actual) for name, expected, actual in zip(names, expected, actual)}

                for name, dev in devs.items():
                    print("sigma(user density - {} density) = {}".format(name, dev))

                for name,dev in devs.items():
                    self.assertTrue(dev < 6e-3, '{} has dev = {}'.format(name, dev))  # empirical value obtained from test prints





    @data(1, 2, 3)
    def test_density_decreases_as_1overSqrtN(self,interp_order):
        import matplotlib
        matplotlib.use("Agg")  # for systems without GUI
        import matplotlib.pyplot as plt
        print("test_density_decreases_as_1overSqrtN, interp_order = {}".format(interp_order))

        nbr_particles = np.asarray([100, 1000, 5000, 10000])
        noise = np.zeros(len(nbr_particles))

        for inbr,nbrpart in enumerate(nbr_particles):

            hier = self.getHierarchy(interp_order, None, "moments",
                                     nbr_part_per_cell=nbrpart,
                                     diag_outputs="phare_outputs_{}".format(nbrpart),
                                     density=lambda x:np.zeros_like(x)+1.,
                                     smallest_patch_size=480,
                                     largest_patch_size=480,
                                     cells=960,
                                     dl=0.0125)

            from pyphare.pharein import global_vars
            model   = global_vars.sim.model
            protons = model.model_dict["protons"]
            density_fn = protons["density"]

            patch       =  hier.level(0).patches[0]
            ion_density = patch.patch_datas["density"].dataset[:]
            x           = patch.patch_datas["density"].x

            layout = patch.patch_datas["density"].layout
            centering = layout.centering["X"][patch.patch_datas["density"].field_name]
            nbrGhosts = layout.nbrGhosts(interp_order, centering)

            expected = density_fn(x[nbrGhosts:-nbrGhosts])
            actual  = ion_density[nbrGhosts:-nbrGhosts]
            noise[inbr] = np.std(expected-actual)
            print("noise is {} for {} particles per cell".format(noise[inbr], nbrpart))

            plt.figure()
            plt.plot(x[nbrGhosts:-nbrGhosts], actual, label="actual")
            plt.plot(x[nbrGhosts:-nbrGhosts], expected, label="expected")
            plt.legend()
            plt.title(r"$\sigma =$ {}".format(noise[inbr]))
            plt.savefig("noise_{}_interp_{}.png".format(nbrpart, interp_order))
            plt.close("all")

        plt.figure()
        plt.plot(nbr_particles, noise/noise[0], label=r"$\sigma/\sigma_0$")
        plt.plot(nbr_particles, 1/np.sqrt(nbr_particles/nbr_particles[0]), label=r"$1/sqrt(nppc/nppc0)$")
        plt.xlabel("nbr_particles")
        plt.legend()
        plt.savefig("noise_nppc_interp_{}.png".format(interp_order))
        plt.close("all")

        noiseMinusTheory = noise/noise[0] - 1/np.sqrt(nbr_particles/nbr_particles[0])
        plt.figure()
        plt.plot(nbr_particles, noiseMinusTheory,
                 label=r"$\sigma/\sigma_0 - 1/sqrt(nppc/nppc0)$")
        plt.xlabel("nbr_particles")
        plt.legend()
        plt.savefig("noise_nppc_minus_theory_interp_{}.png".format(interp_order))
        plt.close("all")
        self.assertTrue(noiseMinusTheory[1:].mean() < 2e-2)




    @data(1, 2, 3)
    def test_nbr_particles_per_cell_is_as_provided(self, interp_order):
        datahier = self.getHierarchy(interp_order, {"L0": {"B0": [(10, ), (20, )]}}, "particles")
        print("test_nbr_particles_per_cell_is_as_provided, interp_order = {}".format(interp_order))
        L0 = datahier.level(0)
        for patch in L0.patches:
            pd = patch.patch_datas["protons_particles"]
            icells = pd.dataset.iCells
            mincell = icells.min()
            # bincount only works for non-negative values
            # but icells could be -1 or -2 for interp order 1 or (2,3)
            # so we artificially add the min (-1 or -2) and count the
            # number of occurence of cell indexes
            # this should be a list of only nbr_part_per_cell
            counts = np.bincount(icells-mincell)
            self.assertTrue(np.all(counts == 100)) #100 is default nbr for maxwellian model.





    @data((1, {"L0": {"B0": [(10, ), (20, )]}}),
          (2, {"L0": {"B0": [(10, ), (20, )]}}),
          (3, {"L0": {"B0": [(10, ), (20, )]}}),
          (1, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}),
          (2, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}),
          (3, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}))
    @unpack
    def test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(self,interp_order, refinement_boxes):
        datahier = self.getHierarchy(interp_order, refinement_boxes, "particles")
        print("test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles")
        print("interporder : {}".format(interp_order))
        overlaps = hierarchy_overlaps(datahier)
        for ilvl, lvl_overlaps in overlaps.items():
            print("level {}".format(ilvl))
            for overlap in lvl_overlaps:

                if ilvl != 0: #only root level tested here
                    continue

                if "particles"  not in overlap["pdatas"][0].quantity :
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



    @data((1, {"L0": {"B0": [(10, ), (20, )]}}),
          (2, {"L0": {"B0": [(10, ), (20, )]}}),
          (3, {"L0": {"B0": [(10, ), (20, )]}}),
          (1, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}),
          (2, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}),
          (3, {"L0": {"B0": [(2, ), (12, )], "B1": [(13,), (25,)]}}))
    @unpack
    def test_overlapped_particledatas_have_identical_particles(self, interp_order, refinement_boxes):
        print("test_overlapped_particledatas_have_identical_particles")
        from copy import copy
        datahier = self.getHierarchy(interp_order, refinement_boxes, "particles")
        print("interporder : {}".format(interp_order))
        overlaps = hierarchy_overlaps(datahier)

        for ilvl, lvl in datahier.patch_levels.items():

            print("testing level {}".format(ilvl))
            for overlap in overlaps[ilvl]:

                pd1, pd2 = overlap["pdatas"]
                box      = overlap["box"]
                offsets  = overlap["offset"]

                self.assertEqual(pd1.quantity, pd2.quantity)

                if "particles" in pd1.quantity:

                    # the following uses 'offset', we need to remember that offset
                    # is the quantity by which a patch has been moved to detect
                    # overlap with the other one.
                    # so shift by +offset when evaluating patch data in overlap box
                    # index space, and by -offset when we want to shift box indexes
                    # to the associated patch index space.

                    # overlap box must be shifted by -offset to select data in the patches
                    part1 = copy(pd1.dataset.select(boxm.shift(box, -offsets[0])))
                    part2 = copy(pd2.dataset.select(boxm.shift(box, -offsets[1])))

                    idx1 = np.argsort(part1.iCells + part1.deltas)
                    idx2 = np.argsort(part2.iCells + part2.deltas)

                    # if there is an overlap, there should be particles
                    # in these cells
                    assert(len(idx1) >0)
                    assert(len(idx2) >0)

                    print("respectively {} and {} in overlaped patchdatas".format(len(idx1), len(idx2)))

                    # particle iCells are in their patch AMR space
                    # so we need to shift them by +offset to move them to the box space
                    np.testing.assert_array_equal(part1.iCells[idx1]+offsets[0], part2.iCells[idx2]+offsets[1])

                    self.assertTrue(np.allclose(part1.deltas[idx1], part2.deltas[idx2], atol=1e-12))
                    self.assertTrue(np.allclose(part1.v[idx1,0], part2.v[idx2,0], atol=1e-12))
                    self.assertTrue(np.allclose(part1.v[idx1,1], part2.v[idx2,1], atol=1e-12))
                    self.assertTrue(np.allclose(part1.v[idx1,2], part2.v[idx2,2], atol=1e-12))



    def _test_levelghostparticles_have_correct_split_from_coarser_particle(self, dim, interp_order, refinement_boxes):
        print("test_levelghostparticles_have_correct_split_from_coarser_particle for dim/interp : {}/{}".format(dim, interp_order))

        datahier = self.getHierarchy(interp_order, refinement_boxes, "particles")
        from pyphare.pharein.global_vars import sim
        assert sim is not None

        def domainParticles_for(ilvl, particles_id):
            patchDatas = [
              patch.patch_datas[particles_id]
                  for patch in datahier.levels()[ilvl].patches
                      if particles_id in patch.patch_datas
            ]
            return aggregate_particles([
              patchData.dataset.select(patchData.box) for patchData in patchDatas
            ]) # including patch ghost particles means duplicates

        for ilvl, particle_gaboxes in level_ghost_boxes(datahier).items():
            assert ilvl > 0 # has no level 0

            for particles_id, gaboxes_list in particle_gaboxes.items():
                coarse_particles = domainParticles_for(ilvl - 1, particles_id)
                coarse_split_particles = coarse_particles.split(sim)
                assert coarse_split_particles.size() == coarse_particles.size() * sim.refined_particle_nbr

                for gabox in gaboxes_list:
                    gabox_patchData = gabox["pdata"]

                    for ghostBox in gabox["boxes"]:
                        part1 = gabox_patchData.dataset.select(ghostBox)
                        part2 = coarse_split_particles.select(ghostBox)

                        idx1 = np.argsort(part1.iCells + part1.deltas)
                        idx2 = np.argsort(part2.iCells + part2.deltas)

                        np.testing.assert_array_equal(part1.iCells[idx1], part2.iCells[idx2])

                        self.assertTrue(np.allclose(part1.deltas[idx1], part2.deltas[idx2], atol=1e-12))
                        self.assertTrue(np.allclose(part1.v[idx1,0], part2.v[idx2,0], atol=1e-12))
                        self.assertTrue(np.allclose(part1.v[idx1,1], part2.v[idx2,1], atol=1e-12))
                        self.assertTrue(np.allclose(part1.v[idx1,2], part2.v[idx2,2], atol=1e-12))


    @data(
       ({"L0": {"B0": Box1D(10, 20)}}),
       ({"L0": {"B0": Box1D(10, 100)}, "L1": {"B0": Box1D(40, 60)}}),
       ({"L0": {"B0": Box1D( 2, 12), "B1": Box1D(13, 25)}}),
    )
    def test_levelghostparticles_have_correct_split_from_coarser_particle(self, refinement_boxes):
        dim = len(refinement_boxes["L0"]["B0"].lower)
        for interp_order in [1, 2, 3]:
            self._test_levelghostparticles_have_correct_split_from_coarser_particle(dim, interp_order, refinement_boxes)

if __name__ == "__main__":
    unittest.main()
