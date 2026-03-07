#


import os
import numpy as np


import pyphare.pharein as ph
from pyphare.core.box import nDBox
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.geometry import level_ghost_boxes
from pyphare.pharesee.hierarchy.hierarchy_utils import merge_particles
from pyphare.pharesee.particles import aggregate as aggregate_particles

from tests.simulator.test_initialization import InitializationTest


class HybridInitializationTest(InitializationTest):
    def getHierarchy(
        self,
        ndim,
        interp_order,
        qty,
        refinement_boxes,
        nbr_part_per_cell=100,
        density=None,
        extra_diag_options=None,
        beam=False,
        time_step_nbr=1,
        smallest_patch_size=None,
        largest_patch_size=10,
        cells=120,
        dl=0.1,
        diag_outputs="",
        **kwargs,
    ):
        if smallest_patch_size is None:
            from pyphare.pharein.simulation import check_patch_size

            _, smallest_patch_size = check_patch_size(
                ndim, interp_order=interp_order, cells=cells
            )

        base_diag_dir = "phare_outputs/init"
        base_diag_dir = (
            os.path.join(base_diag_dir, diag_outputs) if diag_outputs else base_diag_dir
        )
        extra_diag_options = extra_diag_options or dict()
        extra_diag_options["dir"] = base_diag_dir
        extra_diag_options["mode"] = "overwrite"
        sim = self.simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=largest_patch_size,
            time_step_nbr=time_step_nbr,
            final_time=30.0,
            boundary_types=["periodic"] * ndim,
            cells=[cells] * ndim,
            dl=[dl] * ndim,
            interp_order=interp_order,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5", "options": extra_diag_options},
            strict=True,
            **kwargs,
        )
        diag_outputs = sim.diag_options["options"]["dir"]
        L = sim.simulation_domain()

        def _density(*xyz):
            hL = np.array(sim.simulation_domain()) / 2
            return 0.3 + np.exp(
                sum([-((xyz[i] - hL[i]) ** 2) for i in range(len(xyz))])
            )

        def beam_density(*xyz):
            return np.zeros_like(xyz[0]) + 0.3

        def bx(*xyz):
            return 1.0

        def by(*xyz):
            return np.asarray(
                [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
            ).prod(axis=0)

        def bz(*xyz):
            return np.asarray(
                [0.1 * np.sin(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
            ).prod(axis=0)

        def vx(*xyz):
            return np.asarray(
                [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
            ).prod(axis=0)

        def vy(*xyz):
            return np.asarray(
                [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
            ).prod(axis=0)

        def vz(*xyz):
            return np.asarray(
                [0.1 * np.sin(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
            ).prod(axis=0)

        def vth(*xyz):
            return 0.01 + np.zeros_like(xyz[0])

        def vthx(*xyz):
            return vth(*xyz)

        def vthy(*xyz):
            return vth(*xyz)

        def vthz(*xyz):
            return vth(*xyz)

        protons = {
            "charge": 1,
            "density": density or _density,
            "vbulkx": vx,
            "vbulky": vy,
            "vbulkz": vz,
            "vthx": vthx,
            "vthy": vthy,
            "vthz": vthz,
            "nbr_part_per_cell": nbr_part_per_cell,
            "init": {"seed": 1337},
        }

        if beam:
            ph.MaxwellianFluidModel(
                bx=bx,
                by=by,
                bz=bz,
                protons=protons,
                beam={
                    "charge": 1,
                    "density": beam_density,
                    "vbulkx": vx,
                    "vbulky": vy,
                    "vbulkz": vz,
                    "vthx": vthx,
                    "vthy": vthy,
                    "vthz": vthz,
                    "nbr_part_per_cell": nbr_part_per_cell,
                    "init": {"seed": 1337},
                },
            )

        else:
            ph.MaxwellianFluidModel(bx=bx, by=by, bz=bz, protons=protons)

        ph.ElectronModel(closure="isothermal", Te=0.12)

        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(
                quantity=quantity, write_timestamps=np.zeros(time_step_nbr)
            )

        for quantity in ["charge_density", "bulkVelocity"]:
            ph.FluidDiagnostics(
                quantity=quantity, write_timestamps=np.zeros(time_step_nbr)
            )

        poplist = ["protons", "beam"] if beam else ["protons"]
        for pop in poplist:
            for quantity in ["density", "flux"]:
                ph.FluidDiagnostics(
                    quantity=quantity,
                    write_timestamps=np.zeros(time_step_nbr),
                    population_name=pop,
                )

            for quantity in ["domain", "levelGhost"]:
                ph.ParticleDiagnostics(
                    quantity=quantity,
                    write_timestamps=np.zeros(time_step_nbr),
                    population_name=pop,
                )

        Simulator(sim).initialize().reset()

        eb_hier = None
        if qty in ["e", "eb"]:
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_E.h5", hier=eb_hier
            )
        if qty in ["b", "eb"]:
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_B.h5", hier=eb_hier
            )
        if qty in ["e", "b", "eb"]:
            return eb_hier

        if qty == "moments":
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_charge_density.h5"
            )
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_bulkVelocity.h5", hier=mom_hier
            )
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_pop_protons_density.h5", hier=mom_hier
            )
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_pop_protons_flux.h5", hier=mom_hier
            )
            if beam:
                mom_hier = hierarchy_from(
                    h5_filename=diag_outputs + "/ions_pop_beam_density.h5",
                    hier=mom_hier,
                )
                mom_hier = hierarchy_from(
                    h5_filename=diag_outputs + "/ions_pop_beam_flux.h5", hier=mom_hier
                )
            return mom_hier

        # else particles

        assert qty == "particles"

        particle_hier = hierarchy_from(
            h5_filename=diag_outputs + "/ions_pop_protons_domain.h5"
        )
        particle_hier = hierarchy_from(
            h5_filename=diag_outputs + "/ions_pop_protons_levelGhost.h5",
            hier=particle_hier,
        )

        merge_particles(particle_hier)

        return particle_hier

    def _test_nbr_particles_per_cell_is_as_provided(
        self, dim, interp_order, default_ppc=100
    ):
        datahier = self.getHierarchy(
            dim,
            interp_order,
            "particles",
            {"L0": {"B0": nDBox(dim, 10, 20)}},
        )

        for patch in datahier.level(0).patches:
            pd = patch.patch_datas["protons_particles"]
            icells = pd.dataset[patch.box].iCells
            H, _ = np.histogramdd(icells)
            self.assertTrue((H == default_ppc).all())

    def _domainParticles_for(self, datahier, ilvl):
        patch0 = datahier.levels()[ilvl].patches[0]
        pop_names = [
            key for key in patch0.patch_datas.keys() if key.endswith("particles")
        ]
        particlePatchDatas = {k: [] for k in pop_names}
        for patch in datahier.levels()[ilvl].patches:
            for pop_name, patch_data in patch.patch_datas.items():
                particlePatchDatas[pop_name].append(patch_data)
        return {
            pop_name: aggregate_particles(
                [patchData.dataset.select(patchData.box) for patchData in patchDatas]
            )  # including patch ghost particles means duplicates
            for pop_name, patchDatas in particlePatchDatas.items()
        }

    def _test_domainparticles_have_correct_split_from_coarser_particle(
        self, ndim, interp_order, refinement_boxes, **kwargs
    ):
        print(
            "test_domainparticles_have_correct_split_from_coarser_particle for dim/interp : {}/{}".format(
                ndim, interp_order
            )
        )
        datahier = self.getHierarchy(
            ndim,
            interp_order,
            "particles",
            refinement_boxes,
            cells=30,
            **kwargs,
        )

        from pyphare.pharein.global_vars import sim

        assert sim is not None and len(sim.cells) == ndim

        levels = datahier.levels()
        self.assertTrue(len(levels) > 1)

        for ilvl in range(1, len(levels)):
            self.assertTrue(ilvl > 0)  # skip level 0
            level = levels[ilvl]
            coarse_particles = self._domainParticles_for(datahier, ilvl - 1)

            self.assertTrue(
                all([particles.size() > 0 for _, particles in coarse_particles.items()])
            )

            coarse_split_particles = {
                k: particles.split(sim) for k, particles in coarse_particles.items()
            }

            for k, particles in coarse_particles.items():
                self.assertTrue(coarse_split_particles[k].size() > 0)
                self.assertTrue(
                    coarse_split_particles[k].size()
                    == particles.size() * sim.refined_particle_nbr
                )

            for patch in level.patches:
                for pop_name in [
                    key for key in patch.patch_datas.keys() if key.endswith("particles")
                ]:
                    part1 = patch.patch_datas[pop_name].dataset.select(
                        patch.box
                    )  # drop ghosts
                    part2 = coarse_split_particles[pop_name].select(patch.box)
                    self.assertEqual(part1, part2)

    def _test_levelghostparticles_have_correct_split_from_coarser_particle(
        self, datahier
    ):
        dim = datahier.level(0).patches[0].box.ndim

        from pyphare.pharein.global_vars import sim

        assert sim is not None
        assert len(sim.cells) == dim

        particle_level_ghost_boxes_per_level = level_ghost_boxes(datahier, "particles")

        self.assertTrue(len(particle_level_ghost_boxes_per_level.items()) > 0)
        for ilvl, particle_gaboxes in particle_level_ghost_boxes_per_level.items():
            self.assertTrue(ilvl > 0)  # has no level 0

            lvlParticles = self._domainParticles_for(datahier, ilvl - 1)
            for pop_name, gaboxes_list in particle_gaboxes.items():
                coarse_particles = lvlParticles[pop_name]
                self.assertTrue(coarse_particles.size() > 0)

                coarse_split_particles = coarse_particles.split(sim)
                self.assertTrue(coarse_split_particles.size() > 0)
                self.assertTrue(
                    coarse_split_particles.size()
                    == coarse_particles.size() * sim.refined_particle_nbr
                )

                for gabox in gaboxes_list:
                    gabox_patchData = gabox["pdata"]

                    for ghostBox in gabox["boxes"]:
                        part1 = gabox_patchData.dataset.select(ghostBox)
                        part2 = coarse_split_particles.select(ghostBox)
                        self.assertEqual(part1, part2)
