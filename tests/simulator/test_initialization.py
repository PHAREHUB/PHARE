#
#

import unittest
import numpy as np
from ddt import ddt

from pyphare.cpp import cpp_lib
from pyphare.core.box import nDBox
from pyphare.core.phare_utilities import assert_fp_any_all_close
from pyphare.pharein import ElectronModel, MaxwellianFluidModel
from pyphare.pharein.diagnostics import (
    ElectromagDiagnostics,
    FluidDiagnostics,
    ParticleDiagnostics,
)
from pyphare.pharein.simulation import Simulation
from pyphare.pharesee.geometry import level_ghost_boxes
from pyphare.pharesee.hierarchy.hierarchy_utils import merge_particles
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.particles import aggregate as aggregate_particles
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest


cpp = cpp_lib()


@ddt
class InitializationTest(SimulatorTest):
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
        diag_outputs="",
        nbr_part_per_cell=100,
        density=_density,
        extra_diag_options={},
        beam=False,
        time_step_nbr=1,
        smallest_patch_size=None,
        largest_patch_size=10,
        cells=120,
        dl=0.1,
        **kwargs,
    ):
        diag_outputs = self.unique_diag_dir_for_test_case(
            "phare_outputs/init", ndim, interp_order, diag_outputs
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

        def beam_density(*xyz):
            return np.zeros_like(xyz[0]) + 0.3

        def bx(*xyz):
            return 1.0

        def by(*xyz):
            L = sim.simulation_domain()
            _ = lambda i: 0.1 * np.cos(2 * np.pi * xyz[i] / L[i])
            return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

        def bz(*xyz):
            L = sim.simulation_domain()
            _ = lambda i: 0.1 * np.sin(2 * np.pi * xyz[i] / L[i])
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
            _ = lambda i: 0.1 * np.sin(2 * np.pi * xyz[i] / L[i])
            return np.asarray([_(i) for i in range(len(xyz))]).prod(axis=0)

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
            "density": density,
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
            MaxwellianFluidModel(
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
            MaxwellianFluidModel(bx=bx, by=by, bz=bz, protons=protons)

        ElectronModel(closure="isothermal", Te=0.12)

        for quantity in ["E", "B"]:
            ElectromagDiagnostics(
                quantity=quantity, write_timestamps=np.zeros(time_step_nbr)
            )

        for quantity in ["charge_density", "bulkVelocity"]:
            FluidDiagnostics(
                quantity=quantity, write_timestamps=np.zeros(time_step_nbr)
            )

        poplist = ["protons", "beam"] if beam else ["protons"]
        for pop in poplist:
            for quantity in ["density", "flux"]:
                FluidDiagnostics(
                    quantity=quantity,
                    write_timestamps=np.zeros(time_step_nbr),
                    population_name=pop,
                )

            for quantity in ["domain", "levelGhost"]:
                ParticleDiagnostics(
                    quantity=quantity,
                    write_timestamps=np.zeros(time_step_nbr),
                    population_name=pop,
                )

        Simulator(global_vars.sim).initialize().reset()

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

    def _test_B_is_as_provided_by_user(self, dim, interp_order, ppc=100, **kwargs):
        print(
            "test_B_is_as_provided_by_user : dim  {} interp_order : {}".format(
                dim, interp_order
            )
        )
        now = self.datetime_now()
        hier = self.getHierarchy(
            dim,
            interp_order,
            refinement_boxes=None,
            qty="b",
            nbr_part_per_cell=ppc,
            diag_outputs=f"test_b/{dim}/{interp_order}/{self.ddt_test_id()}",
            **kwargs,
        )
        print(
            f"\n{self._testMethodName}_{dim}d init took {self.datetime_diff(now)} seconds"
        )
        now = self.datetime_now()

        from pyphare.pharein import global_vars

        model = global_vars.sim.model

        bx_fn = model.model_dict["bx"]
        by_fn = model.model_dict["by"]
        bz_fn = model.model_dict["bz"]
        for ilvl, level in hier.levels().items():
            self.assertTrue(ilvl == 0)  # only level 0 is expected perfect precision
            print("checking level {}".format(ilvl))
            for patch in level.patches:
                bx_pd = patch.patch_datas["Bx"]
                by_pd = patch.patch_datas["By"]
                bz_pd = patch.patch_datas["Bz"]

                bx = bx_pd.dataset[:]
                by = by_pd.dataset[:]
                bz = bz_pd.dataset[:]

                xbx = bx_pd.x[:]
                xby = by_pd.x[:]
                xbz = bz_pd.x[:]

                if dim == 1:
                    # discrepancy in 1d for some reason : https://github.com/PHAREHUB/PHARE/issues/580
                    assert_fp_any_all_close(bx, bx_fn(xbx), atol=1e-15, rtol=0)
                    assert_fp_any_all_close(by, by_fn(xby), atol=1e-15, rtol=0)
                    assert_fp_any_all_close(bz, bz_fn(xbz), atol=1e-15, rtol=0)

                if dim >= 2:
                    ybx = bx_pd.y[:]
                    yby = by_pd.y[:]
                    ybz = bz_pd.y[:]

                if dim == 2:
                    xbx, ybx = [
                        a.flatten() for a in np.meshgrid(xbx, ybx, indexing="ij")
                    ]
                    xby, yby = [
                        a.flatten() for a in np.meshgrid(xby, yby, indexing="ij")
                    ]
                    xbz, ybz = [
                        a.flatten() for a in np.meshgrid(xbz, ybz, indexing="ij")
                    ]

                    assert_fp_any_all_close(bx, bx_fn(xbx, ybx), atol=1e-16, rtol=0)
                    assert_fp_any_all_close(
                        by, by_fn(xby, yby).reshape(by.shape), atol=1e-16, rtol=0
                    )
                    assert_fp_any_all_close(
                        bz, bz_fn(xbz, ybz).reshape(bz.shape), atol=1e-16, rtol=0
                    )

                if dim == 3:
                    zbx = bx_pd.z[:]
                    zby = by_pd.z[:]
                    zbz = bz_pd.z[:]

                    xbx, ybx, zbx = [
                        a.flatten() for a in np.meshgrid(xbx, ybx, zbx, indexing="ij")
                    ]
                    xby, yby, zby = [
                        a.flatten() for a in np.meshgrid(xby, yby, zby, indexing="ij")
                    ]
                    xbz, ybz, zbz = [
                        a.flatten() for a in np.meshgrid(xbz, ybz, zbz, indexing="ij")
                    ]

                    np.testing.assert_allclose(
                        bx, bx_fn(xbx, ybx, zbx), atol=1e-16, rtol=0
                    )
                    np.testing.assert_allclose(
                        by, by_fn(xby, yby, zby).reshape(by.shape), atol=1e-16, rtol=0
                    )
                    np.testing.assert_allclose(
                        bz, bz_fn(xbz, ybz, zbz).reshape(bz.shape), atol=1e-16, rtol=0
                    )

        print(f"\n{self._testMethodName}_{dim}d took {self.datetime_diff(now)} seconds")

    def _test_bulkvel_is_as_provided_by_user(
        self, dim, interp_order, ppc=100, **kwargs
    ):
        hier = self.getHierarchy(
            dim,
            interp_order,
            {"L0": {"B0": nDBox(dim, 10, 19)}},
            "moments",
            nbr_part_per_cell=ppc,
            beam=True,
            diag_outputs=f"test_bulkV/{dim}/{interp_order}/{self.ddt_test_id()}",
            **kwargs,
        )

        from pyphare.pharein import global_vars

        model = global_vars.sim.model
        # protons and beam have same bulk vel here so take only proton func.
        vx_fn = model.model_dict["protons"]["vx"]
        vy_fn = model.model_dict["protons"]["vy"]
        vz_fn = model.model_dict["protons"]["vz"]
        nprot = model.model_dict["protons"]["density"]
        nbeam = model.model_dict["beam"]["density"]

        for ilvl, level in hier.levels().items():
            print("checking density on level {}".format(ilvl))
            for ip, patch in enumerate(level.patches):
                print("patch {}".format(ip))

                pdatas = patch.patch_datas
                layout = pdatas["protons_Fx"].layout
                centering = layout.centering["X"][pdatas["protons_Fx"].field_name]
                nbrGhosts = layout.nbrGhosts(
                    interp_order, centering
                )  # primal in all directions
                select = tuple([slice(nbrGhosts, -nbrGhosts) for i in range(dim)])

                def domain(patch_data):
                    if dim == 1:
                        return patch_data.dataset[select]
                    return patch_data.dataset[:].reshape(
                        patch.box.shape + (nbrGhosts * 2) + 1
                    )[select]

                ni = domain(pdatas["rho"])
                vxact = (domain(pdatas["protons_Fx"]) + domain(pdatas["beam_Fx"])) / ni
                vyact = (domain(pdatas["protons_Fy"]) + domain(pdatas["beam_Fy"])) / ni
                vzact = (domain(pdatas["protons_Fz"]) + domain(pdatas["beam_Fz"])) / ni

                select = pdatas["protons_Fx"].meshgrid(select)
                vxexp = (
                    nprot(*select) * vx_fn(*select) + nbeam(*select) * vx_fn(*select)
                ) / (nprot(*select) + nbeam(*select))
                vyexp = (
                    nprot(*select) * vy_fn(*select) + nbeam(*select) * vy_fn(*select)
                ) / (nprot(*select) + nbeam(*select))
                vzexp = (
                    nprot(*select) * vz_fn(*select) + nbeam(*select) * vz_fn(*select)
                ) / (nprot(*select) + nbeam(*select))
                for vexp, vact in zip((vxexp, vyexp, vzexp), (vxact, vyact, vzact)):
                    self.assertTrue(np.std(vexp - vact) < 1e-2)

    def _test_density_is_as_provided_by_user(self, ndim, interp_order, **kwargs):
        print(
            f"test_density_is_as_provided_by_user : dim {ndim} interp_order {interp_order}"
        )

        empirical_dim_devs = {
            1: 6e-3,
            2: 3e-2,
            3: 2e-1,
        }
        nbParts = {1: 10000, 2: 1000, 3: 20}
        hier = self.getHierarchy(
            ndim,
            interp_order,
            {"L0": {"B0": nDBox(ndim, 5, 14)}},
            qty="moments",
            nbr_part_per_cell=nbParts[ndim],
            beam=True,
            **kwargs,
        )

        from pyphare.pharein import global_vars

        model = global_vars.sim.model
        proton_density_fn = model.model_dict["protons"]["density"]
        beam_density_fn = model.model_dict["beam"]["density"]

        for ilvl, level in hier.levels().items():
            print("checking density on level {}".format(ilvl))
            for ip, patch in enumerate(level.patches):
                print("patch {}".format(ip))

                ion_density = patch.patch_datas["rho"].dataset[:]
                proton_density = patch.patch_datas["protons_rho"].dataset[:]
                beam_density = patch.patch_datas["beam_rho"].dataset[:]
                x = patch.patch_datas["rho"].x

                layout = patch.patch_datas["rho"].layout
                centering = layout.centering["X"][patch.patch_datas["rho"].field_name]
                nbrGhosts = layout.nbrGhosts(interp_order, centering)
                select = tuple([slice(nbrGhosts, -nbrGhosts) for i in range(ndim)])

                mesh = patch.patch_datas["rho"].meshgrid(select)
                protons_expected = proton_density_fn(*mesh)
                beam_expected = beam_density_fn(*mesh)
                ion_expected = protons_expected + beam_expected

                protons_actual = proton_density[select]
                beam_actual = beam_density[select]
                ion_actual = ion_density[select]

                names = ("ions", "protons", "beam")
                expected = (ion_expected, protons_expected, beam_expected)
                actual = (ion_actual, protons_actual, beam_actual)
                devs = {
                    name: np.std(expected - actual)
                    for name, expected, actual in zip(names, expected, actual)
                }

                for name, dev in devs.items():
                    print(f"sigma(user density - {name} density) = {dev}")
                    self.assertLess(
                        dev, empirical_dim_devs[ndim], f"{name} has dev = {dev}"
                    )

    def _test_density_decreases_as_1overSqrtN(
        self, ndim, interp_order, nbr_particles=None, cells=960
    ):
        import matplotlib.pyplot as plt

        if nbr_particles is None:
            nbr_particles = np.asarray([100, 1000, 5000, 10000])

        print(
            f"test_density_decreases_as_1overSqrtN, interp_order = {interp_order} {nbr_particles}"
        )

        noise = np.zeros(len(nbr_particles))

        for inbr, nbrpart in enumerate(nbr_particles):
            hier = self.getHierarchy(
                ndim,
                interp_order,
                None,
                "moments",
                nbr_part_per_cell=nbrpart,
                diag_outputs=f"1overSqrtN/{ndim}/{interp_order}/{nbrpart}",
                density=lambda *xyz: np.zeros(tuple(_.shape[0] for _ in xyz)) + 1.0,
                smallest_patch_size=int(cells / 2),
                largest_patch_size=int(cells / 2),
                cells=cells,
                dl=0.0125,
            )

            from pyphare.pharein import global_vars

            model = global_vars.sim.model
            density_fn = model.model_dict["protons"]["density"]

            patch = hier.level(0).patches[0]
            layout = patch.patch_datas["rho"].layout

            centering = layout.centering["X"][patch.patch_datas["rho"].field_name]
            nbrGhosts = layout.nbrGhosts(interp_order, centering)
            select = tuple([slice(nbrGhosts, -nbrGhosts) for i in range(ndim)])
            ion_density = patch.patch_datas["rho"].dataset[:]
            mesh = patch.patch_datas["rho"].meshgrid(select)

            expected = density_fn(*mesh)
            actual = ion_density[select]
            noise[inbr] = np.std(expected - actual)
            print(f"noise is {noise[inbr]} for {nbrpart} particles per cell")

            if ndim == 1:
                x = patch.patch_datas["rho"].x
                plt.figure()
                plt.plot(x[select], actual, label="actual")
                plt.plot(x[select], expected, label="expected")
                plt.legend()
                plt.title(r"$\sigma =$ {}".format(noise[inbr]))
                plt.savefig(f"noise_{nbrpart}_interp_{ndim}_{interp_order}.png")
                plt.close("all")

        plt.figure()
        plt.plot(nbr_particles, noise / noise[0], label=r"$\sigma/\sigma_0$")
        plt.plot(
            nbr_particles,
            1 / np.sqrt(nbr_particles / nbr_particles[0]),
            label=r"$1/sqrt(nppc/nppc0)$",
        )
        plt.xlabel("nbr_particles")
        plt.legend()
        plt.savefig(f"noise_nppc_interp_{ndim}_{interp_order}.png")
        plt.close("all")

        noiseMinusTheory = noise / noise[0] - 1 / np.sqrt(
            nbr_particles / nbr_particles[0]
        )
        plt.figure()
        plt.plot(
            nbr_particles,
            noiseMinusTheory,
            label=r"$\sigma/\sigma_0 - 1/sqrt(nppc/nppc0)$",
        )
        plt.xlabel("nbr_particles")
        plt.legend()
        plt.savefig(f"noise_nppc_minus_theory_interp_{ndim}_{interp_order}.png")
        plt.close("all")
        self.assertGreater(3e-2, noiseMinusTheory[1:].mean())

    def _test_nbr_particles_per_cell_is_as_provided(
        self, ndim, interp_order, ppc=100, **kwargs
    ):
        ddt_test_id = self.ddt_test_id()
        datahier = self.getHierarchy(
            ndim,
            interp_order,
            {},
            "particles",
            diag_outputs=f"ppc/{ndim}/{interp_order}/{ddt_test_id}",
            nbr_part_per_cell=ppc,
            **kwargs,
        )

        if cpp.mpi_rank() > 0:
            return

        for pi, patch in enumerate(datahier.level(0).patches):
            pd = patch.patch_datas["protons_particles"]
            icells = pd.dataset[patch.box].iCells
            H, edges = np.histogramdd(icells, bins=patch.box.shape)
            self.assertTrue((H == ppc).all())

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
            refinement_boxes,
            "particles",
            cells=30,
            **kwargs,
        )

        if cpp.mpi_rank() > 0:
            return

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
        ndim = datahier.level(0).patches[0].box.ndim

        from pyphare.pharein.global_vars import sim

        assert sim is not None
        assert len(sim.cells) == ndim

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


if __name__ == "__main__":
    unittest.main()
