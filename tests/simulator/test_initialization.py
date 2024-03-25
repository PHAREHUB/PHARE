from pyphare.cpp import cpp_lib

cpp = cpp_lib()

import unittest

import numpy as np
from ddt import ddt
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
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
from pyphare.pharesee.particles import aggregate as aggregate_particles
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest


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
                quantity=quantity,
                write_timestamps=np.zeros(time_step_nbr),
                compute_timestamps=np.zeros(time_step_nbr),
            )

        for quantity in ["density", "bulkVelocity"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=np.zeros(time_step_nbr),
                compute_timestamps=np.zeros(time_step_nbr),
            )

        poplist = ["protons", "beam"] if beam else ["protons"]
        for pop in poplist:
            for quantity in ["density", "flux"]:
                FluidDiagnostics(
                    quantity=quantity,
                    write_timestamps=np.zeros(time_step_nbr),
                    compute_timestamps=np.zeros(time_step_nbr),
                    population_name=pop,
                )

            for quantity in ["domain", "levelGhost", "patchGhost"]:
                ParticleDiagnostics(
                    quantity=quantity,
                    compute_timestamps=np.zeros(time_step_nbr),
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

        if qty == "particles":
            merge_particles(particle_hier)

        if is_particle_type:
            return particle_hier

        if qty == "moments":
            mom_hier = hierarchy_from(h5_filename=diag_outputs + "/ions_density.h5")
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

    def _test_B_is_as_provided_by_user(self, dim, interp_order, **kwargs):
        print(
            "test_B_is_as_provided_by_user : dim  {} interp_order : {}".format(
                dim, interp_order
            )
        )
        hier = self.getHierarchy(
            dim,
            interp_order,
            refinement_boxes=None,
            qty="b",
            **kwargs,
        )

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
                    raise ValueError("Unsupported dimension")

    def _test_bulkvel_is_as_provided_by_user(self, dim, interp_order):
        hier = self.getHierarchy(
            dim,
            interp_order,
            {"L0": {"B0": nDBox(dim, 10, 19)}},
            "moments",
            nbr_part_per_cell=100,
            beam=True,
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

                layout = patch.patch_datas["protons_Fx"].layout
                centering = layout.centering["X"][
                    patch.patch_datas["protons_Fx"].field_name
                ]
                nbrGhosts = layout.nbrGhosts(interp_order, centering)

                if dim == 1:
                    x = patch.patch_datas["protons_Fx"].x[nbrGhosts:-nbrGhosts]

                    fpx = patch.patch_datas["protons_Fx"].dataset[nbrGhosts:-nbrGhosts]
                    fpy = patch.patch_datas["protons_Fy"].dataset[nbrGhosts:-nbrGhosts]
                    fpz = patch.patch_datas["protons_Fz"].dataset[nbrGhosts:-nbrGhosts]

                    fbx = patch.patch_datas["beam_Fx"].dataset[nbrGhosts:-nbrGhosts]
                    fby = patch.patch_datas["beam_Fy"].dataset[nbrGhosts:-nbrGhosts]
                    fbz = patch.patch_datas["beam_Fz"].dataset[nbrGhosts:-nbrGhosts]

                    ni = patch.patch_datas["rho"].dataset[nbrGhosts:-nbrGhosts]

                    vxact = (fpx + fbx) / ni
                    vyact = (fpy + fby) / ni
                    vzact = (fpz + fbz) / ni

                    vxexp = (nprot(x) * vx_fn(x) + nbeam(x) * vx_fn(x)) / (
                        nprot(x) + nbeam(x)
                    )
                    vyexp = (nprot(x) * vy_fn(x) + nbeam(x) * vy_fn(x)) / (
                        nprot(x) + nbeam(x)
                    )
                    vzexp = (nprot(x) * vz_fn(x) + nbeam(x) * vz_fn(x)) / (
                        nprot(x) + nbeam(x)
                    )

                    for vexp, vact in zip((vxexp, vyexp, vzexp), (vxact, vyact, vzact)):
                        std = np.std(vexp - vact)
                        print("sigma(user v - actual v) = {}".format(std))
                        self.assertTrue(
                            std < 1e-2
                        )  # empirical value obtained from print just above

                def reshape(patch_data, nGhosts):
                    return patch_data.dataset[:].reshape(
                        patch.box.shape + (nGhosts * 2) + 1
                    )

                if dim == 2:
                    xx, yy = np.meshgrid(
                        patch.patch_datas["protons_Fx"].x,
                        patch.patch_datas["protons_Fx"].y,
                        indexing="ij",
                    )

                    density = reshape(patch.patch_datas["rho"], nbrGhosts)

                    protons_Fx = reshape(patch.patch_datas["protons_Fx"], nbrGhosts)
                    protons_Fy = reshape(patch.patch_datas["protons_Fy"], nbrGhosts)
                    protons_Fz = reshape(patch.patch_datas["protons_Fz"], nbrGhosts)

                    beam_Fx = reshape(patch.patch_datas["beam_Fx"], nbrGhosts)
                    beam_Fy = reshape(patch.patch_datas["beam_Fy"], nbrGhosts)
                    beam_Fz = reshape(patch.patch_datas["beam_Fz"], nbrGhosts)

                    x = xx[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]
                    y = yy[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]

                    fpx = protons_Fx[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]
                    fpy = protons_Fy[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]
                    fpz = protons_Fz[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]

                    fbx = beam_Fx[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]
                    fby = beam_Fy[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]
                    fbz = beam_Fz[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]

                    ni = density[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]

                    vxact = (fpx + fbx) / ni
                    vyact = (fpy + fby) / ni
                    vzact = (fpz + fbz) / ni

                    vxexp = (nprot(x, y) * vx_fn(x, y) + nbeam(x, y) * vx_fn(x, y)) / (
                        nprot(x, y) + nbeam(x, y)
                    )
                    vyexp = (nprot(x, y) * vy_fn(x, y) + nbeam(x, y) * vy_fn(x, y)) / (
                        nprot(x, y) + nbeam(x, y)
                    )
                    vzexp = (nprot(x, y) * vz_fn(x, y) + nbeam(x, y) * vz_fn(x, y)) / (
                        nprot(x, y) + nbeam(x, y)
                    )

                    for vexp, vact in zip((vxexp, vyexp, vzexp), (vxact, vyact, vzact)):
                        self.assertTrue(np.std(vexp - vact) < 1e-2)

    def _test_density_is_as_provided_by_user(self, dim, interp_order):
        empirical_dim_devs = {
            1: 6e-3,
            2: 3e-2,
        }

        nbParts = {1: 10000, 2: 1000}
        print(
            "test_density_is_as_provided_by_user : interp_order : {}".format(
                interp_order
            )
        )
        hier = self.getHierarchy(
            dim,
            interp_order,
            {"L0": {"B0": nDBox(dim, 10, 20)}},
            qty="moments",
            nbr_part_per_cell=nbParts[dim],
            beam=True,
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

                if dim == 1:
                    protons_expected = proton_density_fn(x[nbrGhosts:-nbrGhosts])
                    beam_expected = beam_density_fn(x[nbrGhosts:-nbrGhosts])
                    ion_expected = protons_expected + beam_expected

                    ion_actual = ion_density[nbrGhosts:-nbrGhosts]
                    beam_actual = beam_density[nbrGhosts:-nbrGhosts]
                    protons_actual = proton_density[nbrGhosts:-nbrGhosts]

                if dim == 2:
                    y = patch.patch_datas["rho"].y
                    xx, yy = np.meshgrid(x, y, indexing="ij")

                    x0 = xx[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]
                    y0 = yy[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]

                    protons_expected = proton_density_fn(x0, y0)
                    beam_expected = beam_density_fn(x0, y0)
                    ion_expected = protons_expected + beam_expected

                    ion_actual = ion_density[nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts]
                    beam_actual = beam_density[
                        nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts
                    ]
                    protons_actual = proton_density[
                        nbrGhosts:-nbrGhosts, nbrGhosts:-nbrGhosts
                    ]

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
                        dev, empirical_dim_devs[dim], f"{name} has dev = {dev}"
                    )

    def _test_density_decreases_as_1overSqrtN(
        self, dim, interp_order, nbr_particles=None, cells=960
    ):
        import matplotlib.pyplot as plt

        print(f"test_density_decreases_as_1overSqrtN, interp_order = {interp_order}")

        if nbr_particles is None:
            nbr_particles = np.asarray([100, 1000, 5000, 10000])

        noise = np.zeros(len(nbr_particles))

        for inbr, nbrpart in enumerate(nbr_particles):
            hier = self.getHierarchy(
                dim,
                interp_order,
                None,
                "moments",
                nbr_part_per_cell=nbrpart,
                diag_outputs=f"{nbrpart}",
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
            select = tuple([slice(nbrGhosts, -nbrGhosts) for i in range(dim)])
            ion_density = patch.patch_datas["rho"].dataset[:]
            mesh = patch.patch_datas["rho"].meshgrid(select)

            expected = density_fn(*mesh)
            actual = ion_density[select]
            noise[inbr] = np.std(expected - actual)
            print(f"noise is {noise[inbr]} for {nbrpart} particles per cell")

            if dim == 1:
                x = patch.patch_datas["rho"].x
                plt.figure()
                plt.plot(x[nbrGhosts:-nbrGhosts], actual, label="actual")
                plt.plot(x[nbrGhosts:-nbrGhosts], expected, label="expected")
                plt.legend()
                plt.title(r"$\sigma =$ {}".format(noise[inbr]))
                plt.savefig(f"noise_{nbrpart}_interp_{dim}_{interp_order}.png")
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
        plt.savefig(f"noise_nppc_interp_{dim}_{interp_order}.png")
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
        plt.savefig(f"noise_nppc_minus_theory_interp_{dim}_{interp_order}.png")
        plt.close("all")
        self.assertGreater(3e-2, noiseMinusTheory[1:].mean())

    def _test_nbr_particles_per_cell_is_as_provided(
        self, dim, interp_order, default_ppc=100
    ):
        datahier = self.getHierarchy(
            dim,
            interp_order,
            {"L0": {"B0": nDBox(dim, 10, 20)}},
            "particles",
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
            refinement_boxes,
            "particles",
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

    def _test_patch_ghost_on_refined_level_case(self, dim, has_patch_ghost, **kwargs):
        import pyphare.pharein as ph

        out = "phare_outputs"
        refinement_boxes = {"L0": [nDBox(dim, 10, 19)]}
        kwargs["interp_order"] = kwargs.get("interp_order", 1)
        kwargs["diag_outputs"] = f"{has_patch_ghost}"
        datahier = self.getHierarchy(
            ndim=dim,
            refinement_boxes=refinement_boxes,
            qty="particles_patch_ghost",
            **kwargs,
        )

        self.assertTrue(
            any(
                [
                    diagInfo.quantity.endswith("patchGhost")
                    for diagname, diagInfo in ph.global_vars.sim.diagnostics.items()
                ]
            )
        )
        self.assertTrue((1 in datahier.levels()) == has_patch_ghost)

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


if __name__ == "__main__":
    unittest.main()
