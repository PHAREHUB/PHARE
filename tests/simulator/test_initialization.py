#
#  Common base class across hybrid and mhd tests
#   see
#      tests/simulator/initialize/test_init_mhd.py
#      tests/simulator/initialize/test_init_hybrid.py
#


import unittest
import numpy as np
from ddt import ddt

import pyphare.pharein as ph
from pyphare.core.box import nDBox
from pyphare.core.phare_utilities import assert_fp_any_all_close

from tests.simulator import SimulatorTest


@ddt
class InitializationTest(SimulatorTest):
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

        model = ph.global_vars.sim.model

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

        model = ph.global_vars.sim.model
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

        model = ph.global_vars.sim.model
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

            model = ph.global_vars.sim.model
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


if __name__ == "__main__":
    unittest.main()
