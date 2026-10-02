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

from tests.simulator import SimulatorTest, basicSimulatorArgs


def vector_potential_init(ndim, L):
    """
    Magnetic init kwargs in vector-potential mode, for 2D and 3D.

    A is smooth with a periodic curl, but has linear terms (uniform in-plane B)
    that make A itself non-periodic. 2D: az and a direct bz. 3D: ax, ay, az.
    """
    k = 2 * np.pi / np.asarray(L)

    if ndim == 2:

        def az(x, y):
            return (
                0.5 * y
                - 0.3 * x
                + 0.1 * np.cos(k[0] * x) * np.cos(k[1] * y)
                + 0.05 * np.sin(2 * k[0] * x) * np.sin(k[1] * y)
            )

        def bz(x, y):
            return 0.2 + 0.1 * np.sin(k[0] * x) * np.cos(k[1] * y)

        return dict(az=az, bz=bz)

    assert ndim == 3

    def ax(x, y, z):
        return 0.2 * z + 0.1 * np.sin(k[1] * y) * np.cos(k[2] * z)

    def ay(x, y, z):
        return 0.3 * x + 0.1 * np.sin(k[2] * z) * np.cos(k[0] * x)

    def az(x, y, z):
        return 0.5 * y + 0.1 * np.cos(k[0] * x) * np.cos(k[1] * y) * np.sin(k[2] * z)

    return dict(ax=ax, ay=ay, az=az)


def _yee_coords(patch):
    """primal and dual node coordinates per direction, ghosts included, from B"""
    bx, by, bz = (patch.patch_datas[f"B{c}"] for c in "xyz")
    primal = [bx.x, by.y] + ([bz.z] if bx.ndim == 3 else [])
    dual = [by.x, bx.y] + ([bx.z] if bx.ndim == 3 else [])
    return primal, dual


def _curl_of_edge_sampled(patch, vecpot):
    """
    numpy discrete curl of A sampled once on the Yee edge lattice of the patch
    (ghosts included), with the solver's two-point stencil. Returns
    {"Bx", "By", "Bz"} on the B ghost boxes.
    """
    primal, dual = _yee_coords(patch)
    ndim = len(primal)
    dl = patch.patch_datas["Bx"].dl

    def on_edges(comp):  # A_comp is dual along comp, primal elsewhere
        coords = [dual[d] if d == comp else primal[d] for d in range(ndim)]
        mesh = np.meshgrid(*coords, indexing="ij")
        return np.zeros(mesh[0].shape) + vecpot["a" + "xyz"[comp]](*mesh)

    Ax, Ay, Az = (on_edges(c) for c in range(3))

    def D(a, d):
        return np.diff(a, axis=d) / dl[d]

    if ndim == 2:
        return {"Bx": D(Az, 1), "By": -D(Az, 0), "Bz": D(Ay, 0) - D(Ax, 1)}
    return {
        "Bx": D(Az, 1) - D(Ay, 2),
        "By": D(Ax, 2) - D(Az, 0),
        "Bz": D(Ay, 0) - D(Ax, 1),
    }


def _domain(pd, data=None):
    data = pd.dataset[:] if data is None else data
    return data[tuple(slice(g, -g) for g in pd.ghosts_nbr)]


def _divB_domain(patch):
    """numpy discrete div B on the patch domain cells"""
    pds = [patch.patch_datas[f"B{c}"] for c in "xyz"][: patch.box.ndim]
    divB = sum(
        np.diff(pd.dataset[:].astype(np.float64), axis=d) / pd.dl[d]
        for d, pd in enumerate(pds)
    )
    return divB[tuple(slice(g, -g) for g in pds[0].ghosts_nbr)]


@ddt
class InitializationTest(SimulatorTest):
    def _test_B_is_as_provided_by_user(self, dim, interp_order, **kwargs):
        print(
            "test_B_is_as_provided_by_user : dim  {} interp_order : {}".format(
                dim, interp_order
            )
        )
        now = self.datetime_now()
        hier = self.getHierarchy(
            dim,
            interp_order,
            qty="b",
            refinement_boxes=None,
            diag_outputs=f"test_b/{dim}/{interp_order}/{self.ddt_test_id()}",
            **kwargs,
        )
        print(
            f"\n{self._testMethodName}_{dim}d init took {self.datetime_diff(now)} seconds"
        )
        now = self.datetime_now()

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

                    assert_fp_any_all_close(
                        bx, bx_fn(xbx, ybx, zbx), atol=1e-16, rtol=0
                    )
                    assert_fp_any_all_close(
                        by, by_fn(xby, yby, zby).reshape(by.shape), atol=1e-16, rtol=0
                    )
                    assert_fp_any_all_close(
                        bz, bz_fn(xbz, ybz, zbz).reshape(bz.shape), atol=1e-16, rtol=0
                    )

        print(f"\n{self._testMethodName}_{dim}d took {self.datetime_diff(now)} seconds")

    def _test_B_from_vector_potential(
        self, dim, interp_order, refinement_boxes, **kwargs
    ):
        """
        A mode: level 0 B is the discrete curl of A sampled on the Yee edges, and the
        discrete div B is at round-off on level 0 and on the refined level 1.
        """
        print(f"test_B_from_vector_potential : dim {dim} interp_order : {interp_order}")
        hier = self.getHierarchy(
            dim,
            interp_order,
            qty="b",
            refinement_boxes=refinement_boxes,
            diag_outputs=f"test_b_vecpot/{dim}/{interp_order}/{self.ddt_test_id()}",
            vecpot=True,
            **kwargs,
        )

        model = ph.global_vars.sim.model
        self.assertTrue("bx" not in model.model_dict)
        vecpot = model.model_dict["vector_potential"]
        direct = model.model_dict["direct_b"]
        self.assertEqual(sorted(direct.keys()), ["z"] if dim == 2 else [])

        levels = hier.levels()
        self.assertTrue(len(levels) > 1)

        a_max, b_max, curl_err = 0.0, 0.0, 0.0
        for patch in levels[0].patches:
            curl = _curl_of_edge_sampled(patch, vecpot)
            primal, dual = _yee_coords(patch)
            a_max = max(
                a_max,
                *[
                    np.max(np.abs(fn(*np.meshgrid(*primal, indexing="ij"))))
                    for fn in vecpot.values()
                ],
            )

            for name, expected in curl.items():
                pd = patch.patch_datas[name]
                actual = pd.dataset[:].astype(np.float64)
                self.assertEqual(actual.shape, tuple(pd.size))
                if name == "Bz" and "z" in direct:
                    expected = np.zeros(actual.shape) + direct["z"](*pd.meshgrid())
                self.assertEqual(actual.shape, expected.shape)
                b_max = max(b_max, np.max(np.abs(actual)))
                curl_err = max(
                    curl_err, np.max(np.abs(_domain(pd) - _domain(pd, expected)))
                )

        # fp32 diagnostics (no PHARE_DIAG_DOUBLES) are bounded by the dump rounding of B
        eps = np.finfo(levels[0].patches[0].patch_datas["Bx"].dataset.dtype).eps
        dl0 = min(levels[0].patches[0].patch_datas["Bx"].dl)
        curl_atol = max(1e-13 * a_max / dl0, 4 * eps * b_max)
        print(f"L0 max|B - curlA| = {curl_err:.3e} (atol {curl_atol:.3e})")
        self.assertLess(curl_err, curl_atol)

        for ilvl, level in levels.items():
            dl = min(level.patches[0].patch_datas["Bx"].dl)
            div_max = max(np.max(np.abs(_divB_domain(p))) for p in level.patches)
            div_atol = max(1e-12 * a_max / dl**2, 4 * dim * eps * b_max / dl)
            print(f"L{ilvl} max|divB| = {div_max:.3e} (atol {div_atol:.3e})")
            self.assertLess(div_max, div_atol)

    def _test_vector_potential_rejections(self, dim):
        """invalid B / A combinations raise ValueError, for both models"""
        fn = {
            1: lambda x: x * 0,
            2: lambda x, y: x * 0,
            3: lambda x, y, z: x * 0,
        }[dim]

        if dim == 1:
            invalid = [dict(ax=fn), dict(ay=fn), dict(az=fn), dict(bx=fn, az=fn)]
        elif dim == 2:
            invalid = [
                dict(az=fn, bx=fn),
                dict(az=fn, by=fn),
                dict(az=fn, bz=fn, ax=fn),
                dict(bz=fn, ay=fn),
            ]
        else:
            invalid = [dict(az=fn, bx=fn), dict(ax=fn, by=fn), dict(az=fn, bz=fn)]

        for model in (ph.MaxwellianFluidModel, ph.MHDModel):
            for magnetic in invalid:
                ph.global_vars.sim = None
                ph.Simulation(**basicSimulatorArgs(dim, 1))
                with self.assertRaises(ValueError, msg=f"{model.__name__} {magnetic}"):
                    model(**magnetic)
        ph.global_vars.sim = None

    def _test_bulkvel_is_as_provided_by_user(self, dim, interp_order, **kwargs):
        hier = self.getHierarchy(
            dim,
            interp_order,
            "moments",
            {"L0": {"B0": nDBox(dim, 10, 18)}},
            beam=True,
            diag_outputs=f"test_bulkV/{dim}/{interp_order}/{self.ddt_test_id()}",
            **kwargs,
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
        empirical_dim_devs = {
            1: 6e-3,
            2: 3e-2,
            3: 2e-1,
        }
        nbParts = {1: 10000, 2: 1000, 3: 20}
        hier = self.getHierarchy(
            ndim,
            interp_order,
            "moments",
            {"L0": {"B0": nDBox(ndim, 5, 14)}},
            nbr_part_per_cell=nbParts[ndim],
            beam=True,
            **kwargs,
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
                "moments",
                None,
                nbr_part_per_cell=nbrpart,
                diag_outputs=f"{nbrpart}",
                density=lambda *xyz: np.zeros(tuple(_.shape[0] for _ in xyz)) + 1.0,
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
