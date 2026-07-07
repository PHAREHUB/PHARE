#
# MHD magnetic-field initialization tests.
#
# Exercises every user path of the B = B0 + B1 split through MHDModel: direct components
# (bx/by/bz, b0*, b1*) and vector potential (ax/ay/az, a0*, a1*, B = curl(A)). B0 and B1 are read
# back from their own diagnostics (EM_B0 / EM_B1, raw fields, no reconstruction) and compared to
# the expected analytical fields. For the potential init we additionally assert the resulting
# B0 / B1 are discretely divergence-free (div . curl = 0).
#

import os
import numpy as np


import pyphare.pharein as ph
from pyphare.core import phare_utilities as phut
from pyphare.core.phare_utilities import assert_fp_any_all_close
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.run.utils import _divB2D, _divB3D

from tests.simulator.test_initialization import InitializationTest


def _const(c):
    """analytical field function returning the constant c, broadcast to the sample points."""
    return lambda *xyz: c + 0.0 * xyz[0]


def linear_potential(coeffs, ndim):
    """Build a *linear* vector potential A and its (constant) analytical curl B = rot(A).

    A_i = sum_j coeffs[ij] * x_j, with ij in {xy,xz,yx,yz,zx,zy} (diagonal terms drop out of the
    curl, so they are not used). Because A is linear, the discrete curl equals the analytical curl
    exactly, which lets the potential init be checked to machine precision.

        Bx = dAz/dy - dAy/dz = zy - yz
        By = dAx/dz - dAz/dx = xz - zx
        Bz = dAy/dx - dAx/dy = yx - xy

    In 2D there is no z coordinate, so the z-derivative terms (xz, yz) and z-dependence vanish.
    Returns (ax, ay, az) functions and the (Bx, By, Bz) constants.
    """
    c = {k: coeffs.get(k, 0.0) for k in ("xy", "xz", "yx", "yz", "zx", "zy")}
    if ndim == 2:
        c["xz"] = c["yz"] = 0.0  # no z coordinate in 2D

    def ax(*r):
        v = c["xy"] * r[1]
        if ndim == 3:
            v = v + c["xz"] * r[2]
        return v + 0.0 * r[0]

    def ay(*r):
        v = c["yx"] * r[0]
        if ndim == 3:
            v = v + c["yz"] * r[2]
        return v + 0.0 * r[0]

    def az(*r):
        return c["zx"] * r[0] + c["zy"] * r[1] + 0.0 * r[0]

    B = (c["zy"] - c["yz"], c["xz"] - c["zx"], c["yx"] - c["xy"])
    return (ax, ay, az), B


def tilted_potential(ndim, L):
    """Non-separable (tilted) sinusoidal vector potential. Its discrete curl is a non-trivial,
    spatially varying B that is still divergence-free by construction (div . curl = 0). Used only
    for the div-free test (the discrete curl of a non-linear A differs from rot(A) by truncation,
    so it is not value-checked)."""
    ks = [2.0 * np.pi / L * (i + 1) for i in range(ndim)]

    def phase(r):
        return sum(ks[i] * r[i] for i in range(ndim))

    # distinct amplitudes/offsets per component so no component is trivially zero
    def ax(*r):
        return 0.7 * np.sin(phase(r) + 0.3) + 0.0 * r[0]

    def ay(*r):
        return 0.5 * np.sin(phase(r) + 1.1) + 0.0 * r[0]

    def az(*r):
        return 0.9 * np.sin(phase(r) + 2.0) + 0.0 * r[0]

    return (ax, ay, az)


def mhd_component_cases(ndim):
    """All direct-component (b/b0/b1) user paths: (name, model_kwargs, B0, B1)."""
    b0v = (0.3, 0.2, 0.1)
    b1v = (0.05, -0.1, 0.2)
    zero = (0.0, 0.0, 0.0)
    return [
        # (name, model kwargs, expected B0, expected B1)
        ("default", {}, zero, (1.0, 0.0, 0.0)),
        (
            "total_b",
            dict(bx=_const(0.3), by=_const(0.2), bz=_const(0.1)),
            zero,
            (0.3, 0.2, 0.1),
        ),
        (
            "b0_only",
            dict(b0x=_const(b0v[0]), b0y=_const(b0v[1]), b0z=_const(b0v[2])),
            b0v,
            zero,
        ),
        (
            "b0_and_b1",
            dict(
                b0x=_const(b0v[0]),
                b0y=_const(b0v[1]),
                b0z=_const(b0v[2]),
                b1x=_const(b1v[0]),
                b1y=_const(b1v[1]),
                b1z=_const(b1v[2]),
            ),
            b0v,
            b1v,
        ),
        (
            "b1_only",
            dict(b1x=_const(b1v[0]), b1y=_const(b1v[1]), b1z=_const(b1v[2])),
            zero,
            b1v,
        ),
    ]


def mhd_potential_cases(ndim):
    """All vector-potential (a/a0/a1) user paths with linear potentials (exact curl match):
    (name, model_kwargs, B0, B1)."""
    zero = (0.0, 0.0, 0.0)
    # coefficients chosen so every curl component is non-zero (in 3D)
    (a1x, a1y, a1z), B1 = linear_potential(
        {"xz": 0.2, "yx": 0.3, "zy": 0.1}, ndim
    )
    (a0x, a0y, a0z), B0 = linear_potential(
        {"xy": 0.3, "yz": 0.2, "zx": 0.5}, ndim
    )
    return [
        (
            "plain_a",  # plain a folds into the perturbation, B0 = 0
            dict(ax=a1x, ay=a1y, az=a1z),
            zero,
            B1,
        ),
        ("a0_only", dict(a0x=a0x, a0y=a0y, a0z=a0z), B0, zero),
        ("a1_only", dict(a1x=a1x, a1y=a1y, a1z=a1z), zero, B1),
        (
            "a0_and_a1",
            dict(a0x=a0x, a0y=a0y, a0z=a0z, a1x=a1x, a1y=a1y, a1z=a1z),
            B0,
            B1,
        ),
    ]


class MHDInitializationTest(InitializationTest):
    def getHierarchy(
        self,
        ndim,
        interp_order,  # torm?
        qty,
        refinement_boxes={},
        density=None,
        time_step_nbr=1,
        time_step=0.001,
        smallest_patch_size=None,
        largest_patch_size=10,
        cells=120,
        dl=0.1,
        hall=False,
        res=False,
        hyper_res=False,
        extra_diag_options=None,
        timestamps=None,
        diag_outputs="",
        max_mhd_level=1,
        **kwargs,
    ):
        """Legacy total-B init hierarchy (compares EM_B to the total field model_dict[bx/by/bz]),
        used by the shared _test_B_is_as_provided_by_user. New per-path B0/B1 tests use
        getMHDHierarchies instead."""
        from tests.diagnostic import all_timestamps

        if smallest_patch_size is None:
            from pyphare.pharein.simulation import check_patch_size

            _, smallest_patch_size = check_patch_size(
                ndim, interp_order=interp_order, cells=cells
            )

        base_diag_dir = "phare_outputs/init_mhd"
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
            time_step=time_step,
            boundary_types=["periodic"] * ndim,
            cells=phut.np_array_ify(cells, ndim),
            dl=phut.np_array_ify(dl, ndim),
            interp_order=interp_order,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5", "options": extra_diag_options},
            strict=True,
            nesting_buffer=1,
            hyper_mode="spatial",
            eta=0.0,
            nu=0.02,
            gamma=5.0 / 3.0,
            reconstruction="Linear",
            limiter="VanLeer",
            riemann="Rusanov",
            mhd_timestepper="TVDRK2",
            hall=hall,
            res=res,
            hyper_res=hyper_res,
            model_options=["MHDModel"],
            max_mhd_level=max_mhd_level,
        )
        diag_outputs = sim.diag_options["options"]["dir"]
        L = sim.simulation_domain()

        def _density(*xyz):
            hL = np.array(sim.simulation_domain()) / 2
            return 0.3 + np.exp(
                sum([-((xyz[i] - hL[i]) ** 2) for i in range(len(xyz))])
            )

        def bx(*xyz):
            return 1.0

        def by(*xyz):
            return np.asarray(
                [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
            ).prod(axis=0)

        def bz(*xyz):
            return np.asarray(
                [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
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
                [0.1 * np.cos(2 * np.pi * xyz[i] / L[i]) for i in range(len(xyz))]
            ).prod(axis=0)

        def p(*xyz):
            return 1.0

        ph.MHDModel(
            density=density or _density, vx=vx, vy=vy, vz=vz, bx=bx, by=by, bz=bz, p=p
        )

        if timestamps is None:
            timestamps = all_timestamps(ph.global_vars.sim)

        ph.ElectromagDiagnostics(quantity="B", write_timestamps=timestamps)

        for quantity in ["rho", "V", "P"]:
            ph.MHDDiagnostics(quantity=quantity, write_timestamps=timestamps)

        Simulator(sim).initialize().reset()

        eb_hier = None
        if qty in ["b", "eb", "fields"]:
            eb_hier = hierarchy_from(h5_filename=diag_outputs + "/EM_B.h5", hier=eb_hier)
        if qty in ["e", "b", "eb"]:
            return eb_hier

        if qty == "moments" or qty == "fields":
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_charge_density.h5", hier=eb_hier
            )
            mom_hier = hierarchy_from(
                h5_filename=diag_outputs + "/ions_bulkVelocity.h5", hier=mom_hier
            )
            return mom_hier

    def getMHDHierarchiesAfterSteps(
        self,
        ndim,
        model_kwargs,
        diag_outputs,
        time_step,
        time_step_nbr,
        quantities=("B0", "B1"),
        interp_order=1,
        cells=None,
        dl=None,
    ):
        """Like getMHDHierarchies, but advances time_step_nbr steps and dumps the requested raw
        electromag quantities at the final time. Used to exercise the time-dependent B0(x,t) restamp
        and the -dB0/dt split sources."""
        if cells is None:
            cells = {1: 40, 2: 20, 3: 16}[ndim]
        if dl is None:
            dl = 1.0 / cells

        from pyphare.pharein.simulation import check_patch_size

        _, smallest_patch_size = check_patch_size(
            ndim, interp_order=interp_order, cells=cells
        )

        base_diag_dir = os.path.join("phare_outputs/init_mhd", diag_outputs)
        sim = self.simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=10,
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            boundary_types=["periodic"] * ndim,
            cells=phut.np_array_ify(cells, ndim),
            dl=phut.np_array_ify(dl, ndim),
            interp_order=interp_order,
            refinement_boxes={},
            diag_options={
                "format": "phareh5",
                "options": {"dir": base_diag_dir, "mode": "overwrite"},
            },
            strict=True,
            nesting_buffer=1,
            hyper_mode="spatial",
            eta=0.0,
            nu=0.0,
            gamma=5.0 / 3.0,
            reconstruction="Linear",
            limiter="VanLeer",
            riemann="Rusanov",
            mhd_timestepper="TVDRK2",
            hall=False,
            res=False,
            hyper_res=False,
            model_options=["MHDModel"],
            max_mhd_level=1,
        )
        diag_dir = sim.diag_options["options"]["dir"]
        final_time = time_step * time_step_nbr

        fns = dict(
            density=_const(1.0),
            vx=_const(0.0),
            vy=_const(0.0),
            vz=_const(0.0),
            p=_const(1.0),
        )
        fns.update(model_kwargs)
        ph.MHDModel(**fns)

        for q in quantities:
            ph.ElectromagDiagnostics(quantity=q, write_timestamps=[final_time])

        Simulator(sim).run().reset()

        # the diagnostics dump only at final_time, so a plain read returns that single time
        return {
            q: hierarchy_from(h5_filename=diag_dir + f"/EM_{q}.h5") for q in quantities
        }

    def getMHDHierarchies(
        self,
        ndim,
        model_kwargs,
        diag_outputs,
        quantities=("B0", "B1"),
        interp_order=1,
        cells=None,
        dl=None,
        **kwargs,
    ):
        """Build an MHD simulation from model_kwargs, dump the requested electromag quantities
        (raw B0 / B1), initialize (no time stepping), and return {qty: hierarchy}."""
        if cells is None:
            cells = {1: 40, 2: 20, 3: 16}[ndim]
        if dl is None:
            dl = 1.0 / cells  # domain length L = 1 per direction

        from pyphare.pharein.simulation import check_patch_size

        _, smallest_patch_size = check_patch_size(
            ndim, interp_order=interp_order, cells=cells
        )

        base_diag_dir = os.path.join("phare_outputs/init_mhd", diag_outputs)
        sim = self.simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=10,
            time_step_nbr=1,
            time_step=0.001,
            boundary_types=["periodic"] * ndim,
            cells=phut.np_array_ify(cells, ndim),
            dl=phut.np_array_ify(dl, ndim),
            interp_order=interp_order,
            refinement_boxes={},
            diag_options={
                "format": "phareh5",
                "options": {"dir": base_diag_dir, "mode": "overwrite"},
            },
            strict=True,
            nesting_buffer=1,
            hyper_mode="spatial",
            eta=0.0,
            nu=0.02,
            gamma=5.0 / 3.0,
            reconstruction="Linear",
            limiter="VanLeer",
            riemann="Rusanov",
            mhd_timestepper="TVDRK2",
            hall=False,
            res=False,
            hyper_res=False,
            model_options=["MHDModel"],
            max_mhd_level=1,
        )
        diag_dir = sim.diag_options["options"]["dir"]

        fns = dict(
            density=_const(1.0),
            vx=_const(0.0),
            vy=_const(0.0),
            vz=_const(0.0),
            p=_const(1.0),
        )
        fns.update(model_kwargs)
        ph.MHDModel(**fns)

        for q in quantities:
            ph.ElectromagDiagnostics(quantity=q, write_timestamps=[0.0])

        Simulator(sim).initialize().reset()

        return {
            q: hierarchy_from(h5_filename=diag_dir + f"/EM_{q}.h5") for q in quantities
        }

    def _slices(self, ghosts, dim):
        return tuple(
            slice(int(g), -int(g) if int(g) else None) for g in ghosts[:dim]
        )

    def _assert_mhd_field(self, hier, prefix, expected, dim, atol, strip_ghosts):
        """Compare each component of the prefix (B0/B1) field to the expected analytical function
        at its native staggered coordinates. When strip_ghosts, only the domain is compared (the
        potential curl is inexact in the outermost ghost layer)."""
        fns = [expected[i] if callable(expected[i]) else _const(expected[i]) for i in range(3)]
        comps = ["x", "y", "z"]
        for ilvl, level in hier.levels().items():
            self.assertEqual(ilvl, 0)
            for patch in level.patches:
                for ci, c in enumerate(comps):
                    pd = patch.patch_datas[f"{prefix}{c}"]
                    data = pd.dataset[:]
                    coords = [pd.x[:]]
                    if dim >= 2:
                        coords.append(pd.y[:])
                    if dim == 3:
                        coords.append(pd.z[:])
                    if strip_ghosts:
                        sl = self._slices(pd.ghosts_nbr, dim)
                        data = data[sl]
                        coords = [
                            coords[d][
                                slice(
                                    int(pd.ghosts_nbr[d]),
                                    -int(pd.ghosts_nbr[d])
                                    if int(pd.ghosts_nbr[d])
                                    else None,
                                )
                            ]
                            for d in range(dim)
                        ]
                    mesh = [a.flatten() for a in np.meshgrid(*coords, indexing="ij")]
                    exp = np.asarray(fns[ci](*mesh)) + np.zeros(data.size)
                    assert_fp_any_all_close(
                        data.flatten(), exp.flatten(), atol=atol, rtol=0
                    )

    def _assert_mhd_div_free(self, hier, prefix, dim, atol):
        for ilvl, level in hier.levels().items():
            for patch in level.patches:
                bx_pd = patch.patch_datas[f"{prefix}x"]
                by_pd = patch.patch_datas[f"{prefix}y"]
                bx = bx_pd.dataset[:]
                by = by_pd.dataset[:]
                if dim == 2:
                    divB = _divB2D(bx, by, bx_pd.x, by_pd.y)
                elif dim == 3:
                    bz_pd = patch.patch_datas[f"{prefix}z"]
                    divB = _divB3D(
                        bx, by, bz_pd.dataset[:], bx_pd.x, by_pd.y, bz_pd.z
                    )
                else:
                    raise ValueError("div-free test is 2D/3D only")
                sl = self._slices(bx_pd.ghosts_nbr, dim)
                self.assertLess(float(np.max(np.abs(divB[sl]))), atol)

    # -- entry points invoked by the per-dimension test files ---------------------------------

    def _test_mhd_split_is_as_provided(self, dim, name, model_kwargs, b0, b1, potential):
        hiers = self.getMHDHierarchies(
            dim, model_kwargs, diag_outputs=f"split/{dim}/{name}/{self.ddt_test_id()}"
        )
        # component init is exact everywhere; linear-potential curl is exact on the domain but
        # approximate in the outermost ghost, so strip ghosts for potential paths.
        atol = 1e-12 if potential else 1e-14
        self._assert_mhd_field(hiers["B0"], "B0", b0, dim, atol, strip_ghosts=potential)
        self._assert_mhd_field(hiers["B1"], "B1", b1, dim, atol, strip_ghosts=potential)

    def _test_mhd_potential_is_divergence_free(self, dim):
        L = 1.0
        (a0x, a0y, a0z) = tilted_potential(dim, L)
        (a1x, a1y, a1z) = tilted_potential(dim, L)
        model_kwargs = dict(
            a0x=a0x, a0y=a0y, a0z=a0z, a1x=a1x, a1y=a1y, a1z=a1z
        )
        hiers = self.getMHDHierarchies(
            dim, model_kwargs, diag_outputs=f"divfree/{dim}/{self.ddt_test_id()}"
        )
        self._assert_mhd_div_free(hiers["B0"], "B0", dim, atol=1e-11)
        self._assert_mhd_div_free(hiers["B1"], "B1", dim, atol=1e-11)

    def _test_mhd_time_dependent_uniform_B0(self, dim):
        """Manufactured solution for a time-dependent external field B0(x,t).

        Take a spatially UNIFORM B0(t) = B0_0 + t * rate (trivially curl-free) on a plasma at rest
        (V = 0) with periodic boundaries and uniform rho/P. Then J = curl(B1) = 0, E = 0, so the
        total field obeys dB/dt = -curl E = 0 and stays at B0_0. The split therefore predicts, at
        every time t:
            B0(t) = B0_0 + t * rate         (re-stamped once per step by the solver)
            B1(t) = B - B0 = -t * rate      (built up by the -dB0/dt induction source)
        The plasma stays at rest, so this isolates and checks exactly the new restamp + source path.
        """
        B0_0 = (0.3, 0.2, 0.1)
        rate = (0.1, -0.2, 0.05)
        time_step = 0.001
        time_step_nbr = 5
        tf = time_step * time_step_nbr

        def b0(i):
            return lambda *args: B0_0[i] + args[-1] * rate[i] + 0.0 * args[0]

        def db0(i):
            return lambda *args: rate[i] + 0.0 * args[0]

        model_kwargs = dict(
            b0x=b0(0), b0y=b0(1), b0z=b0(2),
            db0x_dt=db0(0), db0y_dt=db0(1), db0z_dt=db0(2),
        )
        hiers = self.getMHDHierarchiesAfterSteps(
            dim,
            model_kwargs,
            diag_outputs=f"tdep_b0/{dim}/{self.ddt_test_id()}",
            time_step=time_step,
            time_step_nbr=time_step_nbr,
        )
        expected_b0 = tuple(B0_0[i] + tf * rate[i] for i in range(3))
        expected_b1 = tuple(-tf * rate[i] for i in range(3))
        # uniform fields: exact to round-off, no ghost stripping needed
        self._assert_mhd_field(hiers["B0"], "B0", expected_b0, dim, atol=1e-12, strip_ghosts=False)
        self._assert_mhd_field(hiers["B1"], "B1", expected_b1, dim, atol=1e-12, strip_ghosts=False)

    def _test_mhd_time_dependent_potential_B0(self, dim):
        """Same manufactured solution as the component test, but B0 comes from a time-dependent
        vector potential A0(x,t) = t * A_lin(x). Since A_lin is linear, B0(t) = curl(A0) = t * B_lin
        is uniform, and dB0/dt = curl(dA0/dt) = curl(A_lin) = B_lin. This exercises the potential
        restamp path (B0 = curl(A0(t)), dB0/dt = curl(dA0/dt(t)) via the same discrete curl). With a
        plasma at rest the total field stays at B0(0) = 0, so B1(tf) = -tf * B_lin.
        """
        (ax, ay, az), B_lin = linear_potential({"xz": 0.2, "yx": 0.3, "zy": 0.1}, dim)
        a_fns = (ax, ay, az)
        time_step = 0.001
        time_step_nbr = 5
        tf = time_step * time_step_nbr

        def a0(i):
            return lambda *args: args[-1] * a_fns[i](*args[:-1]) + 0.0 * args[0]

        def da0(i):
            return lambda *args: a_fns[i](*args[:-1]) + 0.0 * args[0]

        model_kwargs = dict(
            a0x=a0(0), a0y=a0(1), a0z=a0(2),
            da0x_dt=da0(0), da0y_dt=da0(1), da0z_dt=da0(2),
        )
        hiers = self.getMHDHierarchiesAfterSteps(
            dim,
            model_kwargs,
            diag_outputs=f"tdep_b0_pot/{dim}/{self.ddt_test_id()}",
            time_step=time_step,
            time_step_nbr=time_step_nbr,
        )
        expected_b0 = tuple(tf * B_lin[i] for i in range(3))
        expected_b1 = tuple(-tf * B_lin[i] for i in range(3))
        # linear-potential curl is exact on the domain but approximate in the outermost ghost
        self._assert_mhd_field(hiers["B0"], "B0", expected_b0, dim, atol=1e-11, strip_ghosts=True)
        self._assert_mhd_field(hiers["B1"], "B1", expected_b1, dim, atol=1e-11, strip_ghosts=True)
