#

import os
import numpy as np


import pyphare.pharein as ph
from pyphare.core import phare_utilities as phut
from pyphare.core.phare_utilities import assert_fp_any_all_close
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from

from tests.diagnostic import all_timestamps
from tests.simulator.test_initialization import InitializationTest


def _const(c):
    return lambda *xyz: c + 0.0 * xyz[0]


def linear_potential(ndim):
    if ndim == 2:
        potential = (
            lambda x, y: 0.3 * y,
            lambda x, y: 0.2 * x,
            lambda x, y: 0.5 * x + 0.1 * y,
        )
        return potential, (0.1, -0.5, -0.1)
    potential = (
        lambda x, y, z: 0.3 * y + 0.6 * z,
        lambda x, y, z: 0.2 * x + 0.4 * z,
        lambda x, y, z: 0.5 * x + 0.1 * y,
    )
    return potential, (-0.3, 0.1, -0.1)


def linear_potential_growing_in_time(ndim, power=1):
    (ax, ay, az), b = linear_potential(ndim)

    def grows(a):
        if ndim == 2:
            return lambda x, y, t: t**power * a(x, y)
        return lambda x, y, z, t: t**power * a(x, y, z)

    def rate(a):
        if ndim == 2:
            return lambda x, y, t: power * t ** (power - 1) * a(x, y)
        return lambda x, y, z, t: power * t ** (power - 1) * a(x, y, z)

    potential = tuple(grows(a) for a in (ax, ay, az))
    derivative = tuple(rate(a) for a in (ax, ay, az))
    return potential, derivative, b


def mhd_split_cases(ndim):
    potential, b0 = linear_potential(ndim)
    external_field = {"type": "user-defined", "potential": potential}
    total = (0.3, 0.2, 0.1)
    b1 = (0.05, -0.1, 0.2)
    zero = (0.0, 0.0, 0.0)
    total_kwargs = dict(bx=_const(total[0]), by=_const(total[1]), bz=_const(total[2]))
    b1_kwargs = dict(b1x=_const(b1[0]), b1y=_const(b1[1]), b1z=_const(b1[2]))
    return [
        dict(
            name="default",
            model_kwargs={},
            external_field=None,
            b0=zero,
            b1=(1.0, 0.0, 0.0),
        ),
        dict(
            name="total_without_b0",
            model_kwargs=total_kwargs,
            external_field=None,
            b0=zero,
            b1=total,
        ),
        dict(
            name="total_with_b0",
            model_kwargs=total_kwargs,
            external_field=external_field,
            b0=b0,
            b1=tuple(t - b for t, b in zip(total, b0)),
        ),
        dict(
            name="b1_without_b0",
            model_kwargs=b1_kwargs,
            external_field=None,
            b0=zero,
            b1=b1,
        ),
        dict(
            name="b1_with_b0",
            model_kwargs=b1_kwargs,
            external_field=external_field,
            b0=b0,
            b1=b1,
        ),
    ]


class MHDInitializationTest(InitializationTest):
    def getHierarchy(
        self,
        ndim,
        _interp_order,
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
        extra_diag_options=None,
        timestamps=None,
        diag_outputs="",
        max_mhd_level=1,
        **kwargs,
    ):
        # the shared hybrid+MHD base passes this positionally: an MHD-only run
        # has no interp order, so it is accepted and deliberately unused --
        # forwarding it would override pharein's derivation from model_options
        # and select the coupled hybrid+MHD module.
        if smallest_patch_size is None:
            from pyphare.pharein.simulation import check_patch_size

            _, smallest_patch_size = check_patch_size(
                ndim, interp_order=None, cells=cells
            )

        # ----------------------------------------------------------------------
        # simulation setup and running
        # ----------------------------------------------------------------------
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
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5", "options": extra_diag_options},
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
            hall=hall,
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
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_B.h5", hier=eb_hier
            )
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

    def getMHDSplitHierarchies(
        self,
        ndim,
        model_kwargs,
        external_field=None,
        time_step_nbr=1,
        time_step=0.001,
        run=False,
        quantities=("B0", "B1", "B"),
        mhd_timestepper="TVDRK2",
    ):
        cells = {2: 20, 3: 16}[ndim]

        from pyphare.pharein.simulation import check_patch_size

        _, smallest_patch_size = check_patch_size(ndim, interp_order=None, cells=cells)

        sim_kwargs = dict(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=10,
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            boundary_types=["periodic"] * ndim,
            cells=phut.np_array_ify(cells, ndim),
            dl=phut.np_array_ify(1.0 / cells, ndim),
            diag_options={
                "format": "phareh5",
                "options": {"dir": "phare_outputs/init_mhd_split", "mode": "overwrite"},
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
            mhd_timestepper=mhd_timestepper,
            hall=False,
            model_options=["MHDModel"],
            max_mhd_level=1,
        )
        if external_field is not None:
            sim_kwargs["external_field"] = external_field
        sim = self.simulation(**sim_kwargs)
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

        dump_time = time_step * time_step_nbr if run else 0.0
        for q in quantities:
            ph.ElectromagDiagnostics(quantity=q, write_timestamps=[dump_time])

        simulator = Simulator(sim)
        if run:
            simulator.run()
        else:
            simulator.initialize()
        simulator.reset()

        return {
            q: hierarchy_from(h5_filename=diag_dir + f"/EM_{q}.h5") for q in quantities
        }

    def _assert_uniform_field(self, hier, prefix, expected, atol):
        nodes = 0
        for level in hier.levels().values():
            for patch in level.patches:
                for c, value in zip("xyz", expected):
                    data = patch.patch_datas[f"{prefix}{c}"].dataset[:]
                    assert_fp_any_all_close(
                        data.flatten(), np.full(data.size, value), atol=atol, rtol=0
                    )
                    nodes += data.size
        self.assertGreater(nodes, 0)

    def _test_mhd_split_is_as_provided(
        self, ndim, name, model_kwargs, external_field, b0, b1
    ):
        hiers = self.getMHDSplitHierarchies(ndim, model_kwargs, external_field)
        self._assert_uniform_field(hiers["B0"], "B0", b0, atol=1e-13)
        self._assert_uniform_field(hiers["B1"], "B1", b1, atol=1e-13)
        total = tuple(x + y for x, y in zip(b0, b1))
        self._assert_uniform_field(hiers["B"], "B", total, atol=1e-13)

    def _test_mhd_time_dependent_uniform_B0(self, ndim):
        potential, derivative, b_lin = linear_potential_growing_in_time(ndim)
        time_step, time_step_nbr = 0.001, 5
        final_time = time_step * time_step_nbr
        hiers = self.getMHDSplitHierarchies(
            ndim,
            dict(b1x=_const(0.0), b1y=_const(0.0), b1z=_const(0.0)),
            external_field={
                "type": "user-defined",
                "potential": potential,
                "potential_time_derivative": derivative,
            },
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            run=True,
        )
        self._assert_uniform_field(
            hiers["B0"], "B0", tuple(final_time * b for b in b_lin), atol=1e-12
        )
        self._assert_uniform_field(
            hiers["B1"], "B1", tuple(-final_time * b for b in b_lin), atol=1e-12
        )
        self._assert_uniform_field(hiers["B"], "B", (0.0, 0.0, 0.0), atol=1e-12)

    def _test_mhd_B0_follows_the_rk_stages(self, ndim, mhd_timestepper):
        potential, derivative, b_lin = linear_potential_growing_in_time(ndim, power=2)
        time_step, time_step_nbr = 0.001, 5
        final_time = time_step * time_step_nbr
        hiers = self.getMHDSplitHierarchies(
            ndim,
            dict(b1x=_const(0.0), b1y=_const(0.0), b1z=_const(0.0)),
            external_field={
                "type": "user-defined",
                "potential": potential,
                "potential_time_derivative": derivative,
            },
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            run=True,
            mhd_timestepper=mhd_timestepper,
        )
        b0 = tuple(final_time**2 * b for b in b_lin)
        self._assert_uniform_field(hiers["B0"], "B0", b0, atol=1e-14)
        self._assert_uniform_field(hiers["B1"], "B1", tuple(-b for b in b0), atol=1e-14)
        self._assert_uniform_field(hiers["B"], "B", (0.0, 0.0, 0.0), atol=1e-14)
