#

import os
import numpy as np


import pyphare.pharein as ph
from pyphare.core import phare_utilities as phut
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from

from tests.diagnostic import all_timestamps
from tests.simulator.test_initialization import InitializationTest


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
        hall=True,
        res=False,
        hyper_res=True,
        extra_diag_options=None,
        timestamps=None,
        diag_outputs="",
        **kwargs
    ):
        if smallest_patch_size is None:
            from pyphare.pharein.simulation import check_patch_size

            _, smallest_patch_size = check_patch_size(
                ndim, interp_order=interp_order, cells=cells
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
            interp_order=interp_order,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5", "options": extra_diag_options},
            strict=True,
            nesting_buffer=1,
            hyper_mode="spatial",
            eta=0.0,
            nu=0.02,
            gamma=5.0 / 3.0,
            reconstruction="WENOZ",
            limiter="None",
            riemann="Rusanov",
            mhd_timestepper="TVDRK3",
            hall=hall,
            res=res,
            hyper_res=hyper_res,
            model_options=["MHDModel"],
            max_mhd_level=3,
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
