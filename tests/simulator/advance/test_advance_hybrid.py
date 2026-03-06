#
#

import os
import numpy as np

import pyphare.pharein as ph
import pyphare.core.box as boxm
from pyphare.core import phare_utilities as phut
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.pharesee.geometry import hierarchy_overlaps
from pyphare.pharesee.hierarchy.hierarchy_utils import merge_particles

from tests.diagnostic import all_timestamps
from tests.simulator.test_advance import AdvanceTestBase


class HybridAdvanceTest(AdvanceTestBase):
    def getHierarchy(
        self,
        ndim,
        interp_order,
        qty,
        refinement_boxes={},
        nbr_part_per_cell=100,
        density=None,
        smallest_patch_size=None,
        largest_patch_size=20,
        cells=120,
        time_step=0.001,
        model_init=None,
        dl=0.2,
        extra_diag_options=None,
        time_step_nbr=1,
        timestamps=None,
        block_merging_particles=False,
        diag_outputs="",
    ):
        if smallest_patch_size is None:
            from pyphare.pharein.simulation import check_patch_size

            _, smallest_patch_size = check_patch_size(
                ndim, interp_order=interp_order, cells=cells
            )

        model_init = model_init or dict()

        base_diag_dir = "phare_outputs/advance"
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

        def vth(*xyz):
            return 0.01 + np.zeros_like(xyz[0])

        def vthx(*xyz):
            return vth(*xyz)

        def vthy(*xyz):
            return vth(*xyz)

        def vthz(*xyz):
            return vth(*xyz)

        ph.MaxwellianFluidModel(
            bx=bx,
            by=by,
            bz=bz,
            protons={
                "charge": 1,
                "density": density or _density,
                "vbulkx": vx,
                "vbulky": vy,
                "vbulkz": vz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
                "nbr_part_per_cell": nbr_part_per_cell,
                "init": model_init,
            },
        )

        ph.ElectronModel(closure="isothermal", Te=0.12)

        if timestamps is None:
            timestamps = all_timestamps(ph.global_vars.sim)

        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

        for quantity in ["charge_density", "bulkVelocity"]:
            ph.FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
            )

        poplist = ["protons"]
        for pop in poplist:
            for quantity in ["density", "flux"]:
                ph.FluidDiagnostics(
                    quantity=quantity,
                    write_timestamps=timestamps,
                    population_name=pop,
                )

            for quantity in ["domain", "levelGhost"]:
                ph.ParticleDiagnostics(
                    quantity=quantity, write_timestamps=timestamps, population_name=pop
                )

        Simulator(sim).run()

        eb_hier = None
        if qty in ["e", "eb", "fields"]:
            eb_hier = hierarchy_from(
                h5_filename=diag_outputs + "/EM_E.h5", hier=eb_hier
            )
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

        # else particle tests

        assert qty == "particles"

        particle_hier = hierarchy_from(
            h5_filename=diag_outputs + "/ions_pop_protons_domain.h5"
        )
        particle_hier = hierarchy_from(
            h5_filename=diag_outputs + "/ions_pop_protons_levelGhost.h5",
            hier=particle_hier,
        )

        if not block_merging_particles:
            merge_particles(particle_hier)
        return particle_hier

    def _test_overlapped_particledatas_have_identical_particles(
        self, ndim, interp_order, refinement_boxes, ppc=100, **kwargs
    ):
        print(
            "test_overlapped_particledatas_have_identical_particles, interporder : {}".format(
                interp_order
            )
        )
        from copy import copy

        time_step_nbr = 3
        time_step = 0.001

        datahier = self.getHierarchy(
            ndim,
            interp_order,
            qty="particles",
            refinement_boxes=refinement_boxes,
            time_step=time_step,
            time_step_nbr=time_step_nbr,
            nbr_part_per_cell=ppc,
            **kwargs,
        )

        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time = time_step_idx * time_step

            overlaps = hierarchy_overlaps(datahier, coarsest_time)

            for ilvl in datahier.levels():
                print("testing level {}".format(ilvl))
                for overlap in overlaps[ilvl]:
                    pd1, pd2 = overlap["pdatas"]
                    box = overlap["box"]
                    offsets = overlap["offset"]

                    self.assertEqual(pd1.quantity, pd2.quantity)

                    if "particles" in pd1.quantity:
                        # the following uses 'offset', we need to remember that offset
                        # is the quantity by which a patch has been moved to detect
                        # overlap with the other one.
                        # so shift by +offset when evaluating patch data in overlap box
                        # index space, and by -offset when we want to shift box indexes
                        # to the associated patch index space.

                        # overlap box must be shifted by -offset to select data in the patches
                        part1 = copy(
                            pd1.dataset.select(boxm.shift(box, -np.asarray(offsets[0])))
                        )
                        part2 = copy(
                            pd2.dataset.select(boxm.shift(box, -np.asarray(offsets[1])))
                        )

                        # periodic icell overlaps need shifting to be the same
                        part1.iCells = part1.iCells + offsets[0]
                        part2.iCells = part2.iCells + offsets[1]
                        self.assertEqual(part1, part2)

    def _test_L0_particle_number_conservation(self, ndim, interp_order, ppc=100):
        cells = 120
        time_step_nbr = 10
        time_step = 0.001

        n_particles = ppc * (cells**ndim)

        datahier = self.getHierarchy(
            ndim,
            interp_order,
            "particles",
            None,
            time_step=time_step,
            time_step_nbr=time_step_nbr,
            nbr_part_per_cell=ppc,
            cells=cells,
        )

        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time = time_step_idx * time_step
            n_particles_at_t = 0
            for patch in datahier.level(0, coarsest_time).patches:
                n_particles_at_t += (
                    patch.patch_datas["protons_particles"].dataset[patch.box].size()
                )
            self.assertEqual(n_particles, n_particles_at_t)

    def base_test_domain_particles_on_refined_level(self, datahier, new_time=None):
        """
        !! test assumes only domain particle patch_datas are present !!
        """
        times = datahier.times() if new_time is None else [new_time]
        checks = 0
        for coarsest_time in times:
            for patch in datahier.level(1, coarsest_time).patches:
                for pd_key, pd in patch.patch_datas.items():
                    if pd_key.endswith("_domain"):
                        self.assertGreater(pd[pd.box].size(), 0)
                        self.assertEqual(pd[pd.box].size(), pd[pd.ghost_box].size())
                        checks += 1
        self.assertGreater(checks, 0)

    def _test_domain_particles_on_refined_level(
        self, ndim, interp_order, refinement_boxes, **kwargs
    ):
        time_step_nbr = 5
        time_step = 0.001

        self.base_test_domain_particles_on_refined_level(
            self.getHierarchy(
                ndim,
                interp_order,
                qty="particles",
                refinement_boxes=refinement_boxes,
                time_step=time_step,
                time_step_nbr=time_step_nbr,
                block_merging_particles=True,
                **kwargs,
            )
        )
