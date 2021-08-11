from pyphare.cpp import cpp_lib
cpp = cpp_lib()

from pyphare.simulator.simulator import Simulator, startMPI
from pyphare.core.phare_utilities import np_array_ify
from pyphare.pharesee.hierarchy import hierarchy_from, merge_particles
from pyphare.pharein import MaxwellianFluidModel
from pyphare.pharein.diagnostics import ParticleDiagnostics, FluidDiagnostics, ElectromagDiagnostics
from pyphare.pharein import ElectronModel
from pyphare.pharein.simulation import Simulation, supported_dimensions
from pyphare.pharesee.geometry import level_ghost_boxes, hierarchy_overlaps
from pyphare.core.gridlayout import yee_element_is_primal
from pyphare.pharesee.particles import aggregate as aggregate_particles
import pyphare.core.box as boxm
from pyphare.core.box import Box, Box1D
import numpy as np
import unittest
from ddt import ddt, data, unpack
from tests.diagnostic import all_timestamps
from tests.simulator.test_advance import AdvanceTestBase

interp_orders = [1]
def per_interp(dic):
    return [(interp, dic) for interp in interp_orders]

@ddt
class OffloadTest(unittest.TestCase):

    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]


    def _density(*xyz):
        from pyphare.pharein.global_vars import sim
        hL = np.array(sim.simulation_domain()) / 2
        _ = lambda i: -(xyz[i]-hL[i]) ** 2
        return .3 + np.exp(sum([_(i) for i,v in enumerate(xyz)]))



    def getHierarchy(self, interp_order, refinement_boxes, qty,
                     diag_outputs, nbr_part_per_cell=100, density = _density,
                     smallest_patch_size=None, largest_patch_size=20,
                     cells=120, time_step=0.001, model_init={},
                     dl=0.2, extra_diag_options={}, time_step_nbr=1, timestamps=None, ndim=1):
        diag_outputs = f"phare_outputs/offloading/{diag_outputs}"
        from pyphare.pharein import global_vars
        global_vars.sim = None

        if smallest_patch_size is None:
            from pyphare.pharein.simulation import check_patch_size
            _, smallest_patch_size = check_patch_size(ndim, interp_order=interp_order, cells=cells)

        startMPI()
        extra_diag_options["mode"] = "overwrite"
        extra_diag_options["dir"] = diag_outputs
        Simulation(
            smallest_patch_size=smallest_patch_size,
            largest_patch_size=largest_patch_size,
            time_step_nbr=time_step_nbr,
            time_step=time_step,
            boundary_types=["periodic"] * ndim,
            cells=np_array_ify(cells, ndim),
            dl=np_array_ify(dl, ndim),
            interp_order=interp_order,
            refinement_boxes=refinement_boxes,
            diag_options={"format": "phareh5",
                          "options": extra_diag_options},
            strict=True,
            offload='cuda',
        )


        def S(x,x0,l):
            return 0.5*(1+np.tanh((x-x0)/l))

        def bx(*xyz):
            return 1.


        def by(*xyz):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            _ = lambda i: 0.1*np.cos(2*np.pi*xyz[i]/L[i])
            return np.asarray([_(i) for i,v in enumerate(xyz)]).prod(axis=0)

        def bz(*xyz):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            _ = lambda i: 0.1*np.cos(2*np.pi*xyz[i]/L[i])
            return np.asarray([_(i) for i,v in enumerate(xyz)]).prod(axis=0)

        def vx(*xyz):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            _ = lambda i: 0.1*np.cos(2*np.pi*xyz[i]/L[i])
            return np.asarray([_(i) for i,v in enumerate(xyz)]).prod(axis=0)

        def vy(*xyz):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            _ = lambda i: 0.1*np.cos(2*np.pi*xyz[i]/L[i])
            return np.asarray([_(i) for i,v in enumerate(xyz)]).prod(axis=0)

        def vz(*xyz):
            from pyphare.pharein.global_vars import sim
            L = sim.simulation_domain()
            _ = lambda i: 0.1*np.cos(2*np.pi*xyz[i]/L[i])
            return np.asarray([_(i) for i,v in enumerate(xyz)]).prod(axis=0)


        def vth(*xyz):
            return 0.01 + np.zeros_like(xyz[0])

        def vthx(*xyz):
            return vth(*xyz)

        def vthy(*xyz):
            return vth(*xyz)

        def vthz(*xyz):
            return vth(*xyz)


        MaxwellianFluidModel(bx=bx, by=by, bz=bz,
                             protons={"charge": 1,
                                      "density": density,
                                      "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
                                      "vthx": vthx, "vthy": vthy, "vthz": vthz,
                                      "nbr_part_per_cell": nbr_part_per_cell,
                                      "init": model_init,
                                      })

        ElectronModel(closure="isothermal", Te=0)#0.12)

        if timestamps is None:
            timestamps = all_timestamps(global_vars.sim)

        for quantity in ["E", "B"]:
            ElectromagDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps
            )

        for quantity in ["density", "bulkVelocity"]:
            FluidDiagnostics(
                quantity=quantity,
                write_timestamps=timestamps,
                compute_timestamps=timestamps
            )

        poplist = ["protons"]
        for pop in poplist:
            for quantity in ["density", "flux"]:
                FluidDiagnostics(quantity=quantity,
                                 write_timestamps=timestamps,
                                 compute_timestamps=timestamps,
                                 population_name=pop)

            for quantity in ['domain', 'levelGhost', 'patchGhost']:
                ParticleDiagnostics(quantity=quantity,
                                    compute_timestamps=timestamps,
                                    write_timestamps=timestamps,
                                    population_name=pop)

        Simulator(global_vars.sim).run()

        eb_hier = None
        if qty in ["e", "eb", "fields"]:
            eb_hier = hierarchy_from(h5_filename=diag_outputs+"/EM_E.h5", hier=eb_hier)
        if qty in ["b", "eb", "fields"]:
            eb_hier = hierarchy_from(h5_filename=diag_outputs+"/EM_B.h5", hier=eb_hier)
        if qty in ["e", "b", "eb"]:
            return eb_hier

        is_particle_type = qty == "particles" or qty == "particles_patch_ghost"

        if is_particle_type:
            particle_hier = None

        if qty == "particles":
            particle_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_domain.h5")
            particle_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_levelGhost.h5", hier=particle_hier)

        if is_particle_type:
            particle_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_patchGhost.h5", hier=particle_hier)

        if qty == "particles":
            merge_particles(particle_hier)

        if is_particle_type:
            return particle_hier

        if qty == "moments" or qty == "fields":
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_density.h5", hier=eb_hier)
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_bulkVelocity.h5", hier=mom_hier)
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_density.h5", hier=mom_hier)
            mom_hier = hierarchy_from(h5_filename=diag_outputs+"/ions_pop_protons_flux.h5", hier=mom_hier)
            return mom_hier



    @data(
      *per_interp({}),
    )
    @unpack
    def test_offloading(self, interp_order, refinement_boxes):
        ndim = 1
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr=1
        time_step=0.001
        diag_outputs=f"offloader_comparison_{ndim}_{self.ddt_test_id()}"
        
        import random
        rando = random.randint(0, 1e10)
        def _getHier(getter, qty="eb"):
            return getter.getHierarchy(interp_order, refinement_boxes, qty, diag_outputs=diag_outputs,
                                  time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim,
                                  model_init={"seed": rando},)
        
        gpu_datahier = _getHier(self)
        cpu_datahier = _getHier(AdvanceTestBase())        
        
        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time =  time_step_idx * time_step
            cpu_patches = cpu_datahier.level(0, coarsest_time).patches
            gpu_patches = gpu_datahier.level(0, coarsest_time).patches
            for EM in ["E", "B"]:
                for xyz in ["x", "y", "z"]:
                    qty = f"{EM}{xyz}"
                    print(f"{qty} cpu_patches", cpu_patches[0].box, cpu_patches[0].patch_datas[qty].dataset[:])
                    print(f"{qty} gpu_patches", gpu_patches[0].box, gpu_patches[0].patch_datas[qty].dataset[:])

        cpu_datahier = _getHier(AdvanceTestBase(), "moments")
        gpu_datahier = _getHier(self, "moments")
        
        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time =  time_step_idx * time_step
            cpu_patches = cpu_datahier.level(0, coarsest_time).patches
            gpu_patches = gpu_datahier.level(0, coarsest_time).patches
            for qty in list(gpu_patches[0].patch_datas.keys()):
                print(qty, cpu_patches[0].box, cpu_patches[0].patch_datas[qty].dataset[:])
                print(qty, gpu_patches[0].box, gpu_patches[0].patch_datas[qty].dataset[:])
        
        #self._test_overlaped_fields_are_equal(datahier, time_step_nbr, time_step)
        qty = "particles"
        cpu_datahier = _getHier(AdvanceTestBase(), qty)
        gpu_datahier = _getHier(self, qty)
        
        def print_particle(patches, key, idx=0):
            particles = patches[0].patch_datas[pop_key].dataset[patches[0].box]
            print(f"{qty} {key}", patches[0].box, particles.v[idx],
                                                  particles.iCells[idx],
                                                  particles.deltas[idx],
                                                  particles.charges[idx], 
                                                  particles.weights[idx])
        
        for time_step_idx in range(time_step_nbr + 1):
            coarsest_time =  time_step_idx * time_step
            cpu_patches = cpu_datahier.level(0, coarsest_time).patches
            gpu_patches = gpu_datahier.level(0, coarsest_time).patches
            pop_key = list(gpu_patches[0].patch_datas.keys())[0]
            gpu_particles = gpu_patches[0].patch_datas[pop_key].dataset
            print_particle(gpu_patches, "gpu_patches")
            print_particle(cpu_patches, "cpu_patches")
            print_particle(gpu_patches, "gpu_patches", -1)
            print_particle(cpu_patches, "cpu_patches", -1)


if __name__ == "__main__":
    unittest.main()
