
import copy

import unittest
import numpy as np

from ddt import ddt, data, unpack

from pyphare.cpp import cpp_lib
cpp = cpp_lib()

import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest
from tests.diagnostic import dump_all_diags


def permute(dict):
    #from pyphare.pharein.simulation import supported_dimensions # eventually
    dims = [1] # supported_dimensions()
    return [
      [dim, interp, dict] for dim in dims for interp in [1,2,3]
    ]


def setup_model(ppc=100):
    def density(x): return 1.
    def bx(x): return 0.
    def S(x,x0,l): return 0.5*(1+np.tanh((x-x0)/l))
    def by(x):
        L = ph.global_vars.sim.simulation_domain()[0]
        v1, v2= -1, 1.
        return v1 + (v2-v1)*(S(x,L*0.25,1) -S(x, L*0.75, 1))
    def bz(x): return 0.5
    def b2(x): return bx(x)**2 + by(x)**2 + bz(x)**2
    def T(x):
        K = 1
        return 1/density(x)*(K - b2(x)*0.5)
    def vx(x): return 2.
    def vy(x): return 0.
    def vz(x): return 0.
    def vthxyz(x): return T(x)
    vvv = {
        "vbulkx": vx, "vbulky": vy, "vbulkz": vz,
        "vthx": vthxyz, "vthy": vthxyz, "vthz": vthxyz
    }
    model = ph.MaxwellianFluidModel(
        bx=bx, by=by, bz=bz,
        protons={"mass":1, "charge": 1, "density": density, **vvv, "nbr_part_per_cell":ppc, "init": {"seed": 1337}},
        alpha={"mass":4, "charge": 1, "density": density, **vvv, "nbr_part_per_cell":ppc, "init": {"seed": 2334}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    return model



timestep=.001
out = "phare_outputs/restarts/"
simArgs = dict(
  time_step_nbr = 5, # avoid regrid for refinement boxes https://github.com/LLNL/SAMRAI/issues/199
  time_step = timestep,
  boundary_types = "periodic",
  cells = 200,
  dl = 0.3,
  diag_options = dict(format="phareh5", options=dict(dir=out, mode="overwrite")),
  restart_options = dict(dir=out, mode="overwrite", timestamps=[timestep*4])
)

def dup(dic):
    dic.update(copy.deepcopy(simArgs))
    return dic


@ddt
class RestartsTest(SimulatorTest):


    def __init__(self, *args, **kwargs):
        super(RestartsTest, self).__init__(*args, **kwargs)
        self.simulator = None


    def tearDown(self):
        super(RestartsTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None


    def ddt_test_id(self):
        return self._testMethodName.split("_")[-1]


    @data(
      *permute(dup(dict(
          max_nbr_levels=2,
          refinement="tagging",
      ))),
      *permute(dup(dict())), # refinement boxes set later
    )
    @unpack
    def test_restarts(self, dim, interp, simInput):
        print(f"test_restarts dim/interp:{dim}/{interp}")

        simput = copy.deepcopy(simInput)

        for key in ["cells", "dl", "boundary_types"]:
            simput[key] = [simput[key]] * dim

        if "refinement" not in simput:
            b0 = [[10 for i in range(dim)], [19 for i in range(dim)]]
            simput["refinement_boxes"] = {"L0": {"B0": b0}}
        else: # https://github.com/LLNL/SAMRAI/issues/199
            simput["time_step_nbr"] = 10

        # if restart time exists it "loads" from restart file
        #  otherwise just saves restart files based on timestamps
        assert "restart_time" not in simput["restart_options"]

        simput["interp_order"] = interp
        time_step = simput["time_step"]
        time_step_nbr = simput["time_step_nbr"]

        restart_idx = 4
        restart_time=time_step * restart_idx
        timestamps = [time_step * restart_idx, time_step * time_step_nbr]

        # first simulation
        local_out = f"{out}/test/{dim}/{interp}/mpi_n/{cpp.mpi_size()}/id{self.ddt_test_id()}"
        simput["restart_options"]["dir"] = local_out
        simput["diag_options"]["options"]["dir"] = local_out
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        assert "restart_time" not in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))
        Simulator(ph.global_vars.sim).run().reset()
        self.register_diag_dir_for_cleanup(local_out)
        diag_dir0 = local_out


        # second restarted simulation
        local_out = f"{local_out}_n2"
        simput["diag_options"]["options"]["dir"] = local_out
        simput["restart_options"]["restart_time"] = restart_time
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        assert "restart_time" in ph.global_vars.sim.restart_options
        model = setup_model()
        dump_all_diags(model.populations, timestamps=np.array(timestamps))
        Simulator(ph.global_vars.sim).run().reset()
        self.register_diag_dir_for_cleanup(local_out)
        diag_dir1 = local_out



        def check(qty0, qty1, checker):
            checks = 0
            for ilvl, lvl0 in qty0.patch_levels.items():
                patch_level1 = qty1.patch_levels[ilvl]
                for p_idx, patch0 in enumerate(lvl0):
                    patch1 = patch_level1.patches[p_idx]
                    for pd_key, pd0 in patch0.patch_datas.items():
                        pd1 = patch1.patch_datas[pd_key]
                        self.assertNotEqual(id(pd0), id(pd1))
                        checker(pd0, pd1)
                        checks += 1
            return checks

        def check_particles(qty0, qty1):
            return check(qty0, qty1, lambda pd0, pd1: self.assertEqual(pd0.dataset, pd1.dataset))

        def check_field(qty0, qty1):
            return  check(qty0, qty1, lambda pd0, pd1: np.testing.assert_equal(pd0.dataset[:], pd1.dataset[:]))

        pops = model.populations
        for time in timestamps:
            checks = 0

            run0 = Run(diag_dir0)
            run1 = Run(diag_dir1)
            checks += check_particles(run0.GetParticles(time, pops), run1.GetParticles(time, pops))
            checks += check_field(run0.GetB(time), run1.GetB(time))
            checks += check_field(run0.GetE(time), run1.GetE(time))
            checks += check_field(run0.GetNi(time), run1.GetNi(time))
            checks += check_field(run0.GetVi(time), run1.GetVi(time))

            for pop in pops:
                checks += check_field(run0.GetFlux(time, pop), run1.GetFlux(time, pop))
                checks += check_field(run0.GetN(time, pop), run1.GetN(time, pop))

            self.assertGreaterEqual(checks, 14)






    def test_mode_conserve(self, dim = 1, interp = 1 , simput = dup(simArgs)):
        print(f"test_mode_conserve dim/interp:{dim}/{interp}")

        for key in ["cells", "dl", "boundary_types"]:
            simput[key] = [simput[key]] * dim

        # first simulation
        local_out = f"{out}/conserve/{dim}/{interp}/mpi_n/{cpp.mpi_size()}/id{self.ddt_test_id()}"
        self.register_diag_dir_for_cleanup(local_out)

        simput["restart_options"]["dir"] = local_out
        simput["restart_options"]["mode"] = "conserve"
        ph.global_vars.sim = ph.Simulation(**simput)
        self.assertEqual(len(ph.global_vars.sim.restart_options["timestamps"]), 1)
        self.assertEqual(ph.global_vars.sim.restart_options["timestamps"][0], .004)
        model = setup_model()
        Simulator(ph.global_vars.sim).run().reset()

        # second simulation (not restarted)
        ph.global_vars.sim = None
        ph.global_vars.sim = ph.Simulation(**simput)
        self.assertEqual(len(ph.global_vars.sim.restart_options["timestamps"]), 0)







if __name__ == "__main__":
    unittest.main()


