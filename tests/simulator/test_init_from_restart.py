import sys
import copy
import unittest
import subprocess
import numpy as np
import pyphare.pharein as ph

from pyphare.simulator.simulator import Simulator
from pyphare.core import phare_utilities as phut
from pyphare.pharesee.hierarchy.fromh5 import get_all_available_quantities_from_h5
from pyphare.pharesee.particles import single_patch_per_level_per_pop_from
from pyphare.pharesee.hierarchy.hierarchy_utils import flat_finest_field

from tests.simulator import SimulatorTest, test_restarts
from tests.diagnostic import dump_all_diags

timestep = 0.001
first_mpi_size = 4
ppc = 100
cells = 200
first_out = "phare_outputs/reinit/first"
secnd_out = "phare_outputs/reinit/secnd"
timestamps = np.array([timestep * 4])
simInitArgs = dict(
    largest_patch_size=100,
    time_step_nbr=5,
    time_step=timestep,
    cells=cells,
    dl=0.3,
    init_options=dict(dir=f"{first_out}/00000.00400", mpi_size=first_mpi_size),
    diag_options=dict(format="phareh5", options=dict(dir=secnd_out, mode="overwrite")),
)


def setup_model(sim):
    model = ph.MaxwellianFluidModel(
        protons={"mass": 1, "charge": 1, "nbr_part_per_cell": ppc},
        alpha={"mass": 4.0, "charge": 1, "nbr_part_per_cell": ppc},
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    dump_all_diags(model.populations, timestamps)
    return model


class RestartsParserTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(RestartsParserTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(RestartsParserTest, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def test_reinit(self):
        self.register_diag_dir_for_cleanup("phare_outputs/reinit")
        sim = ph.Simulation(**copy.deepcopy(simInitArgs))
        setup_model(sim)
        Simulator(sim).run().reset()
        datahier0 = get_all_available_quantities_from_h5(first_out, timestamps[0])
        datahier1 = get_all_available_quantities_from_h5(secnd_out, timestamps[0])

        for k in "xyz":
            a = flat_finest_field(datahier0, f"B{k}", timestamps[0], 0)
            b = flat_finest_field(datahier1, f"B{k}", timestamps[0], 0)
            phut.assert_fp_any_all_close(a, b)

        def get_merged(hier):
            return single_patch_per_level_per_pop_from(datahier0)

        ds = [get_merged(datahier0), get_merged(datahier1)]
        for key in ["alpha", "protons"]:
            a, b = [d.level(0).patches[0].patch_datas[f"{key}_domain"] for d in ds]
            self.assertGreater(a.size(), (cells - 1) * ppc)
            self.assertEqual(a, b)


def run_first_sim():
    simput = copy.deepcopy(test_restarts.simArgs)
    simput["restart_options"]["dir"] = first_out
    simput["restart_options"]["timestamps"] = timestamps
    simput["diag_options"]["options"]["dir"] = first_out
    sim = ph.Simulation(**simput)
    dump_all_diags(test_restarts.setup_model().populations, timestamps=timestamps)
    Simulator(sim).run()


def launch():
    cmd = f"mpirun -n {first_mpi_size} python3 -O tests/simulator/test_init_from_restart.py lol"
    proc = subprocess.run(cmd.split(" "), check=True, capture_output=True)
    return (proc.stdout + proc.stderr).decode().strip()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        launch()
        unittest.main()
    else:
        run_first_sim()
