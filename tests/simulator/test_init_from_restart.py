import sys
import copy
import unittest
import subprocess
import numpy as np
import pyphare.pharein as ph
from pathlib import Path

from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy.patchdata import FieldData, ParticleData
from pyphare.pharesee.hierarchy.fromh5 import get_all_available_quantities_from_h5
from pyphare.pharesee.hierarchy.hierarchy import format_timestamp
from pyphare.pharesee.hierarchy.hierarchy_utils import single_patch_for_LO
from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare
from tests.simulator import SimulatorTest, test_restarts
from tests.diagnostic import dump_all_diags


time_step = 0.001
time_step_nbr = 5
final_time = time_step_nbr * time_step
first_mpi_size = 5
ppc = 100
cells = 200
first_out = "test_init_from_restart"
secnd_out = "phare_outputs/reinit/secnd"
timestamps = np.arange(0, final_time + time_step, time_step)
restart_idx = Z = 2
simInitArgs = dict(
    time_step_nbr=time_step_nbr,
    time_step=time_step,
    cells=cells,
    dl=0.3,
    init_options=dict(dir=f"{first_out}/00000.00{Z}00", mpi_size=first_mpi_size),
    diag_options=dict(format="phareh5", options=dict(dir=secnd_out, mode="overwrite")),
)


def setup_model(sim):
    model = ph.MaxwellianFluidModel(
        protons={"mass": 1, "charge": 1, "nbr_part_per_cell": ppc},
        alpha={"mass": 4, "charge": 1, "nbr_part_per_cell": ppc},
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)
    dump_all_diags(model.populations, timestamps=timestamps)
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
        fidx, sidx = 4, 2
        datahier0 = get_all_available_quantities_from_h5(first_out, timestamps[fidx])
        datahier0.time_hier = {  # swap times
            format_timestamp(timestamps[sidx]): datahier0.time_hier[
                format_timestamp(timestamps[fidx])
            ]
        }
        datahier1 = get_all_available_quantities_from_h5(secnd_out, timestamps[sidx])
        qties = None
        skip = ["protons_patchGhost", "alpha_patchGhost"]
        ds = [single_patch_for_LO(d, qties, skip) for d in [datahier0, datahier1]]
        self.assertTrue(hierarchy_compare(*ds, atol=1e-12))


if __name__ == "__main__":
    if Path("test_init_from_restart").exists():
        unittest.main()
    else:
        print("No source data found - see tools/test_data_gen.py ")
