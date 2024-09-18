import copy
import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from tests.simulator import test_restarts
from tests.diagnostic import dump_all_diags
from pathlib import Path

time_step = 0.001
time_step_nbr = 5
final_time = time_step_nbr * time_step
timestamps = np.arange(0, final_time + time_step, time_step)
first_out = str(Path.home() / "phare_data" / "tests/simulator/test_from_init")

if __name__ == "__main__":
    simput = copy.deepcopy(test_restarts.simArgs)
    simput["restart_options"]["dir"] = first_out
    simput["restart_options"]["timestamps"] = timestamps
    simput["diag_options"]["options"]["dir"] = first_out
    sim = ph.Simulation(**simput)
    dump_all_diags(test_restarts.setup_model().populations, timestamps=timestamps)
    Simulator(sim).run()
