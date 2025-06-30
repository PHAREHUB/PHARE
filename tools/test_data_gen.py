import sys
import copy
import unittest
import subprocess
import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from tests.simulator import test_restarts
from tests.diagnostic import dump_all_diags

output_dir = "phare_data_gen/"


def run_restartable_data_sim():
    from tests.simulator import test_init_from_restart as restartable_sim

    """uses params from tests_restarts.py"""
    out = output_dir + "/tests/simulator/test_init_from_restart"
    simput = copy.deepcopy(test_restarts.simArgs)
    simput["restart_options"]["dir"] = out
    simput["restart_options"]["timestamps"] = restartable_sim.timestamps
    simput["diag_options"]["options"]["dir"] = out
    sim = ph.Simulation(**simput)
    dump_all_diags(
        test_restarts.setup_model().populations, timestamps=restartable_sim.timestamps
    )
    Simulator(sim).run()


def launch(fn, n=5):
    """Launch secondary process to run first simulation to avoid initalizing MPI now"""
    cmd = f"mpirun -n {n} python3 -O {__file__} {fn}"
    print(cmd)
    try:
        p = subprocess.run(cmd.split(" "), check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print("CalledProcessError", e)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        globals()[sys.argv[1]]()
    else:
        for k in [k for k in list(globals().keys()) if k.startswith("run_")]:
            launch(k)
#
# tar czf phare_data_gen.tar.xz phare_data_gen
# mc put phare_data_gen.tar.xz minio/phare/phare_data_gen.tar.xz --insecure
