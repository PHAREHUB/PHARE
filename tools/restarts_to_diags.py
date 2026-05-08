#
#  Convert samrai restarts to diagnostic files for that timestamp
#
#  Must be run with the same number of cores which *created* the restart files.
#
#  use like:
#   mpirun -n "$CORES" python3 tools/restarts_to_diags.py [$input $format $output]
#

import sys
import numpy as np
from pathlib import Path
from copy import deepcopy


file_path = sys.argv[1] if len(sys.argv) > 1 else "checkpoints"
outformat = sys.argv[2] if len(sys.argv) > 2 else "phareh5"
outpath = sys.argv[3] if len(sys.argv) > 3 else f"phare_outputs/{outformat}"

diag_options = dict(format=outformat, options=dict(dir=outpath))
restart_options = {
    "dir": file_path,
    "mode": "conserve",  # shouldn't matter as we do no advances
    "restart_time": 0,  # overwritten as we go
}


def do(restart_file_load_path: Path):
    import pyphare.pharein as ph
    from pyphare.cpp import cpp_etc_lib
    from pyphare.simulator.simulator import Simulator
    from pyphare.pharein.simulation import deserialize as deserialize_sim

    Path(outpath).mkdir(exist_ok=True, parents=True)
    time = float(restart_file_load_path.name)

    sim = ph.global_vars.sim = deserialize_sim(
        cpp_etc_lib().serialized_simulation_string(str(restart_file_load_path))
    )

    sim.restart_options = deepcopy(restart_options)
    sim.restart_options["restart_time"] = time
    sim.diag_options = deepcopy(diag_options)

    for k, v in sim.diagnostics.items():
        v.write_timestamps = np.asarray([time])

    Simulator(sim).setup().initialize().reset()


def main():
    for item in Path(file_path).iterdir():
        do(item)


if __name__ == "__main__":
    try:
        from pyphare.simulator.simulator import startMPI

        startMPI()
        main()
    except Exception as e:
        import traceback

        print(f"Exception caught in phare/tools/restarts_to_diags.py: \n{e}")
        print(traceback.format_exc())
        sys.exit(1)
