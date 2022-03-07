
import os
import numpy as np
from pyphare.core import phare_utilities


# ------------------------------------------------------------------------------


def validate(sim):

    sim.restart_options["timestamps"] = phare_utilities.np_array_ify(sim.restart_options["timestamps"])
    init_time = sim.start_time()

    timestamps = sim.restart_options["timestamps"]
    if np.any(timestamps < init_time):
        raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be less than simulation.init_time({init_time}))")
    if np.any(timestamps > sim.final_time):
        raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be greater than simulation.final_time({sim.final_time}))")
    if not np.all(np.diff(timestamps) >= 0):
        raise RuntimeError(f"Error: {clazz}.{key} not in ascending order)")
    if not np.all(np.abs(timestamps / sim.time_step - np.rint(timestamps/sim.time_step) < 1e-9)):
        raise RuntimeError(f"Error: {clazz}.{key} is inconsistent with simulation.time_step")

    sim.restart_options["timestamps"] = conserve_existing(sim)

# ------------------------------------------------------------------------------



def conserve_existing(sim):
    """
      trim timestamps from array if files exist for that time and mode is "conserve"
    """

    restart_options = sim.restart_options
    timestamps = restart_options["timestamps"]

    if "mode" in restart_options:

        if restart_options["mode"] == "conserve":
            from pyphare.cpp import cpp_etc_lib

            torm = []

            for i, time in enumerate(timestamps):

                restart_file = cpp_etc_lib().restart_path_for_time(sim.restart_file_path(), time)

                if os.path.exists(restart_file):
                    torm += [i]

            for i in reversed(torm):
                timestamps = np.delete(timestamps, i)

    return timestamps


# ------------------------------------------------------------------------------

