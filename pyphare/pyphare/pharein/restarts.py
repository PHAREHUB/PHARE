import os
import numpy as np
from pathlib import Path

from pyphare.core import phare_utilities


def format_time(time):
    """
    0.006 == "00000.00600"
    """
    return "{:0>11.5f}".format(time)


def dump(simulator, time, time_step):
    new_restart_made = simulator.cpp_sim.dump_restarts(time, time_step)
    if new_restart_made:
        try_delete_obsolete_restarts(simulator)


def try_delete_obsolete_restarts(simulator):
    sim = simulator.simulation
    restart_options = sim.restart_options

    if "keep_last" not in restart_options:
        return  # nothing to do

    keep_last = restart_options["keep_last"]
    directory = restart_options.get("dir", ".")

    dirs = []
    for path_object in Path(directory).iterdir():
        if path_object.is_dir():
            try:
                dirs.append(float(path_object.name))
            except ValueError:
                ...  # skip

    if len(dirs) < keep_last:
        return  # nothing to do

    dirs = sorted(dirs)
    to_rm = len(dirs) - keep_last
    assert to_rm >= 0

    import shutil

    for i in range(to_rm):
        try:
            shutil.rmtree(str(Path(directory) / format_time(dirs[i])))
        except OSError as e:
            import warnings

            warnings.warn(f"Failed to remove restart directory {dirs[i]}: {e}")


def restart_time(restart_options):
    if "restart_time" in restart_options:
        if restart_options["restart_time"] == "auto":
            return find_latest_time_from_restarts(restart_options)
        return restart_options["restart_time"]
    return None


def find_latest_time_from_restarts(restart_options):
    directory = restart_options.get("dir", ".")

    dirs = []
    for path_object in Path(directory).iterdir():
        if path_object.is_dir():
            try:
                dirs.append(float(path_object.name))
            except ValueError:
                ...  # skipped

    return None if len(dirs) == 0 else sorted(dirs)[-1]


# ------------------------------------------------------------------------------


def validate(sim):
    restart_options = sim.restart_options

    if "elapsed_timestamps" in restart_options:
        import datetime

        restart_options["elapsed_timestamps"] = [
            int(ts.total_seconds()) if isinstance(ts, datetime.timedelta) else ts
            for ts in phare_utilities.np_array_ify(
                restart_options["elapsed_timestamps"]
            )
        ]

        if not np.all(np.diff(restart_options["elapsed_timestamps"]) >= 0):
            raise RuntimeError(
                "Error: restart_options elapsed_timestamps not in ascending order)"
            )

    if "timestamps" in restart_options:
        restart_options["timestamps"] = phare_utilities.np_array_ify(
            restart_options["timestamps"]
        )
        init_time = sim.start_time()

        timestamps = restart_options["timestamps"]
        if np.any(timestamps < init_time):
            raise RuntimeError(
                f"Error: timestamp({sim.time_step_nbr}) cannot be less than simulation.init_time({init_time}))"
            )
        if np.any(timestamps > sim.final_time):
            raise RuntimeError(
                f"Error: timestamp({sim.time_step_nbr}) cannot be greater than simulation.final_time({sim.final_time}))"
            )
        if not np.all(np.diff(timestamps) >= 0):
            raise RuntimeError(
                "Error: restart_options timestamps not in ascending order)"
            )
        if not np.all(
            np.abs(
                timestamps / sim.time_step - np.rint(timestamps / sim.time_step) < 1e-9
            )
        ):
            raise RuntimeError(
                "Error: restart_options timestamps is inconsistent with simulation.time_step"
            )

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
                restart_file = cpp_etc_lib().restart_path_for_time(
                    sim.restart_file_path(), time
                )

                if os.path.exists(restart_file):
                    torm += [i]

            for i in reversed(torm):
                timestamps = np.delete(timestamps, i)

    return timestamps


# ------------------------------------------------------------------------------


def is_restartable_compared_to(curr_sim, prev_sim):
    import operator

    failed = []

    def _do(op, keys):
        for key in keys:
            curr = getattr(curr_sim, key)
            prev = getattr(prev_sim, key)
            if any([op(curr, prev)]):
                failed.append((key, curr, prev, op))

    # use negative operator for printing
    _do(operator.ne, ["cells", "dl", "max_nbr_levels"])

    if failed:
        print("ERROR: Simulation not restartable - variable mismatch")
        for key, curr, prev, op in failed:
            print(f"{key} current({curr}) {op.__name__} previous({prev})")

    return not failed
