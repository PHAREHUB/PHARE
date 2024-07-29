import os
import numpy as np
from pyphare.core import phare_utilities


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

        # seconds_in_an_hour = 60 ** 2
        # for cmp_idx, ref_ts in enumerate(restart_options["elapsed_timestamps"][1:]):
        #     cmp_ts = restart_options["elapsed_timestamps"][cmp_idx]
        #     if ref_ts - cmp_ts < seconds_in_an_hour:
        #         raise RuntimeError("Error: time betweeen restart_options elapsed_timestamps must be at least one hour)")

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
