
import os
import numpy as np

from ..core import phare_utilities
from . import global_vars

# ------------------------------------------------------------------------------


def restarts_checker(func):
    def wrapper(restarts_object, **kwargs):

        mandatory_keywords = ['write_timestamps']

        # check if some mandatory keywords are not missing
        missing_mandatory_kwds = phare_utilities.check_mandatory_keywords(mandatory_keywords, **kwargs)
        if len(missing_mandatory_kwds) > 0:
            raise RuntimeError("Error: missing mandatory parameters : " + ', '.join(missing_mandatory_kwds))

        accepted_keywords = mandatory_keywords

        # check that all passed keywords are in the accepted keyword list
        wrong_kwds = phare_utilities.not_in_keywords_list(accepted_keywords, **kwargs)
        if len(wrong_kwds) > 0:
            raise RuntimeError("Error: invalid arguments - " + " ".join(wrong_kwds))

        return func(restarts_object, **kwargs)

    return wrapper


# ------------------------------------------------------------------------------


def validate_timestamps(clazz, **kwargs):
    sim = global_vars.sim

    for key in ["write_timestamps"]:
        timestamps = kwargs[key]

        if np.any(timestamps < sim.init_time):
            raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be less than simulation.init_time({sim.init_time}))")
        if np.any(timestamps > sim.final_time):
            raise RuntimeError(f"Error: timestamp({sim.time_step_nbr}) cannot be greater than simulation.final_time({sim.final_time}))")
        if not np.all(np.diff(timestamps) >= 0):
            raise RuntimeError(f"Error: {clazz}.{key} not in ascending order)")
        if not np.all(np.abs(timestamps / sim.time_step - np.rint(timestamps/sim.time_step) < 1e-9)):
            raise RuntimeError(f"Error: {clazz}.{key} is inconsistent with simulation.time_step")


# ------------------------------------------------------------------------------


class Restarts(object):

    @restarts_checker
    def __init__(self, **kwargs):

        if global_vars.sim is None:
            raise RuntimeError("A simulation must be created before adding restarts")

        validate_timestamps(self.__class__.__name__, **kwargs)
        self.write_timestamps = self._conserve_existing(kwargs['write_timestamps'])
        global_vars.sim.add_restarts(self)



    def _conserve_existing(self, timestamps):
        """
          trim timestamps from array if files exist for that time and mode is "conserve"
        """

        sim = global_vars.sim
        restart_options = sim.restart_options

        if "mode" in restart_options:

            if restart_options["mode"] == "conserve":
                from pyphare.cpp import cpp_etc_lib

                torm = []

                for i, time in enumerate(timestamps):
                    restart_file = cpp_etc_lib().restart_path_for_time(sim.restart_file_path(), time)
                    if os.path.exists(restart_file):
                        torm += [i]

                for i in reversed(torm):
                    del timestamps[i]

        return timestamps


# ------------------------------------------------------------------------------

