#
#
#

import os
import json
import importlib
from . import validate

__all__ = ["validate"]

_libs = {}


def simulator_id(sim):
    parts = [str(sim.ndim)]

    if sim.interp_order:
        parts += [str(sim.interp_order), str(sim.refined_particle_nbr)]

    if sim.mhd_timestepper:
        hall_active = "true" if sim.hall else "false"
        res_active = "true" if sim.res else "false"
        hyper_res_active = "true" if sim.hyper_res else "false"
        parts += [
            sim.mhd_timestepper,
            sim.reconstruction,
            sim.limiter,
            sim.riemann,
            hall_active,
            res_active,
            hyper_res_active,
        ]

    return "_".join(parts)


def cpp_lib(sim):
    global _libs

    mod_str = f"pybindlibs.cpp_{simulator_id(sim)}"
    if mod_str not in _libs:
        _libs[mod_str] = importlib.import_module(mod_str)
    return _libs[mod_str]


def cpp_etc_lib():
    return importlib.import_module("pybindlibs.cpp_etc")


def build_config():
    return cpp_etc_lib().phare_build_config()


def build_config_as_json():
    return json.dumps(build_config())


def splitter_type(sim):
    return getattr(cpp_lib(sim), "Splitter")


def split_pyarrays_fn(sim):
    return getattr(cpp_lib(sim), "split_pyarray_particles")


def mpi_rank():
    return getattr(cpp_etc_lib(), "mpi_rank")()


def mpi_size():
    return getattr(cpp_etc_lib(), "mpi_size")()


def mpi_barrier():
    return getattr(cpp_etc_lib(), "mpi_barrier")()


def mpi_initialized():
    return getattr(cpp_etc_lib(), "mpi_initialized")()


def print_rank0(*args, **kwargs):
    def should_print():
        try:
            if mpi_initialized():
                return mpi_rank() == 0
        except ImportError:
            # missing module or mpi not initialized
            ...
        envs = ["OMPI_COMM_WORLD_RANK", "SLURM_PROCID"]
        for env in envs:
            if env in os.environ:
                return int(os.environ[env]) == 0
        return True  # FALL BACK ALWAYS PRINT

    if should_print():
        print(*args, **kwargs)
