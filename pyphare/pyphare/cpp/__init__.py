#
#
#


import json
import importlib
from pyphare.cpp import validate

__all__ = ["validate"]

_libs = {}


def cpp_lib(sim_str, override=None):
    global _libs
    if sim_str in _libs:
        return _libs[sim_str]

    _libs[sim_str] = importlib.import_module(f"pybindlibs.cpp_{sim_str}")
    return _libs[sim_str]


def simulator_id(dim, interp, nbrRefinedPart):
    return f"{dim}_{interp}_{nbrRefinedPart}"


def cpp_etc_lib():
    import importlib

    return importlib.import_module("pybindlibs.cpp_etc")


def build_config():
    return cpp_etc_lib().phare_build_config()


def build_config_as_json():
    return json.dumps(build_config())


def splitter_type(dim, interp, nbrRefinedPart):
    return getattr(cpp_lib(simulator_id(dim, interp, nbrRefinedPart)), "Splitter")


def create_splitter(dim, interp, nbrRefinedPart):
    return splitter_type(dim, interp, nbrRefinedPart)()


def split_pyarrays_fn(dim, interp, nbrRefinedPart):
    return getattr(
        cpp_lib(simulator_id(dim, interp, nbrRefinedPart)), "split_pyarray_particles"
    )


def mpi_rank():
    return getattr(cpp_etc_lib(), "mpi_rank")()


def mpi_size():
    return getattr(cpp_etc_lib(), "mpi_size")()


def mpi_barrier():
    return getattr(cpp_etc_lib(), "mpi_barrier")()
