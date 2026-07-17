#
#
#

import json
import importlib
from . import validate

__all__ = ["validate"]

_libs = {}


def simulator_id(sim, layout=None, allocator=None):
    layout = str(layout) if layout else "AoSMapped"
    allocator = str(allocator) if allocator else "CPU"

    if getattr(sim, "mhd_timestepper", None):
        Hall = "true" if sim.hall else "false"
        Res = "true" if sim.res else "false"
        Hyper_Res = "true" if sim.hyper_res else "false"
        return (
            f"{sim.ndim}_{sim.interp_order}_{sim.refined_particle_nbr}_{layout}_{allocator}_"
            f"{sim.mhd_timestepper}_{sim.reconstruction}_{sim.limiter}_"
            f"{sim.riemann}_{Hall}_{Res}_{Hyper_Res}"
        )

    if type(layout) is str and "." in layout:
        layout = str(layout).split(".")[1] if layout else None
    return "_".join(
        str(s)
        for s in [
            sim.ndim,
            sim.interp_order,
            sim.refined_particle_nbr,
            layout,
            allocator,
        ]
    )


def cpp_lib(sim, layout=None, allocator=0):
    global _libs

    mod_str = f"pybindlibs.cpp_{simulator_id(sim, layout, allocator)}"
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


def mpi_initialized():
    return getattr(cpp_etc_lib(), "mpi_initialized")()


def mpi_rank():
    return getattr(cpp_etc_lib(), "mpi_rank")()


def mpi_size():
    return getattr(cpp_etc_lib(), "mpi_size")()


def mpi_barrier():
    return getattr(cpp_etc_lib(), "mpi_barrier")()


def layout_mode_enum():
    return getattr(cpp_etc_lib(), "LayoutMode")()


def supported_particle_layouts():
    return getattr(cpp_etc_lib(), "supported_layouts")()
