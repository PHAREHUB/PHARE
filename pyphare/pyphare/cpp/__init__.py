#
#
#

import json
import importlib
from . import validate

__all__ = ["validate"]

_libs = {}


def simulator_id(sim):
    """Generate a unique identifier string for a simulator configuration.
    
    Args:
        sim: Simulator object with ndim, interp_order, and refined_particle_nbr attributes
        
    Returns:
        str: Identifier in format "{ndim}_{interp_order}_{refined_particle_nbr}"
             e.g., "1_2_3" for 1D, order 2 interpolation, 3 refined particles
    """
    return f"{sim.ndim}_{sim.interp_order}_{sim.refined_particle_nbr}"


def cpp_lib(sim):
    """Load and cache the C++ library module for a specific simulator configuration.
    
    Each simulator permutation (ndim, interp_order, refined_particle_nbr) has its own
    compiled C++ module. This function loads the appropriate module and caches it
    to avoid repeated imports.
    
    Args:
        sim: Simulator object defining the configuration
        
    Returns:
        module: The imported pybind11 module for this simulator configuration
    """
    global _libs

    mod_str = f"pybindlibs.cpp_{simulator_id(sim)}"
    # Cache modules to avoid repeated imports of the same configuration
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


# MPI utility functions - thin wrappers around cpp_etc_lib MPI bindings
# These are independent of simulator configuration and use the cpp_etc module


def mpi_rank():
    """Get the MPI rank of the current process."""
    return getattr(cpp_etc_lib(), "mpi_rank")()


def mpi_size():
    """Get the total number of MPI processes."""
    return getattr(cpp_etc_lib(), "mpi_size")()


def mpi_barrier():
    """Synchronize all MPI processes at a barrier."""
    return getattr(cpp_etc_lib(), "mpi_barrier")()
