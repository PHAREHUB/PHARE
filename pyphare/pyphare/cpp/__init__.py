"""
Python interface to PHARE C++ bindings.

This module manages dynamic loading of simulator-specific compiled C++ modules.
Each simulator configuration (dimension, interpolation order, refined particle number)
has its own dedicated compiled module to reduce binary size and compilation time.
"""

import json
import importlib
from . import validate

__all__ = ["validate"]

# Cache for loaded C++ simulator modules to avoid re-importing
_libs = {}


def simulator_id(sim):
    """
    Generate unique identifier for a simulator configuration.
    
    Returns a string matching compiled module names (e.g., '1_1_2' for a 1D simulator
    with interpolation order 1 and 2 refined particles).
    
    Args:
        sim: Simulator object with ndim, interp_order, and refined_particle_nbr attributes
        
    Returns:
        String identifier in format "{ndim}_{interp_order}_{refined_particle_nbr}"
    """
    return f"{sim.ndim}_{sim.interp_order}_{sim.refined_particle_nbr}"


def cpp_lib(sim):
    """
    Load and cache the C++ simulator module for a given simulator configuration.
    
    Each simulator template permutation (dimension × interpolation order × refined particles)
    is compiled into a separate pybind11 module. This function lazy-loads the appropriate
    module and caches it to avoid repeated imports.
    
    Args:
        sim: Simulator object with ndim, interp_order, and refined_particle_nbr attributes
        
    Returns:
        The imported C++ module corresponding to the simulator configuration
    """
    global _libs

    mod_str = f"pybindlibs.cpp_{simulator_id(sim)}"
    if mod_str not in _libs:
        _libs[mod_str] = importlib.import_module(mod_str)
    return _libs[mod_str]


def cpp_etc_lib():
    """
    Load the C++ utilities module (non-simulator-specific bindings).
    
    This module contains shared functionality like MPI utilities, hierarchy management,
    and version information that doesn't depend on simulator template parameters.
    
    Returns:
        The cpp_etc pybind11 module
    """
    return importlib.import_module("pybindlibs.cpp_etc")


def build_config():
    return cpp_etc_lib().phare_build_config()


def build_config_as_json():
    return json.dumps(build_config())


def splitter_type(sim):
    """Get the Splitter class for the given simulator configuration."""
    return getattr(cpp_lib(sim), "Splitter")


def split_pyarrays_fn(sim):
    """Get the split_pyarray_particles function for the given simulator configuration."""
    return getattr(cpp_lib(sim), "split_pyarray_particles")


def mpi_rank():
    return getattr(cpp_etc_lib(), "mpi_rank")()


def mpi_size():
    return getattr(cpp_etc_lib(), "mpi_size")()


def mpi_barrier():
    return getattr(cpp_etc_lib(), "mpi_barrier")()
