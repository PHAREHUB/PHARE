#
#
#

import os
import json
import importlib
from pathlib import Path

_libs = {}
_pybindlibs_checked = False


def _simulator_id_parts(sim):
    """The build permutation fields, in the order the res/sim files list them."""
    parts = [str(sim.ndim)]

    if sim.interp_order:
        parts += [str(sim.interp_order), str(sim.refined_particle_nbr)]

    if sim.mhd_timestepper:
        hall_active = "true" if sim.hall else "false"
        parts += [
            sim.reconstruction,
            sim.limiter,
            sim.riemann,
            hall_active,
        ]

    return parts


def simulator_id(sim):
    return "_".join(_simulator_id_parts(sim))


def permutation_line(sim):
    """The res/sim permutation file line that builds this simulation's module.
    Not derivable from simulator_id by splitting: field values contain '_'
    (SSPRK4_5)."""
    return ",".join(_simulator_id_parts(sim))


def _pythonpath_pybindlibs():
    for entry in os.environ.get("PYTHONPATH", "").split(os.pathsep):
        if entry and (Path(entry) / "pybindlibs" / "__init__.py").is_file():
            return (Path(entry) / "pybindlibs").resolve()
    return None


def _check_pybindlibs_not_shadowed(pybindlibs):
    expected = _pythonpath_pybindlibs()
    if expected is None:
        return
    if pybindlibs.__file__:
        loaded = Path(pybindlibs.__file__).resolve().parent
    else:
        loaded = Path(next(iter(pybindlibs.__path__))).resolve()
    if loaded != expected:
        raise ImportError(
            f"pybindlibs was imported from {loaded}, not from {expected} found on PYTHONPATH.\n"
            "An editable install of PHARE in this environment (uv sync, uv run or pip install -e)"
            " takes precedence over PYTHONPATH.\n"
            "Remove it from this environment (uv sync --no-install-project, or"
            " pip uninstall pyphare), or unset PYTHONPATH to use the installed build."
        )


def _pybindlibs():
    global _pybindlibs_checked

    import pybindlibs

    if not _pybindlibs_checked:
        _check_pybindlibs_not_shadowed(pybindlibs)
        _pybindlibs_checked = True
    return pybindlibs


def built_simulator_ids():
    """Permutations present in the loaded build, or None if pybindlibs itself is
    not importable (no build on PYTHONPATH at all)."""
    import pkgutil

    try:
        pybindlibs = _pybindlibs()
    except ModuleNotFoundError:
        return None

    return sorted(
        name[len("cpp_") :]
        for _, name, _ in pkgutil.iter_modules(pybindlibs.__path__)
        if name.startswith("cpp_") and name != "cpp_etc"
    )


def _pybindlibs_location():
    try:
        pybindlibs = _pybindlibs()
    except ModuleNotFoundError:
        return "pybindlibs not importable"

    return ", ".join(pybindlibs.__path__)


def _no_such_module_message(sim, mod_str):
    built = built_simulator_ids()
    if built is None:
        listing = "    pybindlibs not importable -- no build on PYTHONPATH"
    elif built:
        listing = "\n".join(f"    {name}" for name in built)
    else:
        listing = "    (none)"
    return (
        f"No module named '{mod_str}'.\n"
        f"This run needs the build permutation\n"
        f"    {permutation_line(sim)}\n"
        f"pybindlibs loaded from: {_pybindlibs_location()}\n"
        f"To add it, either:\n"
        f"  - CMake build: add the line to a res/sim permutation file (res/sim/all.txt) and rebuild\n"
        f"  - installed package: put the line in a file and reinstall with\n"
        f"    pip install . -C cmake.define.PHARE_PERMUTATIONS=<file>\n"
        f"Permutations currently built:\n{listing}"
    )


def cpp_lib(sim):
    global _libs

    mod_str = f"pybindlibs.cpp_{simulator_id(sim)}"
    if mod_str not in _libs:
        try:
            _pybindlibs()
            _libs[mod_str] = importlib.import_module(mod_str)
        except ModuleNotFoundError as e:
            raise ModuleNotFoundError(_no_such_module_message(sim, mod_str)) from e
    return _libs[mod_str]


def cpp_etc_lib():
    _pybindlibs()
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
