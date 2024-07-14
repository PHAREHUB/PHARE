#
#
#


import json

# continue to use override if set
_cpp_lib_override = None


def cpp_lib(override=None):
    import importlib

    global _cpp_lib_override
    if override is not None:
        _cpp_lib_override = override
    if _cpp_lib_override is not None:
        return importlib.import_module(_cpp_lib_override)

    if not __debug__:
        return importlib.import_module("pybindlibs.cpp")
    try:
        return importlib.import_module("pybindlibs.cpp_dbg")
    except ImportError:
        return importlib.import_module("pybindlibs.cpp")


def cpp_etc_lib():
    import importlib

    return importlib.import_module("pybindlibs.cpp_etc")


def build_config():
    return cpp_etc_lib().phare_build_config()


def env_vars():
    return cpp_etc_lib().phare_env_vars()


def print_env_vars_info():
    # see: src/core/env.hpp
    for env_var_name, env_var in cpp_etc_lib().phare_env_vars().items():
        print(f"{env_var_name}: {env_var.desc}")
        print("Options:")
        for option_key, option_val in env_var.options:
            print(f"  {option_key}: {option_val}")
        print("")


def build_config_as_json():
    return json.dumps(build_config())


def splitter_type(dim, interp, n_particles):
    return getattr(cpp_lib(), f"Splitter_{dim}_{interp}_{n_particles}")


def create_splitter(dim, interp, n_particles):
    return splitter_type(dim, interp, n_particles)()


def split_pyarrays_fn(dim, interp, n_particles):
    return getattr(cpp_lib(), f"split_pyarray_particles_{dim}_{interp}_{n_particles}")
