#
#
#


def cpp_lib():
    import importlib

    return importlib.import_module("pybindlibs.cpp")


def cpp_etc_lib():
    import importlib

    return importlib.import_module("pybindlibs.cpp_etc")


def build_config():
    return cpp_etc_lib().phare_build_config()


def build_config_as_json():
    import json

    return json.dumps(build_config())


def splitter_type(dim, interp, n_particles):
    return getattr(cpp_lib(), f"Splitter_{dim}_{interp}_{n_particles}")


def create_splitter(dim, interp, n_particles):
    return splitter_type(dim, interp, n_particles)()


def split_pyarrays_fn(dim, interp, n_particles):
    return getattr(cpp_lib(), f"split_pyarray_particles_{dim}_{interp}_{n_particles}")
