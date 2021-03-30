

def cpp_lib():
    import importlib
    if not __debug__:
        return importlib.import_module("pybindlibs.cpp")
    try:
        return importlib.import_module("pybindlibs.cpp_dbg")
    except ImportError as err:
        return importlib.import_module("pybindlibs.cpp")


def splitter_type(dim, interp, n_particles):
    return getattr(cpp_lib(), f"Splitter_{dim}_{interp}_{n_particles}")

def create_splitter(dim, interp, n_particles):
    return splitter_type(dim, interp, n_particles)()


def split_pyarrays_fn(dim, interp, n_particles):
    return getattr(cpp_lib(), f"split_pyarray_particles_{dim}_{interp}_{n_particles}")
