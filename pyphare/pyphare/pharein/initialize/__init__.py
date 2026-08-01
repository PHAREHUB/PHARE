import numpy as np

from pyphare.core.phare_utilities import is_scalar

from . import general, hybrid, mhd

__all__ = ["DictPopulator", "general", "hybrid", "mhd"]


class DictPopulator:
    """
    Thin stateless typed wrapper around pybindlibs.dictator.
    Knows nothing about simulations, models, or what a given path means -
    that logic lives in the populateDict() chunks of general/hybrid/mhd.
    """

    def __init__(self):
        self.pp = dictator()

    def add_int(self, path, val):
        self.pp.add_int(path, int(val))

    def add_bool(self, path, val):
        self.pp.add_bool(path, bool(val))

    def add_double(self, path, val):
        self.pp.add_double(path, float(val))

    def add_size_t(self, path, val):
        casted = int(val)
        if casted < 0:
            raise RuntimeError("DictPopulator::add_size_t received negative value")
        self.pp.add_size_t(path, casted)

    def add_vector_int(self, path, val):
        self.pp.add_vector_int(path, list(val))

    def add_vector_string(self, path, val):
        self.pp.add_vector_string(path, list(val))

    def add_string(self, path, val):
        self.pp.add_string(path, val)

    def add_array_as_vector(self, path, val):
        self.pp.add_array_as_vector(path, val)

    def add_optional_size_t(self, path, val):
        self.pp.add_optional_size_t(path, val)

    def add_init_function(self, ndim, path, fn):
        adder = getattr(self.pp, "addInitFunction{:d}".format(ndim) + "D")
        adder(path, fn_wrapper(fn))


# converts scalars to array of expected size
# converts lists to arrays
class py_fn_wrapper:
    def __init__(self, fn):
        self.fn = fn

    def __call__(self, *xyz):
        args = [np.asarray(arg) for arg in xyz]
        ret = self.fn(*args)
        if isinstance(ret, list):
            ret = np.asarray(ret)
        if is_scalar(ret):
            ret = np.full(len(args[-1]), ret)
        return ret


# Wrap calls to user init functions to turn C++ vectors to ndarrays,
#  and returned ndarrays to C++ span
class fn_wrapper(py_fn_wrapper):
    def __init__(self, fn):
        super().__init__(fn)

    def __call__(self, *xyz):
        from pyphare.cpp import cpp_etc_lib

        # convert numpy array to C++ SubSpan
        # couples vector init functions to C++
        return cpp_etc_lib().makePyArrayWrapper(super().__call__(*xyz))


def dictator():
    # doesn't exist until the project is built!
    import pybindlibs.dictator as pp

    return pp
