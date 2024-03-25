import math
import numpy as np


def debug_print(*args):
    if __debug__:
        print(*args)


def all_iterables(*args):
    """
    return true if all arguments are either lists or tuples
    """
    return all([isinstance(arg, list) or isinstance(arg, tuple) for arg in args])


def none_iterable(*args):
    """
    return true if none of the arguments are either lists or tuples
    """
    return all(
        [not isinstance(arg, list) and not isinstance(arg, tuple) for arg in args]
    )


def equal_size(*args):
    """
    return true if all arguments are of equal length
    """
    sizes = [len(arg) for arg in args]
    first = sizes[0]
    return all(x == first for x in sizes)


def check_iterables(*args):
    if not (all_iterables(*args) or none_iterable(*args)):
        raise ValueError


def check_equal_size(*args):
    if all_iterables(*args):
        if not equal_size(*args):
            raise ValueError


def listify(arg):
    if isinstance(arg, np.ndarray):
        return arg
    if none_iterable(arg):
        return [arg]
    return arg


def is_scalar(arg):
    return not isinstance(arg, (list, tuple)) and not is_nd_array(arg)


def is_nd_array(arg):
    return isinstance(arg, np.ndarray)


def np_array_ify(arg, size=1):
    if is_scalar(arg):
        return np.asarray([arg] * size)
    if not is_nd_array(arg):
        return np.asarray(arg)
    return arg


refinement_ratio = 2


def not_in_keywords_list(kwd_list, **kwargs):
    """
    return the list of kwargs keys that are not in 'kwd_list'
    """
    keys = [k for k in kwargs.keys()]
    isIn = [k in kwd_list for k in keys]
    wrong_kwds = [keys[i] for i, wrong in enumerate(isIn) if not wrong]
    return wrong_kwds


def check_mandatory_keywords(mandatory_kwd_list, **kwargs):
    """
    return those of mandatory_kwd_list not found in the kwargs keys
    """
    keys = [k for k in kwargs.keys()]
    check = [(mk, mk in keys) for mk in mandatory_kwd_list]
    return [mk[0] for mk in check if mk[1] is False]


def fp_equal(a, b, atol=1e-6):
    return math.isclose(a, b, abs_tol=atol)


def fp_less_equal(a, b, atol=1e-6):
    return fp_equal(a, b, atol=atol) or a < b


def fp_gtr_equal(a, b, atol=1e-6):
    return fp_equal(a, b, atol=atol) or a > b


class FloatingPoint_comparator:
    def __init__(self, fp, atol=1e-6):
        self.fp = fp
        self.atol = atol

    def __eq__(self, other):
        return fp_equal(self.fp, other.fp, self.atol)

    def __lt__(self, other):
        return self.fp < other.fp

    def __le__(self, other):
        return fp_less_equal(self.fp, other.fp, self.atol)

    def __gt__(self, other):
        return self.fp > other.fp

    def __ge__(self, other):
        return fp_gtr_equal(self.fp, other.fp, self.atol)


def is_fp32(item):
    if is_nd_array(item):
        return item.dtype == np.single
    return isinstance(item, float)


def assert_fp_any_all_close(a, b, atol=1e-16, rtol=0, atol_fp32=None):
    if any([is_fp32(el) for el in [a, b]]):
        atol = atol_fp32 if atol_fp32 else atol * 1e8
    np.testing.assert_allclose(a, b, atol=atol, rtol=rtol)


def decode_bytes(input, errors="ignore"):
    return input.decode("ascii", errors=errors)


def run_cli_cmd(cmd, shell=True, capture_output=True, check=False, print_cmd=False):
    """
    https://docs.python.org/3/library/subprocess.html
    """
    import subprocess

    if print_cmd:
        print(f"running: {cmd}")
    try:
        return subprocess.run(
            cmd, shell=shell, capture_output=capture_output, check=check
        )
    except subprocess.CalledProcessError as e:  # only triggers on failure if check=True
        raise RuntimeError(decode_bytes(e.stderr))


def git_hashes(N=1):
    return decode_bytes(
        run_cli_cmd(f"git log -{N} --pretty=format:%h").stdout
    ).splitlines()


def top_git_hash():
    hashes = git_hashes(1)
    if len(hashes) > 0:
        return hashes[0]
    return "master"  # github actions fails?


def print_trace():
    import sys, traceback

    _, _, tb = sys.exc_info()
    traceback.print_tb(tb)


def deep_copy(item, memo, excludes=[]):
    from copy import deepcopy

    clazz = item.__class__
    that = clazz.__new__(clazz)
    memo[id(item)] = that
    for key, value in item.__dict__.items():
        # some objects are not picklable so cannot be deepcopied eg. h5py datasets
        if key in excludes:
            setattr(that, key, value)
        else:
            setattr(that, key, deepcopy(value, memo))
    return that
