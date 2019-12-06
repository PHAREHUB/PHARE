def all_iterables(*args):
    """
    return true if all arguments are either lists or tuples
    """
    return all([isinstance(arg, list) or isinstance(arg, tuple) for arg in args])


def none_iterable(*args):
    """
    return true if none of the arguments are either lists or tuples
    """
    return all([not isinstance(arg, list) and not isinstance(arg, tuple) for arg in args])


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
    if none_iterable(arg):
        return [arg]
    else:
        return arg


def not_in_keywords_list(kwd_list,**kwargs):
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
    keys  = [k for k in kwargs.keys()]
    check = [(mk, mk in keys) for mk in mandatory_kwd_list]
    return [mk[0] for mk in check if mk[1] is False]



