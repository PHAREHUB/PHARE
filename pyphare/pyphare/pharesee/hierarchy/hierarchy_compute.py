#
#
#

import operator
from copy import deepcopy

from .hierarchy import PatchHierarchy


def rename(hierarchy, names):
    from .hierarchy_utils import compute_hier_from

    return compute_hier_from(compute_rename, hierarchy, new_names=names)


def compute_rename(patch, **kwargs):
    new_names = kwargs["new_names"]
    pd_attrs = []

    for new_name, pd_name in zip(new_names, patch.patch_datas):
        pd_attrs.append(patch[pd_name].copy_as(name=new_name))

    return tuple(pd_attrs)


def compute_mul(patch, **kwargs):
    return _compute_copy_op(patch, operator.__mul__, **kwargs)


def compute_add(patch, **kwargs):
    return _compute_copy_op(patch, operator.__add__, **kwargs)


def compute_sub(patch, **kwargs):
    return _compute_copy_op(patch, operator.__sub__, **kwargs)


def compute_truediv(patch, **kwargs):
    return _compute_copy_op(patch, operator.__truediv__, **kwargs)


def compute_rtruediv(patch, **kwargs):
    return _compute_copy_rop(patch, operator.__truediv__, **kwargs)


def _compute_copy_do(patch_data, λ):
    return patch_data.copy_as(λ(patch_data.dataset[:]))


def drop_ghosts(patch, **kwargs):
    pd_attrs = []
    ghosts_nbr = [0] * patch.box.ndim
    for name, pd in patch.patch_datas.items():
        data = pd[patch.box] if any(pd.ghosts_nbr) else pd[:]
        pd_attrs.append(pd.copy_as(data, ghosts_nbr=ghosts_nbr))
    return tuple(pd_attrs)


class DataAccessor:
    """
    Resolves the right-hand operand of a patch-wise binary operation.

    `accessor` is a HierarchyAccessor locating the patch currently being
    computed (hierarchy, time, level, patch index). `operand` is either
    another hierarchy, in which case indexing by quantity name returns its
    dataset at that same patch location, or anything else operable against
    a dataset (usually a scalar), which is returned unchanged regardless of
    the key.
    """

    def __init__(self, accessor, operand):
        self.accessor = accessor
        self.operand = operand

    def __getitem__(self, key):
        if isinstance(self.operand, PatchHierarchy):
            return self.operand[self.accessor][key].dataset[:]
        return self.operand


def _compute_copy_op(patch, op, accessor, operand, reverse=False):
    def _(a, b):
        return op(b, a) if reverse else op(a, b)

    data = DataAccessor(accessor, operand)
    return tuple(
        _compute_copy_do(pd, lambda ds: _(ds, data[name]))
        for name, pd in patch.patch_datas.items()
    )


def _compute_copy_rop(patch, op, accessor, operand):
    return _compute_copy_op(patch, op, accessor, operand, reverse=True)
