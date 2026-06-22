#
#
#

import operator
import numpy as np
from copy import deepcopy


def rename(hierarchy, names):
    from .hierarchy_utils import compute_hier_from

    return compute_hier_from(compute_rename, hierarchy, new_names=names)


def compute_rename(patch, **kwargs):
    new_names = kwargs["new_names"]
    pd_attrs = []

    for new_name, pd_name in zip(new_names, patch.patch_datas):
        pd_attrs.append({"name": new_name, "data": patch[pd_name]})

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
    new_patch_data = deepcopy(patch_data)
    new_patch_data.dataset = λ(patch_data.dataset[:])
    return new_patch_data


def drop_ghosts(patch, **kwargs):
    pd_attrs = []
    ghosts_nbr = [0] * patch.box.ndim
    for name, pd in patch.patch_datas.items():
        data = pd[patch.box] if any(pd.ghosts_nbr) else pd[:]
        pd_attrs.append(
            {
                "name": name,
                "data": pd.copy_as(data, ghosts_nbr=ghosts_nbr),
            }
        )
    return tuple(pd_attrs)


class DataAccessor:
    def __init__(self, hinfo, other):
        self.hinfo = hinfo
        self.other = other

    def __getitem__(self, key):
        hinfo = self.hinfo
        hier_base = next(
            c for c in reversed(type(hinfo.hier).__mro__) if c is not object and hasattr(c, "level")
        )
        if isinstance(self.other, hier_base):
            return self.other.level(hinfo.ilvl, hinfo.time)[hinfo.patch_idx][
                key
            ].dataset[:]
        return self.other


def _compute_copy_op(patch, op, hinfo, other, reverse=False):
    def _(a, b):
        return op(b, a) if reverse else op(a, b)

    data = DataAccessor(hinfo, other)
    return tuple(
        {
            "name": name,
            "data": _compute_copy_do(pd, lambda ds: _(ds, data[name])),
        }
        for name, pd in patch.patch_datas.items()
    )


def _compute_copy_rop(patch, op, hinfo, other):
    return _compute_copy_op(patch, op, hinfo, other, reverse=True)


def _compute_gaussian_filter_on_scalarfield(patch, **kwargs):
    from scipy.ndimage import gaussian_filter

    ref_pd = patch["value"]
    ndim = ref_pd.box.ndim
    nb_ghosts = ref_pd.ghosts_nbr[0]
    sigma = kwargs["sigma"]
    ds = np.asarray(ref_pd[:])

    ds_ = np.full(list(ds.shape), np.nan)
    gf_ = gaussian_filter(ds, sigma=sigma)
    select = tuple([slice(nb_ghosts, -nb_ghosts) for _ in range(ndim)])
    ds_[select] = np.asarray(gf_[select])

    return ({"name": "value", "data": ref_pd.copy_as(ds_)},)


def _compute_gaussian_filter_on_vectorfield(patch, **kwargs):
    from scipy.ndimage import gaussian_filter

    ref_name = next(iter(patch.patch_datas.keys()))
    ref_pd = patch[ref_name]
    ndim = ref_pd.box.ndim
    nb_ghosts = ref_pd.ghosts_nbr[0]
    sigma = kwargs["sigma"]
    select = tuple([slice(nb_ghosts, -nb_ghosts) for _ in range(ndim)])

    results = []
    for comp in ("x", "y", "z"):
        ds = np.asarray(patch[comp][:])
        ds_ = np.full(list(ds.shape), np.nan)
        ds_[select] = np.asarray(gaussian_filter(ds, sigma=sigma)[select])
        results.append({"name": comp, "data": patch[comp].copy_as(ds_)})

    return tuple(results)
