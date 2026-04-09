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
    new_patch_data.dataset = λ(new_patch_data[:])
    return new_patch_data


def drop_patchdata_ghosts(pd, **kwargs):
    ghosts_nbr = [0] * pd.box.ndim
    data = pd[pd.box] if any(pd.ghosts_nbr) else pd[:]
    return pd.copy_as(data, ghosts_nbr=ghosts_nbr)


def drop_patch_ghosts(patch, **kwargs):
    pd_attrs = []
    for name, pd in patch.patch_datas.items():
        pd_attrs.append({"name": name, "data": drop_patchdata_ghosts(pd)})
    return tuple(pd_attrs)


def drop_ghosts(input, **kwargs):
    """accepts patch or specific FieldData"""
    if hasattr(input, "patch_datas"):
        return drop_patch_ghosts(input)
    return drop_patchdata_ghosts(input)


class DataAccessor:
    def __init__(self, hinfo, other):
        self.hinfo = hinfo
        self.other = other

    def get(self, key, key_map):
        hinfo = self.hinfo
        if issubclass(type(self.other), type(hinfo.hier)):
            opatch = self.other.level(hinfo.ilvl, hinfo.time)[hinfo.patch_idx]
            if key_map[key] in opatch.patch_datas:
                return opatch[key_map[key]][:]
            if key in opatch.patch_datas:
                return opatch[key][:]
            for v in key_map.values():
                for pd_key in opatch.patch_datas.keys():
                    if key.endswith(v) and pd_key.endswith(v):
                        return opatch[pd_key][:]
            raise ValueError(opatch.patch_datas.keys())

        return self.other


def _compute_copy_op(patch, op, hinfo, other, reverse=False, key_map=None):
    def _(a, b):
        return op(b, a) if reverse else op(a, b)

    key_map = key_map or {name: name for name in patch.patch_datas.keys()}
    data = DataAccessor(hinfo, other)

    return tuple(
        {
            "name": key_map[name],
            "data": _compute_copy_do(pd, lambda ds: _(ds, data.get(name, key_map))),
        }
        for name, pd in patch.patch_datas.items()
    )


def _compute_copy_rop(patch, op, hinfo, other, key_map=None):
    return _compute_copy_op(patch, op, hinfo, other, reverse=True, key_map=key_map)


def _compute_gaussian_filter_on_scalarfield(patch, sigma, qty, hinfo):
    from scipy.ndimage import gaussian_filter

    ndim = patch[qty].box.ndim
    nb_ghosts = patch[qty].ghosts_nbr[0]
    ds = np.asarray(patch[qty][:])

    ds_ = np.full(list(ds.shape), np.nan)

    gf_ = gaussian_filter(ds, sigma=sigma)
    select = tuple([slice(nb_ghosts, -nb_ghosts) for _ in range(ndim)])
    ds_[select] = np.asarray(gf_[select])

    return ({"name": qty, "data": patch[qty].copy_as(ds_[select])},)


def _compute_gaussian_filter_on_vectorfield(patch, sigma, hinfo):
    from scipy.ndimage import gaussian_filter

    ndim = patch["x"].box.ndim
    nb_ghosts = patch["x"].ghosts_nbr[0]
    ds_x = np.asarray(patch["x"][:])
    ds_y = np.asarray(patch["y"][:])
    ds_z = np.asarray(patch["z"][:])

    dsx_ = np.full(list(ds_x.shape), np.nan)
    dsy_ = np.full(list(ds_y.shape), np.nan)
    dsz_ = np.full(list(ds_z.shape), np.nan)

    gfx_ = gaussian_filter(ds_x, sigma=sigma)
    gfy_ = gaussian_filter(ds_y, sigma=sigma)
    gfz_ = gaussian_filter(ds_z, sigma=sigma)

    select = tuple([slice(nb_ghosts, -nb_ghosts) for _ in range(ndim)])

    dsx_[select] = np.asarray(gfx_[select])
    dsy_[select] = np.asarray(gfy_[select])
    dsz_[select] = np.asarray(gfz_[select])

    return (
        {"name": "x", "data": patch["x"].copy_as(dsx_)},
        {"name": "y", "data": patch["y"].copy_as(dsy_)},
        {"name": "z", "data": patch["z"].copy_as(dsz_)},
    )
