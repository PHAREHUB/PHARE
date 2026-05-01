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
    select = tuple([slice(nb_ghosts or None, -nb_ghosts or None) for _ in range(ndim)])
    ds_[select] = np.asarray(gf_[select])

    return ({"name": qty, "data": patch[qty].copy_as(ds_)},)


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

    select = tuple([slice(nb_ghosts or None, -nb_ghosts or None) for _ in range(ndim)])

    dsx_[select] = np.asarray(gfx_[select])
    dsy_[select] = np.asarray(gfy_[select])
    dsz_[select] = np.asarray(gfz_[select])

    return (
        {"name": "x", "data": patch["x"].copy_as(dsx_)},
        {"name": "y", "data": patch["y"].copy_as(dsy_)},
        {"name": "z", "data": patch["z"].copy_as(dsz_)},
    )


# ---- operators (formerly pyphare.core.operators) ----


def _compute_dot_product(patch0, hinfo, other, **kwargs):
    patch1 = other.level(hinfo.ilvl, hinfo.time)[hinfo.patch_idx]
    ref_name = next(iter(patch0.patch_datas.keys()))
    dset = (
        patch0["x"][:] * patch1["x"][:]
        + patch0["y"][:] * patch1["y"][:]
        + patch0["z"][:] * patch1["z"][:]
    )

    return ({"name": "value", "data": patch0[ref_name].copy_as(dset)},)


def _compute_sqrt(patch, **kwargs):
    ref_name = next(iter(patch.patch_datas.keys()))

    dset = np.sqrt(patch[ref_name][:])

    return ({"name": "value", "data": patch[ref_name].copy_as(dset)},)


def _compute_cross_product(patch0, hinfo, other, **kwargs):
    patch1 = other.level(hinfo.ilvl, hinfo.time)[hinfo.patch_idx]
    ref_name = next(iter(patch0.patch_datas.keys()))

    dset_x = patch0["y"][:] * patch1["z"][:] - patch0["z"][:] * patch1["y"][:]
    dset_y = patch0["z"][:] * patch1["x"][:] - patch0["x"][:] * patch1["z"][:]
    dset_z = patch0["x"][:] * patch1["y"][:] - patch0["y"][:] * patch1["x"][:]

    return (
        {"name": "x", "data": patch0[ref_name].copy_as(dset_x)},
        {"name": "y", "data": patch0[ref_name].copy_as(dset_y)},
        {"name": "z", "data": patch0[ref_name].copy_as(dset_z)},
    )


def _compute_grad(patch, **kwargs):
    ref_name = next(iter(patch.patch_datas.keys()))
    ndim = patch[ref_name].box.ndim
    nb_ghosts = kwargs["nb_ghosts"]
    ds = patch[ref_name].dataset

    ds_shape = list(ds.shape)

    ds_x = np.full(ds_shape, np.nan)
    ds_y = np.full(ds_shape, np.nan)
    ds_z = np.full(ds_shape, np.nan)

    grad_ds = np.gradient(ds)
    select = tuple([slice(nb_ghosts, -nb_ghosts) for _ in range(ndim)])
    if ndim == 2:
        ds_x[select] = np.asarray(grad_ds[0][select])
        ds_y[select] = np.asarray(grad_ds[1][select])
        ds_z[select].fill(0.0)  # TODO at 2D, gradient is null in z dir

    else:
        raise RuntimeError("dimension not yet implemented")

    return (
        {"name": "x", "data": patch[ref_name].copy_as(ds_x)},
        {"name": "y", "data": patch[ref_name].copy_as(ds_y)},
        {"name": "z", "data": patch[ref_name].copy_as(ds_z)},
    )


def dot(hier0, hier1, **kwargs):
    from pyphare.pharesee.hierarchy import scalarfield, vectorfield
    from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

    if not all(isinstance(hier, vectorfield.VectorField) for hier in [hier0, hier1]):
        raise RuntimeError("type of hierarchy not yet considered")

    return scalarfield.ScalarField.FROM(
        compute_hier_from(
            _compute_dot_product, hier0, **vectorfield.copy_kwargs(hier0, other=hier1)
        )
    )


def cross(hier0, hier1, **kwargs):
    from pyphare.pharesee.hierarchy import vectorfield
    from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

    if not all(isinstance(hier, vectorfield.VectorField) for hier in [hier0, hier1]):
        raise RuntimeError("type of hierarchy not yet considered")

    return vectorfield.VectorField.FROM(
        compute_hier_from(_compute_cross_product, hier0, other=hier1)
    )


def sqrt(hier, **kwargs):
    from pyphare.pharesee.hierarchy import scalarfield
    from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

    return scalarfield.ScalarField.FROM(
        compute_hier_from(_compute_sqrt, hier, **scalarfield.copy_kwargs(hier))
    )


def modulus(hier):
    from pyphare.pharesee.hierarchy import vectorfield

    assert isinstance(hier, vectorfield.VectorField)

    return sqrt(dot(hier, hier))


def grad(hier, **kwargs):
    from pyphare.pharesee.hierarchy import scalarfield, vectorfield
    from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

    assert isinstance(hier, scalarfield.ScalarField)
    nb_ghosts = list(hier.level(0).patches[0].patch_datas.values())[0].ghosts_nbr[0]
    h = compute_hier_from(_compute_grad, hier, nb_ghosts=nb_ghosts)

    return vectorfield.VectorField.FROM(h)
