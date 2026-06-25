import numpy as np

from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from


def _compute_dot_product(patch0, hinfo, other, **kwargs):
    patch1 = other.level(hinfo.ilvl, hinfo.time)[hinfo.patch_idx]
    ref_name = next(iter(patch0.patch_datas.keys()))
    dset = (
        patch0["x"][:] * patch1["x"][:]
        + patch0["y"][:] * patch1["y"][:]
        + patch0["z"][:] * patch1["z"][:]
    )

    return (patch0[ref_name].copy_as(dset, name="value"),)


def _compute_sqrt(patch, **kwargs):
    ref_name = next(iter(patch.patch_datas.keys()))

    dset = np.sqrt(patch["value"][:])

    return (patch[ref_name].copy_as(dset, name="value"),)


def _compute_cross_product(patch0, hinfo, other, **kwargs):
    patch1 = other.level(hinfo.ilvl, hinfo.time)[hinfo.patch_idx]
    ref_name = next(iter(patch0.patch_datas.keys()))

    dset_x = patch0["y"][:] * patch1["z"][:] - patch0["z"][:] * patch1["y"][:]
    dset_y = patch0["z"][:] * patch1["x"][:] - patch0["x"][:] * patch1["z"][:]
    dset_z = patch0["x"][:] * patch1["y"][:] - patch0["y"][:] * patch1["x"][:]

    return (
        patch0[ref_name].copy_as(dset_x, name="x"),
        patch0[ref_name].copy_as(dset_y, name="y"),
        patch0[ref_name].copy_as(dset_z, name="z"),
    )


def _compute_grad(patch, **kwargs):
    ndim = patch["value"].box.ndim
    nb_ghosts = kwargs["nb_ghosts"]
    ds = patch["value"].dataset

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
        patch["value"].copy_as(ds_x, name="x"),
        patch["value"].copy_as(ds_y, name="y"),
        patch["value"].copy_as(ds_z, name="z"),
    )


def dot(hier0, hier1, **kwargs):
    if not all(isinstance(hier, VectorField) for hier in [hier0, hier1]):
        raise RuntimeError("type of hierarchy not yet considered")

    return ScalarField(compute_hier_from(_compute_dot_product, hier0, other=hier1))


def cross(hier0, hier1, **kwargs):
    if not all(isinstance(hier, VectorField) for hier in [hier0, hier1]):
        raise RuntimeError("type of hierarchy not yet considered")

    return VectorField(compute_hier_from(_compute_cross_product, hier0, other=hier1))


def sqrt(hier, **kwargs):
    return ScalarField(compute_hier_from(_compute_sqrt, hier))


def modulus(hier):
    assert isinstance(hier, VectorField)

    return sqrt(dot(hier, hier))


def grad(hier, **kwargs):
    assert isinstance(hier, ScalarField)
    nb_ghosts = list(hier.level(0).patches[0].patch_datas.values())[0].ghosts_nbr[0]
    h = compute_hier_from(_compute_grad, hier, nb_ghosts=nb_ghosts)

    return VectorField(h)
