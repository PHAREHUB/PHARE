import numpy as np

from pyphare.pharesee.hierarchy import ScalarField, VectorField  # TensorField
from pyphare.pharesee.hierarchy import compute_hier_from
from pyphare.pharesee.hierarchy import rename


def _compute_dot_product(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset = (
        patch_datas["left_x"].dataset[:] * patch_datas["right_x"].dataset[:]
        + patch_datas["left_y"].dataset[:] * patch_datas["right_y"].dataset[:]
        + patch_datas["left_z"].dataset[:] * patch_datas["right_z"].dataset[:]
    )

    return (
        {"name": "value", "data": dset, "centering": patch_datas[ref_name].centerings},
    )


def _compute_sqrt(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset = np.sqrt(patch_datas["value"].dataset[:])

    return (
        {"name": "value", "data": dset, "centering": patch_datas[ref_name].centerings},
    )


def _compute_cross_product(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset_x = (
        patch_datas["left_y"].dataset[:] * patch_datas["right_z"].dataset[:]
        - patch_datas["left_z"].dataset[:] * patch_datas["right_y"].dataset[:]
    )
    dset_y = (
        patch_datas["left_z"].dataset[:] * patch_datas["right_x"].dataset[:]
        - patch_datas["left_x"].dataset[:] * patch_datas["right_z"].dataset[:]
    )
    dset_z = (
        patch_datas["left_x"].dataset[:] * patch_datas["right_y"].dataset[:]
        - patch_datas["left_y"].dataset[:] * patch_datas["right_x"].dataset[:]
    )

    return (
        {"name": "x", "data": dset_x, "centering": patch_datas[ref_name].centerings},
        {"name": "y", "data": dset_y, "centering": patch_datas[ref_name].centerings},
        {"name": "z", "data": dset_z, "centering": patch_datas[ref_name].centerings},
    )


def _compute_grad(patch_data, **kwargs):
    ndim = patch_data["value"].box.ndim
    nb_ghosts = kwargs["nb_ghosts"]
    ds = patch_data["value"].dataset

    ds_shape = list(ds.shape)

    ds_x = np.full(ds_shape, np.nan)
    ds_y = np.full(ds_shape, np.nan)
    ds_z = np.full(ds_shape, np.nan)

    grad_ds = np.gradient(ds)

    if ndim == 2:
        ds_x[nb_ghosts:-nb_ghosts, nb_ghosts:-nb_ghosts] = np.asarray(
            grad_ds[0][nb_ghosts:-nb_ghosts, nb_ghosts:-nb_ghosts]
        )
        ds_y[nb_ghosts:-nb_ghosts, nb_ghosts:-nb_ghosts] = np.asarray(
            grad_ds[1][nb_ghosts:-nb_ghosts, nb_ghosts:-nb_ghosts]
        )
        ds_z[nb_ghosts:-nb_ghosts, nb_ghosts:-nb_ghosts].fill(
            0.0
        )  # TODO at 2D, gradient is null in z dir

    else:
        raise RuntimeError("dimension not yet implemented")

    return (
        {"name": "x", "data": ds_x, "centering": patch_data["value"].centerings},
        {"name": "y", "data": ds_y, "centering": patch_data["value"].centerings},
        {"name": "z", "data": ds_z, "centering": patch_data["value"].centerings},
    )


def dot(hier_left, hier_right, **kwargs):
    if isinstance(hier_left, VectorField) and isinstance(hier_right, VectorField):
        names_left = ["left_x", "left_y", "left_z"]
        names_right = ["right_x", "right_y", "right_z"]

    else:
        raise RuntimeError("type of hierarchy not yet considered")

    hl = rename(hier_left, names_left)
    hr = rename(hier_right, names_right)

    h = compute_hier_from(
        _compute_dot_product,
        (hl, hr),
    )

    return ScalarField(h)


def cross(hier_left, hier_right, **kwargs):
    if isinstance(hier_left, VectorField) and isinstance(hier_right, VectorField):
        names_left = ["left_x", "left_y", "left_z"]
        names_right = ["right_x", "right_y", "right_z"]

    else:
        raise RuntimeError("type of hierarchy not yet considered")

    hl = rename(hier_left, names_left)
    hr = rename(hier_right, names_right)

    h = compute_hier_from(
        _compute_cross_product,
        (hl, hr),
    )

    return VectorField(h)


def sqrt(hier, **kwargs):
    h = compute_hier_from(
        _compute_sqrt,
        hier,
    )

    return ScalarField(h)


def modulus(hier):
    assert isinstance(hier, VectorField)

    return sqrt(dot(hier, hier))


def grad(hier, **kwargs):
    assert isinstance(hier, ScalarField)
    nb_ghosts = list(hier.level(0).patches[0].patch_datas.values())[0].ghosts_nbr[0]
    h = compute_hier_from(_compute_grad, hier, nb_ghosts=nb_ghosts)

    return VectorField(h)
