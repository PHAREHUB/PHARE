import numpy as np

from pyphare.core import box as boxm
from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from
from pyphare.pharesee.hierarchy.hierarchy_utils import rename


def _compute_dot_product(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset = (
        patch_datas["left_x"][:] * patch_datas["right_x"][:]
        + patch_datas["left_y"][:] * patch_datas["right_y"][:]
        + patch_datas["left_z"][:] * patch_datas["right_z"][:]
    )

    return (
        {"name": "value", "data": dset, "centering": patch_datas[ref_name].centerings},
    )


def _compute_sqrt(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset = np.sqrt(patch_datas["value"][:])

    return (
        {"name": "value", "data": dset, "centering": patch_datas[ref_name].centerings},
    )


def _compute_cross_product(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset_x = (
        patch_datas["left_y"][:] * patch_datas["right_z"][:]
        - patch_datas["left_z"][:] * patch_datas["right_y"][:]
    )
    dset_y = (
        patch_datas["left_z"][:] * patch_datas["right_x"][:]
        - patch_datas["left_x"][:] * patch_datas["right_z"][:]
    )
    dset_z = (
        patch_datas["left_x"][:] * patch_datas["right_y"][:]
        - patch_datas["left_y"][:] * patch_datas["right_x"][:]
    )

    return (
        {"name": "x", "data": dset_x, "centering": patch_datas[ref_name].centerings},
        {"name": "y", "data": dset_y, "centering": patch_datas[ref_name].centerings},
        {"name": "z", "data": dset_z, "centering": patch_datas[ref_name].centerings},
    )


def _compute_grad(patch_data, **kwargs):
    ndim = patch_data["value"].box.ndim
    ds = patch_data["value"].dataset

    ds_shape = list(ds.shape)

    ds_x = np.full(ds_shape, np.nan)
    ds_y = np.full(ds_shape, np.nan)
    ds_z = np.full(ds_shape, np.nan)

    grad_ds = np.gradient(ds)
    pd = patch_data["value"]
    lbox = boxm.amr_to_local(pd.box, pd.ghost_box)
    lbox.upper += pd.primal_directions()
    if ndim == 2:
        boxm.DataSelector(ds_x)[lbox] = boxm.select(grad_ds[0], lbox)
        boxm.DataSelector(ds_y)[lbox] = boxm.select(grad_ds[1], lbox)
        boxm.DataSelector(ds_z)[lbox] = 0.0  # TODO at 2D, gradient is null in z dir

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
    h = compute_hier_from(_compute_grad, hier)

    return VectorField(h)
