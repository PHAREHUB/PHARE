import numpy as np

from pyphare.core import box as boxm
from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from
from pyphare.pharesee.hierarchy.patchdata import align_ghosts


def _compute_dot_product(patch0, accessor, operand, **kwargs):
    ref_name = next(iter(patch0.patch_datas.keys()))
    pd0, pd1 = _aligned_xyz(patch0, operand[accessor])
    dset = (
        pd0["x"][:] * pd1["x"][:]
        + pd0["y"][:] * pd1["y"][:]
        + pd0["z"][:] * pd1["z"][:]
    )

    return (pd0[ref_name].copy_as(dset, name="value"),)


def _compute_sqrt(patch, **kwargs):
    ref_name = next(iter(patch.patch_datas.keys()))

    dset = np.sqrt(patch["value"][:])

    return (patch[ref_name].copy_as(dset, name="value"),)


def _compute_cross_product(patch0, accessor, operand, **kwargs):
    ref_name = next(iter(patch0.patch_datas.keys()))
    pd0, pd1 = _aligned_xyz(patch0, operand[accessor])

    dset_x = pd0["y"][:] * pd1["z"][:] - pd0["z"][:] * pd1["y"][:]
    dset_y = pd0["z"][:] * pd1["x"][:] - pd0["x"][:] * pd1["z"][:]
    dset_z = pd0["x"][:] * pd1["y"][:] - pd0["y"][:] * pd1["x"][:]

    return (
        pd0[ref_name].copy_as(dset_x, name="x"),
        pd0[ref_name].copy_as(dset_y, name="y"),
        pd0[ref_name].copy_as(dset_z, name="z"),
    )


def _compute_grad(patch, **kwargs):
    ndim = patch["value"].box.ndim
    ds = patch["value"].dataset

    ds_shape = list(ds.shape)

    ds_x = np.full(ds_shape, np.nan)
    ds_y = np.full(ds_shape, np.nan)
    ds_z = np.full(ds_shape, np.nan)

    grad_ds = np.gradient(ds)
    pd = patch["value"]
    lbox = boxm.amr_to_local(pd.box, pd.ghost_box)
    lbox.upper += pd.primal_directions()
    if ndim == 2:
        boxm.DataSelector(ds_x)[lbox] = boxm.select(grad_ds[0], lbox)
        boxm.DataSelector(ds_y)[lbox] = boxm.select(grad_ds[1], lbox)
        boxm.DataSelector(ds_z)[lbox] = 0.0  # TODO at 2D, gradient is null in z dir

    else:
        raise RuntimeError("dimension not yet implemented")

    return (
        patch["value"].copy_as(ds_x, name="x"),
        patch["value"].copy_as(ds_y, name="y"),
        patch["value"].copy_as(ds_z, name="z"),
    )


def dot(hier0, hier1, **kwargs):
    if not all(isinstance(hier, VectorField) for hier in [hier0, hier1]):
        raise TypeError("type of hierarchy not yet considered")

    return ScalarField.FROM(
        compute_hier_from(_compute_dot_product, hier0, operand=hier1)
    )


def cross(hier0, hier1, **kwargs):
    if not all(isinstance(hier, VectorField) for hier in [hier0, hier1]):
        raise TypeError("type of hierarchy not yet considered")

    return VectorField.FROM(
        compute_hier_from(_compute_cross_product, hier0, operand=hier1)
    )


def sqrt(hier, **kwargs):
    return ScalarField.FROM(compute_hier_from(_compute_sqrt, hier))


def modulus(hier):
    assert isinstance(hier, VectorField)

    return sqrt(dot(hier, hier))


def grad(hier, **kwargs):
    assert isinstance(hier, ScalarField)
    h = compute_hier_from(_compute_grad, hier)

    return VectorField.FROM(h)


def _aligned_xyz(patch0, patch1):
    """
    vector components with different ghosts can only be combined on their domain
    all six are aligned together so every product has the same shape
    """
    xyz = ("x", "y", "z")
    aligned = align_ghosts([patch0[c] for c in xyz] + [patch1[c] for c in xyz])
    return dict(zip(xyz, aligned[:3])), dict(zip(xyz, aligned[3:]))
