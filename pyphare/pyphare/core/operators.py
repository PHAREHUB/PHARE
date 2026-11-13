import numpy as np

from pyphare.pharesee.hierarchy import PatchHierarchy, ScalarField, VectorField
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from, rename

# ---------------------------------------------------------------------------- new fields
# Fields of pyphare.pharesee.amr (new reader) are handled by amr.ops at the top of each
# function below; the code after that dispatch is master's, for legacy hierarchies.
# The test is done first because amr fields are registered as virtual subclasses of
# ScalarField / VectorField (isinstance would send them to the legacy code).


def _is_amr(obj):
    """True for a field of pyphare.pharesee.amr (detected by module, nothing imported)."""
    return any(
        getattr(c, "__module__", "").startswith("pyphare.pharesee.amr")
        for c in type(obj).__mro__
    )


def _amr_apply(what, kernel, operands, expected):
    """``kernel`` (an amr.ops function on Snapshots) time by time -> new field.

    expected: required component names of every operand (("x","y","z") or ("value",)).
    TypeError when an operand is a legacy hierarchy (with advice) or of the wrong kind.
    """
    from pyphare.pharesee.amr import _protocol

    for op in operands:
        if not _is_amr(op):
            if isinstance(op, PatchHierarchy):
                raise TypeError(_protocol._LEGACY_MSG.format(what=what))
            raise TypeError(f"{what}: type of hierarchy not yet considered ({type(op).__name__})")
        if tuple(op.names) != expected:
            kind = "VectorField" if expected == ("x", "y", "z") else "ScalarField"
            raise TypeError(
                f"{what}: type of hierarchy not yet considered (expected {kind} "
                f"operands, got components {tuple(op.names)})"
            )
    return _protocol._map(kernel, list(operands), what)


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
    nb_ghosts = kwargs["nb_ghosts"]
    ds = patch_data["value"].dataset

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
        {"name": "x", "data": ds_x, "centering": patch_data["value"].centerings},
        {"name": "y", "data": ds_y, "centering": patch_data["value"].centerings},
        {"name": "z", "data": ds_z, "centering": patch_data["value"].centerings},
    )


def dot(hier_left, hier_right, **kwargs):
    """Scalar product of two VectorFields -> ScalarField (new fields: amr.ops.arith.dot)."""
    if _is_amr(hier_left) or _is_amr(hier_right):
        from pyphare.pharesee.amr.ops import arith

        return _amr_apply("dot", arith.dot, (hier_left, hier_right), ("x", "y", "z"))
    if isinstance(hier_left, VectorField) and isinstance(hier_right, VectorField):
        names_left = ["left_x", "left_y", "left_z"]
        names_right = ["right_x", "right_y", "right_z"]

    else:
        raise TypeError("type of hierarchy not yet considered")

    hl = rename(hier_left, names_left)
    hr = rename(hier_right, names_right)

    h = compute_hier_from(
        _compute_dot_product,
        (hl, hr),
    )

    return ScalarField.FROM(h)


def cross(hier_left, hier_right, **kwargs):
    """Vector product of two VectorFields -> VectorField (new fields: amr.ops.arith.cross)."""
    if _is_amr(hier_left) or _is_amr(hier_right):
        from pyphare.pharesee.amr.ops import arith

        return _amr_apply("cross", arith.cross, (hier_left, hier_right), ("x", "y", "z"))
    if isinstance(hier_left, VectorField) and isinstance(hier_right, VectorField):
        names_left = ["left_x", "left_y", "left_z"]
        names_right = ["right_x", "right_y", "right_z"]

    else:
        raise TypeError("type of hierarchy not yet considered")

    hl = rename(hier_left, names_left)
    hr = rename(hier_right, names_right)

    h = compute_hier_from(
        _compute_cross_product,
        (hl, hr),
    )

    return VectorField.FROM(h)


def sqrt(hier, **kwargs):
    """Square root of a ScalarField -> ScalarField (new fields: np.sqrt, elementwise)."""
    if _is_amr(hier):
        from pyphare.pharesee.amr.ops import arith

        return _amr_apply("sqrt", lambda s: arith.unary(np.sqrt, s), (hier,), ("value",))
    h = compute_hier_from(
        _compute_sqrt,
        hier,
    )

    return ScalarField.FROM(h)


def modulus(hier):
    """Norm of a VectorField -> ScalarField (new fields: amr.ops.arith.modulus)."""
    if _is_amr(hier):
        from pyphare.pharesee.amr.ops import arith

        return _amr_apply("modulus", arith.modulus, (hier,), ("x", "y", "z"))
    assert isinstance(hier, VectorField)

    return sqrt(dot(hier, hier))


def grad(hier, **kwargs):
    """Gradient of a ScalarField -> VectorField.

    New fields: amr.ops.deriv.grad, level-aware, divided by the mesh size (master's legacy
    code below divides by nothing, issue #1446), in 1D, 2D and 3D.
    """
    if _is_amr(hier):
        from pyphare.pharesee.amr.ops import deriv

        return _amr_apply("grad", deriv.grad, (hier,), ("value",))
    assert isinstance(hier, ScalarField)
    nb_ghosts = next(iter(hier.level(0).patches[0].patch_datas.values())).ghosts_nbr[0]
    h = compute_hier_from(_compute_grad, hier, nb_ghosts=nb_ghosts)

    return VectorField.FROM(h)
