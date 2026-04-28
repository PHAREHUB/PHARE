#
#
#

from .hierarchy_utils import (
    compute_hier_from,
    compute_rename,
    rename,
    _compute_mul,
    _compute_add,
    _compute_sub,
    _compute_truediv,
    _compute_scalardiv,
)
from .scalarfield import ScalarField

from . import tensorfield


class VectorField(tensorfield.AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = ["x", "y", "z"]

    @classmethod
    def FROM(cls, hier):
        renamed_hier = compute_hier_from(
            compute_rename, hier, new_names=("x", "y", "z")
        )
        return super().FROM(cls, renamed_hier)

    def __mul__(self, other):
        assert isinstance(other, (int, float))
        h = compute_hier_from(_compute_mul, self, names=["x", "y", "z"], other=other)
        return VectorField.FROM(h)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        names_self_kept = self.quantities()
        names_other_kept = other.quantities()

        if isinstance(other, VectorField):
            names_self = ["self_x", "self_y", "self_z"]
            names_other = ["other_x", "other_y", "other_z"]
        else:
            raise RuntimeError("type of hierarchy not yet considered")

        h_self = rename(self, names_self)
        h_other = rename(other, names_other)

        h = compute_hier_from(
            _compute_add,
            (h_self, h_other),
        )

        self = rename(h_self, names_self_kept)  # needed ?
        other = rename(h_other, names_other_kept)

        return VectorField.FROM(h)

    def __sub__(self, other):
        names_self_kept = self.quantities()
        names_other_kept = other.quantities()

        if isinstance(other, VectorField):
            names_self = ["self_x", "self_y", "self_z"]
            names_other = ["other_x", "other_y", "other_z"]
        else:
            raise RuntimeError("type of hierarchy not yet considered")

        h_self = rename(self, names_self)
        h_other = rename(other, names_other)

        h = compute_hier_from(
            _compute_sub,
            (h_self, h_other),
        )

        self = rename(h_self, names_self_kept)
        other = rename(h_other, names_other_kept)

        return VectorField.FROM(h)

    def __truediv__(self, other):
        if not isinstance(other, (ScalarField, int, float)):
            raise RuntimeError("type of operand not considered")

        if isinstance(other, ScalarField):
            return VectorField.FROM(
                compute_hier_from(
                    _compute_truediv, (self, other), res_names=("x", "y", "z")
                )
            )
        elif isinstance(other, (int, float)):
            return VectorField.FROM(
                compute_hier_from(
                    _compute_scalardiv, (self,), res_names=("x", "y", "z"), scalar=other
                )
            )
