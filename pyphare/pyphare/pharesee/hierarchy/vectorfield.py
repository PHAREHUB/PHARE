#
#
#

from . import tensorfield

from . import compute as hc
from . import hierarchy_utils as hootils


class VectorField(tensorfield.AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = ["x", "y", "z"]

    @classmethod
    def FROM(cls, hier):
        return super().FROM(cls, hier)

    def __mul__(self, other):
        if type(other) is VectorField:
            raise ValueError(
                "VectorField * VectorField is ambiguous, use pyphare.pharesee.hierarchy.compute.dot or .prod"
            )
        return VectorField.FROM(
            hootils.compute_hier_from(hc.compute_mul, self, **copy_kwargs(self, other))
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return VectorField.FROM(
            hootils.compute_hier_from(hc.compute_add, self, **copy_kwargs(self, other))
        )

    def __iadd__(self, other):
        self = VectorField.FROM(
            hootils.compute_hier_from(hc.compute_add, self, **copy_kwargs(self, other))
        )
        return self

    def __sub__(self, other):
        return VectorField.FROM(
            hootils.compute_hier_from(hc.compute_sub, self, **copy_kwargs(self, other))
        )

    def __truediv__(self, other):
        return VectorField.FROM(
            hootils.compute_hier_from(
                hc.compute_truediv, self, **copy_kwargs(self, other)
            )
        )


def copy_kwargs(vector_field, other=None):
    if other:
        return {"other": other, "key_map": key_map(vector_field, other)}
    return {"key_map": key_map(vector_field)}


def key_map(vector_field, other=None):
    if len(vector_field.quantities()) != 3:
        raise ValueError("Invalid VectorField 1!")

    both_vecfields = other and type(vector_field) is type(other)
    if both_vecfields and len(other.quantities()) != 3:
        raise ValueError("Invalid VectorField 2!")

    def _(field):
        return {
            field.quantities()[i]: field.names[i]
            for i in range(len(field.quantities()))
        }

    this = _(vector_field)
    that = _(other) if both_vecfields else {}
    return {**this, **that}
