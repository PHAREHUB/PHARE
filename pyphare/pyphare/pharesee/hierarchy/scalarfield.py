#
#
#


from . import tensorfield
from . import compute as hc
from . import hierarchy_utils as hootils


class ScalarField(tensorfield.AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def FROM(cls, hier):
        return super().FROM(cls, hier)

    def __add__(self, other):
        return ScalarField.FROM(
            hootils.compute_hier_from(hc.compute_add, self, **copy_kwargs(self, other))
        )

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return ScalarField.FROM(
            hootils.compute_hier_from(hc.compute_sub, self, **copy_kwargs(self, other))
        )

    def __mul__(self, other):
        return ScalarField.FROM(
            hootils.compute_hier_from(hc.compute_mul, self, **copy_kwargs(self, other))
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return ScalarField.FROM(
            hootils.compute_hier_from(
                hc.compute_truediv, self, **copy_kwargs(self, other)
            )
        )

    def __rtruediv__(self, other):
        return ScalarField.FROM(
            hootils.compute_hier_from(
                hc.compute_rtruediv, self, **copy_kwargs(self, other)
            )
        )

    def __neg__(self):
        return self * -1


def copy_kwargs(scalar_field, other=None):
    if other:
        return {"other": other, "key_map": key_map(scalar_field, other)}
    return {"key_map": key_map(scalar_field)}


def key_map(scalar_field, other=None):
    if len(scalar_field.quantities()) != 1:
        raise ValueError("Invalid ScalarField 1!")

    both_fields = other and type(scalar_field) is type(other)
    if both_fields and len(other.quantities()) != 1:
        raise ValueError("Invalid ScalarField 2!")

    def _(field):
        return {field.quantities()[i]: "value" for i in range(len(field.quantities()))}

    this = _(scalar_field)
    that = _(other) if both_fields else {}
    return {**this, **that}
