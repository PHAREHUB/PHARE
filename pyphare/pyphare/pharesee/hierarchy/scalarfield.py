from . import hierarchy_compute as hc
from . import tensorfield
from .hierarchy_utils import compute_hier_from


class ScalarField(tensorfield.AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def FROM(cls, hier):
        renamed_hier = compute_hier_from(hc.compute_rename, hier, new_names=("value",))
        return super().FROM(cls, renamed_hier)

    def __add__(self, other):
        return ScalarField.FROM(compute_hier_from(hc.compute_add, self, operand=other))

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return ScalarField.FROM(compute_hier_from(hc.compute_sub, self, operand=other))

    def __mul__(self, other):
        return ScalarField.FROM(compute_hier_from(hc.compute_mul, self, operand=other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return ScalarField.FROM(
            compute_hier_from(hc.compute_truediv, self, operand=other)
        )

    def __rtruediv__(self, other):
        return ScalarField.FROM(
            compute_hier_from(hc.compute_rtruediv, self, operand=other)
        )

    def __neg__(self):
        return self * -1
