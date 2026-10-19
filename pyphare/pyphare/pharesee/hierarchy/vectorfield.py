from . import hierarchy_compute as hc
from . import tensorfield
from .hierarchy_utils import compute_hier_from


class VectorField(tensorfield.AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = ["x", "y", "z"]

    @classmethod
    def FROM(cls, hier):
        renamed_hier = compute_hier_from(
            hc.compute_rename, hier, new_names=("x", "y", "z")
        )
        return super().FROM(cls, renamed_hier)

    def __mul__(self, other):
        if type(other) is VectorField:
            raise ValueError(
                "VectorField * VectorField is ambiguous, use pyphare.core.operators.dot or .prod"
            )
        return VectorField.FROM(compute_hier_from(hc.compute_mul, self, operand=other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return VectorField.FROM(compute_hier_from(hc.compute_add, self, operand=other))

    def __sub__(self, other):
        return VectorField.FROM(compute_hier_from(hc.compute_sub, self, operand=other))

    def __truediv__(self, other):
        return VectorField.FROM(
            compute_hier_from(hc.compute_truediv, self, operand=other)
        )
