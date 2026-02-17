#
#
#


from . import tensorfield
from . import hierarchy_compute as hc
from .hierarchy_utils import compute_hier_from


class ScalarField(tensorfield.AnyTensorField):
    def __init__(self, hier):
        super().__init__(
            compute_hier_from(hc.compute_rename, hier, new_names=("value",))
        )
        # renamed_hier =

    def __add__(self, other):
        return ScalarField(compute_hier_from(hc.compute_add, self, other=other))

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return ScalarField(compute_hier_from(hc.compute_sub, self, other=other))

    def __mul__(self, other):
        return ScalarField(compute_hier_from(hc.compute_mul, self, other=other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return ScalarField(compute_hier_from(hc.compute_truediv, self, other=other))

    def __rtruediv__(self, other):
        return ScalarField(compute_hier_from(hc.compute_rtruediv, self, other=other))

    def __neg__(self):
        return self * -1
