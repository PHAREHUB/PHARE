#
#
#


from .hierarchy import PatchHierarchy
from . import hierarchy_compute as hc
from .hierarchy_utils import compute_hier_from


class AnyTensorField(PatchHierarchy):
    def __init__(self, hier):
        super().__init__(
            hier.patch_levels,
            hier.domain_box,
            hier.refinement_ratio,
            hier.times(),
            hier.data_files,
        )


class TensorField(AnyTensorField):
    def __init__(self, hier):
        self.names = ["xx", "xy", "xz", "yy", "yz", "zz"]
        super().__init__(
            compute_hier_from(hc.compute_rename, hier, new_names=self.names)
        )

    def __mul__(self, other):
        if type(other) is TensorField:
            raise ValueError(
                "TensorField * TensorField is ambiguous, use pyphare.core.operators.dot or .prod"
            )
        return TensorField(compute_hier_from(hc.compute_mul, self, other=other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return TensorField(compute_hier_from(hc.compute_add, self, other=other))

    def __sub__(self, other):
        return TensorField(compute_hier_from(hc.compute_sub, self, other=other))

    def __truediv__(self, other):
        return TensorField(compute_hier_from(hc.compute_truediv, self, other=other))
