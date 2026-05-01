#
#
#


from . import compute as hc
from .hierarchy import PatchHierarchy
from . import hierarchy_utils as hootils


class AnyTensorField(PatchHierarchy):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def FROM(cls, hier):
        return cls(
            hier.patch_levels,
            hier.domain_box,
            hier.refinement_ratio,
            hier.times(),
            hier.data_files,
        )

    def __getitem__(self, input):
        if input in self.__dict__:
            return self.__dict__[input]
        raise IndexError("AnyTensorField.__getitem__ cannot handle input", input)


class TensorField(AnyTensorField):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.names = ["xx", "xy", "xz", "yy", "yz", "zz"]

    @classmethod
    def FROM(cls, hier):
        return super().FROM(cls, hier)

    def __mul__(self, other):
        if type(other) is TensorField:
            raise ValueError(
                "TensorField * TensorField is ambiguous, use pyphare.pharesee.hierarchy.compute.dot or .prod"
            )
        return TensorField(hootils.compute_hier_from(hc.compute_mul, self, other=other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return TensorField.FROM(
            hootils.compute_hier_from(hc.compute_add, self, other=other)
        )

    def __sub__(self, other):
        return TensorField.FROM(
            hootils.compute_hier_from(hc.compute_sub, self, other=other)
        )

    def __truediv__(self, other):
        return TensorField.FROM(
            hootils.compute_hier_from(hc.compute_truediv, self, other=other)
        )
