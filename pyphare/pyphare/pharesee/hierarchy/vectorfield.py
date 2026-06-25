from . import hierarchy_compute as hc
from .hierarchy import PatchHierarchy
from .hierarchy_utils import compute_hier_from


class VectorField(PatchHierarchy):
    def __init__(self, hier):
        renamed_hier = compute_hier_from(
            hc.compute_rename, hier, new_names=("x", "y", "z")
        )
        patch_levels = renamed_hier.patch_levels
        domain_box = renamed_hier.domain_box
        refinement_ratio = renamed_hier.refinement_ratio
        data_files = renamed_hier.data_files

        self.names = ["x", "y", "z"]

        super().__init__(
            patch_levels, domain_box, refinement_ratio, renamed_hier.times(), data_files
        )

    def __mul__(self, other):
        if type(other) is VectorField:
            raise ValueError(
                "VectorField * VectorField is ambiguous, use pyphare.core.operators.dot or .prod"
            )
        return VectorField(compute_hier_from(hc.compute_mul, self, other=other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return VectorField(compute_hier_from(hc.compute_add, self, other=other))

    def __sub__(self, other):
        return VectorField(compute_hier_from(hc.compute_sub, self, other=other))

    def __truediv__(self, other):
        return VectorField(compute_hier_from(hc.compute_truediv, self, other=other))
