from . import hierarchy_compute as hc
from .hierarchy import PatchHierarchy
from .hierarchy_utils import compute_hier_from


class ScalarField(PatchHierarchy):
    def __init__(self, hier):
        renamed_hier = compute_hier_from(hc.compute_rename, hier, new_names=("value",))
        patch_levels = renamed_hier.patch_levels
        domain_box = renamed_hier.domain_box
        refinement_ratio = renamed_hier.refinement_ratio
        data_files = renamed_hier.data_files

        super().__init__(
            patch_levels, domain_box, refinement_ratio, renamed_hier.times(), data_files
        )

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
