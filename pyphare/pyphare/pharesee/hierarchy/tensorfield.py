#
#
#

import numpy as np
from copy import deepcopy


from pyphare.core.box import Box
from .hierarchy import PatchHierarchy

from .patch import Patch
from .patchdata import FieldData
from . import hierarchy_compute as hc
from . import hierarchy_utils as hootils


class AnyTensorField(PatchHierarchy):
    def __init__(self, hier):
        super().__init__(
            hier.patch_levels,
            hier.domain_box,
            hier.refinement_ratio,
            hier.times(),
            hier.data_files,
        )

    def __getitem__(self, input):
        cls = type(input)
        if cls is Box or cls is slice:
            return get_interpolated_selection_from(self, input)
        return self.__dict__[input]


class TensorField(AnyTensorField):
    def __init__(self, hier):
        self.names = ["xx", "xy", "xz", "yy", "yz", "zz"]
        super().__init__(
            hootils.compute_hier_from(hc.compute_rename, hier, new_names=self.names)
        )

    def __mul__(self, other):
        if type(other) is TensorField:
            raise ValueError(
                "TensorField * TensorField is ambiguous, use pyphare.core.operators.dot or .prod"
            )
        return TensorField(hootils.compute_hier_from(hc.compute_mul, self, other=other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        return TensorField(hootils.compute_hier_from(hc.compute_add, self, other=other))

    def __sub__(self, other):
        return TensorField(hootils.compute_hier_from(hc.compute_sub, self, other=other))

    def __truediv__(self, other):
        return TensorField(
            hootils.compute_hier_from(hc.compute_truediv, self, other=other)
        )


def GetDomainSize(hier, **kwargs):
    root_cell_width = hier.level(0).cell_width
    domain_box = hier.domain_box
    return (domain_box.upper + 1) * root_cell_width


def GetDl(hier, time, level="finest"):
    level = hier.finest_level(time) if level == "finest" else level
    return hier.level(level).cell_width


def get_interpolated_selection_from(hier: AnyTensorField, input, interp="nearest"):
    if type(hier) is PatchHierarchy:
        raise ValueError("PatchHierarchy not supported, must be AnyTensorField")

    times = hier.times()
    if len(hier.times()) > 1:
        raise ValueError(
            "AnyTensorField interpoloation does not support multiple times"
        )

    from pyphare.pharesee.run import utils as rutils

    time = times[0]

    dl = GetDl(hier, time)
    domain = GetDomainSize(hier)

    cpy = deepcopy(hier)
    levels = cpy.levels(time)
    level0 = levels[0]
    patch0 = level0.patches[0]

    patch0.layout.box = hier.level_domain_box(len(levels) - 1)  # hax todo etc

    new_patch0 = Patch({}, patch0.id, patch0.layout)

    nbrGhosts = list(hier.level(0).patches[0].patch_datas.values())[0].ghosts_nbr
    for qty in hier.quantities():
        data, coords = hootils.flat_finest_field(hier, qty, time=time)
        interpolator, finest_coords = rutils.make_interpolator(
            data, coords, interp, domain, dl, qty, nbrGhosts
        )

        mesh = np.meshgrid(*finest_coords, indexing="ij")

        pdata = FieldData(
            patch0.layout,
            patch0[qty].name,
            interpolator(*mesh),
            ghosts_nbr=patch0[qty].ghosts_nbr,
            centering=patch0[qty].centerings,
        )
        new_patch0.patch_datas[qty] = pdata

    level0.patches = [new_patch0]
    cpy.time_hier[time] = {0: level0}

    return cpy
