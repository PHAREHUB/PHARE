#
#
#

import numpy as np
from copy import deepcopy


from pyphare.core.box import Box
from .hierarchy import PatchHierarchy

from .patch import Patch
from . import hierarchy_compute as hc
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
        cls = type(input)
        if cls is Box or cls is slice:
            return get_interpolated_selection_from(self, input)
        if input in self.__dict__:
            return self.__dict__[input]
        raise ValueError("AnyTensorField.__getitem__ cannot handle input", input)


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
                "TensorField * TensorField is ambiguous, use pyphare.core.operators.dot or .prod"
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


def GetDomainSize(hier, **kwargs):
    root_cell_width = hier.level(0).cell_width
    domain_box = hier.domain_box
    return (domain_box.upper + 1) * root_cell_width


def GetDl(hier, time, level="finest"):
    level = hier.finest_level(time) if level == "finest" else level
    return hier.level(level).cell_width


def resolve_mutators_in(hier):
    for time in hier.time_hier:
        for ilvl, lvl in hier.levels(time).items():
            for patch in lvl:
                for qty in hier.quantities():
                    patch[qty].dataset = patch[qty][:]
    return hier


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

        pdata = patch0[qty].copy_as(
            layout=patch0.layout,
            name=patch0[qty].name,
            data=interpolator(*mesh),
            ghosts_nbr=nbrGhosts,
            centering=patch0[qty].centerings,
        )
        new_patch0.patch_datas[qty] = pdata

    level0.patches = [new_patch0]
    cpy.time_hier[time] = {0: level0}

    return cpy
