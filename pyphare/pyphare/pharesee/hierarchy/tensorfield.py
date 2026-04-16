#
#
#

import numpy as np
from copy import deepcopy


from pyphare.core.box import Box
from .hierarchy import PatchHierarchy
from .hierarchy import ValueHierarchy
from .hierarchy import IndexHierarchy

from .patch import Patch
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
        if cls is IndexHierarchy:
            return get_from_index_hierarchy(self, input)
        if input in self.__dict__:
            return self.__dict__[input]
        raise IndexError("AnyTensorField.__getitem__ cannot handle input", input)

    def peak_coordinates(self, qty=None, **kwargs):
        return peak_coordinates(self, qty, **kwargs)

    def find_peaks(self, qty=None, **kwargs):
        return find_peaks(self, qty, **kwargs)

    def as_peaks(self, qty=None, **kwargs):
        return as_peaks(self, qty, **kwargs)


def GetDomainSize(hier, **kwargs):
    root_cell_width = hier.level(0).cell_width
    domain_box = hier.domain_box
    return (domain_box.upper + 1) * root_cell_width


def GetDl(hier, time, level="finest"):
    level = hier.finest_level(time) if level == "finest" else level
    return hier.level(level, time).cell_width


def peak_coordinates(hier: AnyTensorField, qty, **kwargs):
    domain = GetDomainSize(hier)
    peaks, _ = find_peaks(hier, qty, **kwargs)
    dl = GetDl(hier, hier.times()[0])
    return peaks * dl[0]  # only works for 1D!


def find_peaks(hier: AnyTensorField, qty, **kwargs):
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html

    if hier.ndim > 1:
        raise ValueError("AnyTensorField::find_peaks only supports 1d simulations")
    if not qty or qty not in hier.quantities():
        raise ValueError(f"AnyTensorField::find_peaks unknown quantity: {qty}")

    from scipy.signal import find_peaks

    times = hier.times()
    if len(times) > 1:
        raise ValueError("AnyTensorField::find_peaks does not support multiple times")

    levels = hier.levels(times[0])
    if len(levels) > 1:
        raise ValueError("AnyTensorField::find_peaks does not support multiple levels")
    if len(next(iter(levels.values())).patches) > 1:
        raise ValueError("AnyTensorField::find_peaks does not support multiple patches")

    patch = next(iter(levels.values())).patches[0]
    for k, v in patch.patch_datas.items():
        if qty == k:
            peaks, props = find_peaks(v.dataset[:], **kwargs)
            if "peak_heights" in props:
                return peaks, props["peak_heights"]
            return peaks, props

    raise ValueError("AnyTensorField::find_peaks found nothing")


def as_peaks(hier: AnyTensorField, qty=None, **kwargs):
    if hier.ndim > 1:
        raise ValueError("AnyTensorField::find_peaks only supports 1d simulations")
    if qty and qty not in hier.quantities():
        raise ValueError(f"AnyTensorField::find_peaks unknown quantitiy: {qty}")

    from scipy.signal import find_peaks

    times = hier.times()
    if len(times) > 1:
        raise ValueError("AnyTensorField::find_peaks does not support multiple times")

    indexes = {}
    for ilvl, lvl in hier.levels(times[0]).items():
        indexes[ilvl] = []
        for patch in lvl:
            peaks = {}
            for k, v in patch.patch_datas.items():
                if qty is None or qty == k:
                    peaks[k] = find_peaks(v.dataset[:], **kwargs)[0]
            indexes[ilvl].append(peaks)
    return IndexHierarchy(times[0], hier, indexes)


def get_from_index_hierarchy(hier, index_hier):
    if hier != index_hier.hier:
        raise ValueError(
            "AnyTensorField::get_from_index_hierarchy cannot operate on two different hierarchies"
        )

    values = {}
    for ilvl, lvl in index_hier.indexes.items():
        patch_level = hier.level(ilvl, index_hier.time)
        values[ilvl] = []
        for idx, indexes in enumerate(lvl):
            data = {}
            patch = patch_level[idx]
            for k, v in indexes.items():
                if k in patch.patch_datas:
                    data[k] = patch[k].dataset[v]
            values[ilvl].append(data)
    return ValueHierarchy(index_hier.time, hier, values)


def get_interpolated_selection_from(hier: AnyTensorField, input, interp="nearest"):
    if type(hier) is PatchHierarchy:
        raise ValueError("PatchHierarchy not supported, must be AnyTensorField")

    times = hier.times()
    if len(times) > 1:
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
    patch0.layout.dl = dl
    patch0.layout.origin = dl * 0

    new_patch0 = Patch({}, patch0.id, patch0.layout)

    nbrGhosts = list(hier.level(0).patches[0].patch_datas.values())[0].ghosts_nbr
    for qty in hier.quantities():
        data, coords = hootils.flat_finest_field(hier, qty, time=time)
        interpolator, finest_coords = rutils.make_interpolator(
            data, coords, interp, domain, dl, qty, nbrGhosts
        )

        mesh = np.meshgrid(*finest_coords, indexing="ij")

        pdata = type(patch0[qty])(
            layout=patch0.layout,
            field_name=patch0[qty].name,
            data=interpolator(*mesh),
            ghosts_nbr=nbrGhosts,
            centering=patch0[qty].centerings,
        )
        new_patch0.patch_datas[qty] = pdata

    level0.patches = [new_patch0]
    cpy.time_hier[time] = {0: level0}

    return cpy
