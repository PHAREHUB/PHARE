#
#
#


import numpy as np
from copy import deepcopy


def GetDomainSize(hier, **kwargs):
    root_cell_width = hier.level(0).cell_width
    domain_box = hier.domain_box
    return (domain_box.upper + 1) * root_cell_width


def GetDl(hier, time, level="finest"):
    level = hier.finest_level(time) if level == "finest" else level
    return hier.level(level, time).cell_width


def GetTime(hier):
    times = hier.times()
    if len(times) > 1:
        raise ValueError("Error: more than 1 time found in hierarchy!")
    return times[0]


def GetFinest(hier, time, qty=None, interp="nearest"):
    from pyphare.pharesee.hierarchy import uniformgrid as uniform

    if not hier.ephemerals:
        hier.ephemerals = {}
    if time not in hier.ephemerals:
        hier.ephemerals[time] = {}

    finest = "finest"
    if finest in hier.ephemerals[time]:
        if not qty:
            return hier.ephemerals[time][finest]
        if qty and qty in hier.ephemerals[time]:
            return hier.ephemerals[time][qty]
    else:
        hier.ephemerals[time][finest] = uniform.UniformGrids({})

    selection_box = None  # todo
    datas = get_interpolated_selection_from(hier, selection_box, qty, interp)
    for k, v in datas.items():
        hier.ephemerals[time][finest][k] = v
    return hier.ephemerals[time][finest]


def get_interpolated_selection_from(hier, selection, quantity=None, interp="nearest"):
    """selection to become selection box or slice etc"""

    times = hier.times()
    if len(times) > 1:
        raise ValueError("Error: interpolation does not support multiple times")

    from pyphare.pharesee.hierarchy import uniformgrid as uniform
    from ..interpolation import make_interpolator, flat_finest_field

    time = times[0]
    if 0 not in hier.levels(time):
        raise ValueError("Error: interpolation only supports coarse timesteps")

    dl = GetDl(hier, time)
    domain = GetDomainSize(hier)

    cpy = deepcopy(hier)
    levels = cpy.levels(time)
    level0 = levels[0]
    patch0 = level0.patches[0]
    patch0.layout.box = hier.level_domain_box(len(levels) - 1)  # hax todo etc
    patch0.layout.dl = dl
    patch0.layout.origin = dl * 0
    nbrGhosts = list(hier.level(0).patches[0].patch_datas.values())[0].ghosts_nbr

    datas = {}

    for qty in [quantity] if quantity else hier.quantities():
        flat = flat_finest_field(hier, qty, time=time)
        interpolator, finest_coords = make_interpolator(flat, interp, domain, dl)

        mesh = np.meshgrid(*finest_coords, indexing="ij")

        datas[qty] = uniform.UniformGrid(
            layout=patch0.layout,
            field_name=patch0[qty].name,
            data=interpolator(*mesh),
            ghosts_nbr=nbrGhosts,
            centering=patch0[qty].centerings,
        )

    return datas
