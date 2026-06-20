#
#
#


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
    from pyphare.pharesee.run import utils as rutils

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

    grids = rutils.interpolate_hierarchy(hier, quantity=qty, interp=interp)
    for k, v in grids.items():
        hier.ephemerals[time][finest][k] = v
    return hier.ephemerals[time][finest]
