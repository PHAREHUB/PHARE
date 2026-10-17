def GetDomainSize(hier):
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


def GetFinest(hier, time=None, qty=None, interp="nearest"):
    """
    Returns UniformGrids of qty, or of all quantities if qty is None
     interpolated grids are cached per time and interp
    """
    from pyphare.pharesee.hierarchy import uniformgrid as uniform
    from pyphare.pharesee.hierarchy.hierarchy import format_timestamp
    from pyphare.pharesee.run import utils as rutils

    time = format_timestamp(GetTime(hier) if time is None else time)
    if hier.ephemerals is None:
        hier.ephemerals = {}
    cache = hier.ephemerals.setdefault(time, {})
    cache = cache.setdefault(("finest", interp), uniform.UniformGrids({}))

    qties = [qty] if qty else hier.quantities()
    for q in [q for q in qties if q not in cache]:
        for k, v in rutils.interpolate_hierarchy(
            hier, quantity=q, interp=interp
        ).items():
            cache[k] = v
    return uniform.UniformGrids({q: cache[q] for q in qties})
