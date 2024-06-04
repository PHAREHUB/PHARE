from .hierarchy import PatchHierarchy
from .fromh5 import hierarchy_fromh5
from .fromsim import hierarchy_from_sim


def hierarchy_from(
    simulator=None, qty=None, pop="", h5_filename=None, times=None, hier=None, **kwargs
):
    from .fromh5 import hierarchy_fromh5
    from .fromsim import hierarchy_from_sim

    """
    this function reads an HDF5 PHARE file and returns a PatchHierarchy from
    which data is accessible.
    if 'time' is None, all times in the file will be read, if a time is given
    then only that time will be read
    if 'hier' is None, then a new hierarchy will be created, if not then the
    given hierarchy 'hier' will be filled.

    The function fails if the data is already in hierarchy
    """

    if simulator is not None and h5_filename is not None:
        raise ValueError("cannot pass both a simulator and a h5 file")

    if h5_filename is not None:
        return hierarchy_fromh5(h5_filename, times, hier, **kwargs)

    if simulator is not None and qty is not None:
        return hierarchy_from_sim(simulator, qty, pop=pop)

    raise ValueError("can't make hierarchy")
