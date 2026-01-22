import os
import numpy as np

from .patch import Patch
from .patchlevel import PatchLevel
from .patchdata import FieldData, ParticleData
from ..particles import Particles
from .hierarchy import PatchHierarchy
from .hierarchy import format_timestamp
from ...core.box import Box
from ...core.phare_utilities import (
    refinement_ratio,
    none_iterable,
    all_iterables,
)
from ...core.gridlayout import GridLayout
from .hierarchy_utils import field_qties

from pathlib import Path
from pyphare.core.phare_utilities import listify


h5_time_grp_key = "t"
particle_files_patterns = ("domain", "levelGhost")


def get_all_available_quantities_from_h5(filepath, time=0, exclude=["tags"], hier=None):
    time = format_timestamp(time)
    path = Path(filepath)
    for h5 in path.glob("*.h5"):
        if h5.parent == path and h5.stem not in exclude:
            hier = hierarchy_fromh5(str(h5), time, hier)
    return hier


def make_layout(h5_patch_grp, cell_width, interp_order):
    origin = h5_patch_grp.attrs["origin"]
    upper = h5_patch_grp.attrs["upper"]
    lower = h5_patch_grp.attrs["lower"]
    return GridLayout(Box(lower, upper), origin, cell_width, interp_order=interp_order)


def is_pop_fluid_file(basename):
    return (is_particle_file(basename) is False) and "pop" in basename


def particle_dataset_name(basename):
    """
    return "alpha_domain" from ions_pop_alpha_domain.h5
    """
    popname = basename.strip(".h5").split("_")[-2]
    part_type = basename.strip(".h5").split("_")[-1]
    dataset_name = popname + "_" + part_type

    return dataset_name


def is_particle_file(filename):
    return any([pattern in filename for pattern in particle_files_patterns])


def pop_name(basename):
    return basename.strip(".h5").split("_")[2]


def add_to_patchdata(patch_datas, h5_patch_grp, basename, layout):
    """
    adds data in the h5_patch_grp in the given PatchData dict
    returns True if valid h5 patch found
    """

    if is_particle_file(basename):
        v = np.asarray(h5_patch_grp["v"])
        s = v.size
        v = v[:].reshape(int(s / 3), 3)
        nbrParts = v.shape[0]
        dl = np.zeros((nbrParts, layout.ndim))
        for i in range(layout.ndim):
            dl[:, i] = layout.dl[i]

        particles = Particles(
            icells=h5_patch_grp["iCell"],
            deltas=h5_patch_grp["delta"],
            v=v,
            weights=h5_patch_grp["weight"],
            charges=h5_patch_grp["charge"],
            dl=dl,
        )

        pdname = particle_dataset_name(basename)
        if pdname in patch_datas:
            raise ValueError("error - {} already in patchdata".format(pdname))

        patch_datas[pdname] = ParticleData(layout, particles, pop_name(basename))

    else:
        for dataset_name in h5_patch_grp.keys():
            dataset = h5_patch_grp[dataset_name]

            if dataset_name not in field_qties:
                raise RuntimeError(
                    "invalid dataset name : {} is not in {}".format(
                        dataset_name, field_qties
                    )
                )

            pdata = FieldData(layout, field_qties[dataset_name], dataset)

            pdata_name = field_qties[dataset_name]

            if is_pop_fluid_file(basename):
                pdata_name = pop_name(basename) + "_" + pdata_name

            if pdata_name in patch_datas:
                raise ValueError("error - {} already in patchdata".format(dataset_name))

            patch_datas[pdata_name] = pdata

    return True  # valid patch assumed


def patch_has_datasets(h5_patch_grp):
    return len(h5_patch_grp.keys()) > 0


def h5_filename_from(diagInfo):
    # diagInfo.quantity starts with a / , hence   [1:]
    return (diagInfo.quantity + ".h5").replace("/", "_")[1:]


def get_times_from_h5(filepath, as_float=True):
    import h5py  # see doc/conventions.md section 2.1.1

    with h5py.File(filepath, "r") as f:
        times = list(f[h5_time_grp_key].keys())
        if as_float:
            times = np.array(sorted([float(s) for s in times]))
        return times


def create_from_all_times(time, hier):
    return time is None and hier is None


def create_from_times(times, hier):
    return times is not None and hier is None


def load_all_times(time, hier):
    return time is None and hier is not None


def load_one_time(time, hier):
    return time is not None and hier is not None


def patch_levels_from_h5(h5f, time, selection_box=None):
    """
    creates a dictionary of PatchLevels from a given time in a h5 file
    {ilvl: PatchLevel}
    """

    root_cell_width = h5f.attrs["cell_width"]
    interp_order = h5f.attrs["interpOrder"]
    basename = os.path.basename(h5f.filename)

    patch_levels = {}

    for lvl_key, lvl in h5f[h5_time_grp_key][time].items():
        ilvl = int(lvl_key[2:])  # pl1-->1
        lvl_cell_width = root_cell_width / refinement_ratio**ilvl

        patches = []

        for h5_patch in lvl.values():
            lower = h5_patch.attrs["lower"]
            upper = h5_patch.attrs["upper"]
            origin = h5_patch.attrs["origin"]

            patch_box = Box(lower, upper)

            pos_upper = [
                orig + shape * dl
                for orig, shape, dl in zip(origin, patch_box.shape, lvl_cell_width)
            ]
            pos_patch_box = Box(origin, pos_upper)

            intersect = None
            if selection_box is not None:
                intersect = selection_box * pos_patch_box

            if intersect is not None or selection_box is None:
                patch_datas = {}
                layout = make_layout(h5_patch, lvl_cell_width, interp_order)
                if patch_has_datasets(h5_patch):
                    # we only add to patchdatas is there are datasets
                    # in the hdf5 patch group.
                    # but we do create a Patch (below) anyway since
                    # we want empty patches to be there as well to access their attributes.
                    add_to_patchdata(patch_datas, h5_patch, basename, layout)

                patches.append(
                    Patch(
                        patch_datas,
                        h5_patch.name.split("/")[-1],
                        layout=layout,
                        attrs={k: v for k, v in h5_patch.attrs.items()},
                    )
                )

        patch_levels[ilvl] = PatchLevel(ilvl, patches)
    return patch_levels


def add_time_from_h5(hier, filepath, time, **kwargs):
    # add times to 'hier'
    # we may have a different selection box for that time as for already existing times
    # but we need to keep them, per time
    import h5py  # see doc/conventions.md section 2.1.1

    h5f = h5py.File(filepath, "r")
    selection_box = kwargs.get("selection_box", None)
    if hier.has_time(time):
        add_data_from_h5(hier, filepath, time)

    else:
        patch_levels = patch_levels_from_h5(h5f, time, selection_box=selection_box)
        hier.add_time(time, patch_levels, h5f, selection_box=selection_box)

    return hier


def add_data_from_h5(hier, filepath, time):
    """
    adds new PatchDatas to an existing PatchHierarchy for an existing time

    Data will be added from the given filepath.
    Data will be extracted from the selection box of the hierarchy at that time.
    """
    if not hier.has_time(time):
        raise ValueError("time does not exist in hierarchy")

    import h5py  # see doc/conventions.md section 2.1.1

    h5f = h5py.File(filepath, "r")

    # force using the hierarchy selection box at that time if existing
    if hier.selection_box is not None:
        selection_box = hier.selection_box[time]
    else:
        selection_box = None
    patch_levels = patch_levels_from_h5(h5f, time, selection_box=selection_box)

    for ilvl, lvl in hier.levels(time).items():
        for ip, patch in enumerate(lvl.patches):
            patch.patch_datas.update(patch_levels[ilvl].patches[ip].patch_datas)

    return hier


def new_from_h5(filepath, times, **kwargs):
    # create a patchhierarchy from a given time and optional selection box
    # loads all datasets from the filepath h5 file as patchdatas
    # we authorize user to pass only one selection box for all times
    # but in this case they're all the same
    import h5py  # see doc/conventions.md section 2.1.1

    selection_box = kwargs.get("selection_box", [None] * len(times))
    if none_iterable(selection_box) and all_iterables(times):
        selection_box = [selection_box] * len(times)

    patch_levels_per_time = []

    h5f = h5py.File(filepath, "r")
    for it, time in enumerate(times):
        if isinstance(time, float):
            time = f"{time:.10f}"
        patch_levels = patch_levels_from_h5(h5f, time, selection_box=selection_box[it])
        patch_levels_per_time.append(patch_levels)

    dim = len(h5f.attrs["domain_box"])
    domain_box = Box([0] * dim, h5f.attrs["domain_box"])

    # in hierarchy, we need to keep selection_box now
    # because we want that operations involving several hierarchies will need to check
    # that each time has the same patch layout.
    hier = PatchHierarchy(
        patch_levels_per_time,
        domain_box,
        refinement_ratio,
        times,
        h5f,
        selection_box=selection_box,
    )

    return hier


def hierarchy_fromh5(h5_filename, time=None, hier=None, silent=True, **kwargs):
    """
    creates a PatchHierarchy from a given time in a h5 file
    if hier is None, a new hierarchy is created
    if hier is not None, data is added to the hierarchy
    """
    if create_from_times(time, hier):
        time = listify(time)
        return new_from_h5(h5_filename, time, **kwargs)

    if create_from_all_times(time, hier):
        times = get_times_from_h5(h5_filename)
        return new_from_h5(h5_filename, times, **kwargs)

    if load_one_time(time, hier):
        time = listify(time)
        assert len(time) == 1
        return add_time_from_h5(hier, h5_filename, time[0], **kwargs)

    if load_all_times(time, hier):
        times = get_times_from_h5(h5_filename, as_float=False)
        for t in times:
            add_time_from_h5(hier, h5_filename, t, **kwargs)
        return hier

    assert False
