import os
import numpy as np

from .patch import Patch
from .patchlevel import PatchLevel
from .patchdata import FieldData, ParticleData
from ..particles import Particles
from .hierarchy import PatchHierarchy
from ...core.box import Box
from ...core.phare_utilities import refinement_ratio
from ...core.gridlayout import GridLayout
from .hierarchy_utils import field_qties
import h5py


particle_files_patterns = ("domain", "patchGhost", "levelGhost")


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

            if dataset_name in patch_datas:
                raise ValueError("error - {} already in patchdata".format(dataset_name))

            patch_datas[pdata_name] = pdata

    return True  # valid patch assumed


def patch_has_datasets(h5_patch_grp):
    return len(h5_patch_grp.keys()) > 0


def h5_filename_from(diagInfo):
    # diagInfo.quantity starts with a / , hence   [1:]
    return (diagInfo.quantity + ".h5").replace("/", "_")[1:]


def get_times_from_h5(filepath):

    f = h5py.File(filepath, "r")
    times = np.array(sorted([float(s) for s in list(f["t"].keys())]))
    f.close()
    return times


h5_time_grp_key = "t"


def create_from_all_times(time, hier):
    return time is None and hier is None


def create_from_one_time(time, hier):
    return time is not None and hier is None


def load_all_times(time, hier):
    return time is None and hier is not None


def load_one_time(time, hier):
    return time is not None and hier is not None


def patch_levels_from_h5(filepath, time, selection_box=None):
    """
    creates a dictionary of PatchLevels from a given time in a h5 file
    {ilvl: PatchLevel}
    """
    import os

    h5f = h5py.File(filepath, "r")
    root_cell_width = h5f.attrs["cell_width"]
    interp_order = h5f.attrs["interpOrder"]
    basename = os.path.basename(filepath)

    patch_levels = {}

    for lvl_key, lvl in h5f["t"][time].items():

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
    return patch_levels, h5f


def add_time_from_h5(hier, filepath, time, selection_box=None):
    # add times to 'hier'
    # we may have a different selection box for that time as for already existing times
    # but we need to keep them, per time
    if hier.has_time(time):
        raise ValueError("time already exists in hierarchy")

    patch_levels, h5f = patch_levels_from_h5(
        filepath, time, selection_box=selection_box
    )

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

    # force using the hierarchy selection box at that time if existing
    patch_levels, h5f = patch_levels_from_h5(
        filepath, time, selection_box=hier.selection_box[time]
    )

    for ilvl, lvl in hier.levels(time).items():
        for ip, patch in enumerate(lvl.patches):
            patch.patch_datas.update(patch_levels[ilvl].patches[ip].patch_datas)

    return hier


def new_from_h5(filepath, time, selection_box=None):
    # create a patchhierarchy from a given time and optional selection box
    # loads all datasets from the filepath h5 file as patchdatas

    patch_levels, h5f = patch_levels_from_h5(
        filepath, time, selection_box=selection_box
    )
    dim = len(h5f.attrs["domain_box"])
    domain_box = Box([0] * dim, h5f.attrs["domain_box"])

    # in hierarchy, we need to keep selection_box now
    # because we want that operations involving several hierarchies will need to check
    # that each time has the same patch layout.
    hier = PatchHierarchy(
        patch_levels,
        domain_box,
        refinement_ratio,
        time,
        h5f,
        selection_box=selection_box,
    )

    return hier


def hierarchy_fromh5_(h5_filename, time=None, hier=None, silent=True):
    """
    creates a PatchHierarchy from a given time in a h5 file
    if hier is None, a new hierarchy is created
    if hier is not None, data is added to the hierarchy
    """
    if create_from_one_time(time, hier):
        return new_from_h5(h5_filename, time)

    if create_from_all_times(time, hier):
        times = get_times_from_h5(h5_filename)
        h = new_from_h5(h5_filename, times[0])
        for time in times[1:]:
            add_time_from_h5(h, h5_filename, time)

    if load_one_time(time, hier):
        return add_time_from_h5(hier, h5_filename, time)

    return add_data_from_h5(hier, h5_filename, time)


def hierarchy_fromh5(h5_filename, time, hier, silent=True):

    data_file = h5py.File(h5_filename, "r")
    basename = os.path.basename(h5_filename)

    root_cell_width = np.asarray(data_file.attrs["cell_width"])
    interp = data_file.attrs["interpOrder"]
    domain_box = Box(
        [0] * len(data_file.attrs["domain_box"]), data_file.attrs["domain_box"]
    )

    if create_from_all_times(time, hier):
        # first create from first time
        # then add all other times
        if not silent:
            print("creating hierarchy from all times in file")
        times = list(data_file[h5_time_grp_key].keys())
        hier = hierarchy_fromh5(h5_filename, time=times[0], hier=hier, silent=silent)
        if len(times) > 1:
            for t in times[1:]:
                hierarchy_fromh5(h5_filename, t, hier, silent=silent)
        return hier

    if create_from_one_time(time, hier):
        if not silent:
            print("creating hierarchy from time {}".format(time))
        t = time

        h5_time_grp = data_file[h5_time_grp_key][time]
        patch_levels = {}

        for plvl_key in h5_time_grp.keys():
            h5_patch_lvl_grp = h5_time_grp[plvl_key]
            ilvl = int(plvl_key[2:])
            lvl_cell_width = root_cell_width / refinement_ratio**ilvl
            patches = {}

            if ilvl not in patches:
                patches[ilvl] = []

            for pkey in h5_patch_lvl_grp.keys():
                h5_patch_grp = h5_patch_lvl_grp[pkey]

                patch_datas = {}
                layout = make_layout(h5_patch_grp, lvl_cell_width, interp)
                if patch_has_datasets(h5_patch_grp):
                    add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)

                patches[ilvl].append(
                    Patch(
                        patch_datas,
                        h5_patch_grp.name.split("/")[-1],
                        layout=layout,
                        attrs={k: v for k, v in h5_patch_grp.attrs.items()},
                    )
                )

            patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

        diag_hier = PatchHierarchy(
            patch_levels, domain_box, refinement_ratio, t, data_files=data_file
        )

        return diag_hier

    if load_one_time(time, hier):
        if not silent:
            print("loading data at time {} into existing hierarchy".format(time))
        h5_time_grp = data_file[h5_time_grp_key][time]
        t = time

        if t in hier.time_hier:
            if not silent:
                print("time already exist, adding data...")

            # time already exists in the hierarchy
            # all we need to do is adding the data
            # as patchDatas in the appropriate patches
            # and levels, if data compatible with hierarchy

            patch_levels = hier.time_hier[t]

            for plvl_key in h5_time_grp.keys():
                ilvl = int(plvl_key[2:])
                lvl_cell_width = root_cell_width / refinement_ratio**ilvl

                for ipatch, pkey in enumerate(h5_time_grp[plvl_key].keys()):
                    h5_patch_grp = h5_time_grp[plvl_key][pkey]

                    if patch_has_datasets(h5_patch_grp):
                        hier_patch = patch_levels[ilvl].patches[ipatch]
                        origin = h5_time_grp[plvl_key][pkey].attrs["origin"]
                        upper = h5_time_grp[plvl_key][pkey].attrs["upper"]
                        lower = h5_time_grp[plvl_key][pkey].attrs["lower"]
                        file_patch_box = Box(lower, upper)

                        assert file_patch_box == hier_patch.box
                        assert (abs(origin - hier_patch.origin) < 1e-6).all()
                        assert (abs(lvl_cell_width - hier_patch.dl) < 1e-6).all()

                        layout = make_layout(h5_patch_grp, lvl_cell_width, interp)
                        add_to_patchdata(
                            hier_patch.patch_datas, h5_patch_grp, basename, layout
                        )

            return hier

        if not silent:
            print("adding data to new time")
        # time does not exist in the hierarchy
        # we have to create a brand new set of patchLevels
        # containing patches, and load data in their patchdatas

        patch_levels = {}

        for plvl_key in h5_time_grp.keys():
            ilvl = int(plvl_key[2:])

            lvl_cell_width = root_cell_width / refinement_ratio**ilvl
            lvl_patches = []

            for ipatch, pkey in enumerate(h5_time_grp[plvl_key].keys()):
                h5_patch_grp = h5_time_grp[plvl_key][pkey]

                layout = make_layout(h5_patch_grp, lvl_cell_width, interp)
                patch_datas = {}
                if patch_has_datasets(h5_patch_grp):
                    add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)
                lvl_patches.append(
                    Patch(
                        patch_datas,
                        h5_patch_grp.name.split("/")[-1],
                        layout=layout,
                        attrs={k: v for k, v in h5_patch_grp.attrs.items()},
                    )
                )

            patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

        hier.time_hier[t] = patch_levels
        return hier

    if load_all_times(time, hier):
        if not silent:
            print("loading all times in existing hier")
        for time in data_file[h5_time_grp_key].keys():
            hier = hierarchy_fromh5(h5_filename, time, hier, silent=silent)

        return hier
