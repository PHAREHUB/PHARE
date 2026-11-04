import numpy as np
from pathlib import Path

from .patch import Patch
from ...core.box import Box
from .patchlevel import PatchLevel
from .hierarchy import PatchHierarchy
from .hierarchy import format_timestamp
from .hierarchy_utils import field_qties
from ...core.gridlayout import GridLayout
from ...core import phare_utilities as phut
from pyphare.core.phare_utilities import listify
from .for_vtk.patchdata import VtkFieldData as FieldData


base_path = "VTKHDF"
level_base_path = base_path + "/Level"
step_level_path = base_path + "/Steps/Level"
h5_time_ds_path = "/VTKHDF/Steps/Values"

_suffixes_per_ncols = {
    1: [""],
    3: ["_x", "_y", "_z"],
    6: ["_xx", "_xy", "_xz", "_yy", "_yz", "_zz"],
}


def get_path_from(h5file, string):
    keys = [v for v in string.split("/") if v]
    node = h5file
    for key in keys:
        node = node[key]
    return node


class VtkFile:
    def __init__(self, path):
        if not Path(path).exists():
            raise ValueError(f"ERROR: VTKHDF file does not exist: {path}")

        import h5py  # see doc/conventions.md section 2.1.1

        self.filepath = path
        self.file = h5py.File(path, "r")
        self._times = None
        self._step_times = None  # one per written step, duplicates included
        self._checked_levels = set()
        self._patches = {}  # {(ilvl, time_idx): (boxes, first data row of each patch)}
        self._domain_box = None
        self._level_spacing = {}

    def level_spacing(self, ilvl=0):
        if ilvl in self._level_spacing:
            return self._level_spacing[ilvl]
        level_group = self._get_path_from(level_base_path + str(ilvl))
        self._level_spacing[ilvl] = level_group.attrs["Spacing"][:]
        return self._level_spacing[ilvl]

    def max_num_levels(self):
        n = 0
        for i in range(11):
            if f"Level{i}" not in self.file["VTKHDF"]:
                break
            n += 1
        return n

    def time_idx(self, time):
        """index of the step holding time: a restart dumps its restart time again, the last
        copy is the one to read. Raises IndexError if time is not in the file"""
        self.times()
        steps = np.flatnonzero(np.isclose(self._step_times, float(time), 1e-10))
        if len(steps) == 0:
            raise IndexError(f"time {time} not in {self.filepath}")
        return steps[-1]

    def _data_set(self, ilvl):
        return self.file["VTKHDF"][f"Level{ilvl}"]["PointData"]["data"]

    @property
    def dimension(self):
        return self.file.attrs["dimension"]

    @property
    def interp_order(self):
        # only written for hybrid runs, as for phareh5 files 0 when absent (MHD)
        return (
            int(self.file.attrs["interpOrder"])
            if "interpOrder" in self.file.attrs
            else 0
        )

    @property
    def domain_box(self):
        if self._domain_box is None:
            dbox_attr = self.file.attrs["domain_box"]
            self._domain_box = Box([0] * len(dbox_attr), dbox_attr)
        return self._domain_box

    def has_time(self, time):
        try:
            return self.time_idx(time) is not None
        except IndexError:
            ...  # expected if time doesn't exist
        return False

    def times(self, reset=False):
        """the times in the file, once each. Only the first NSteps values of Steps/Values were
        written in full, a run that stopped during a dump can list one more"""
        if reset:
            self._times = None
        if self._times is None:
            values = self._get_path_from(h5_time_ds_path)[:]
            steps = self._get_path_from(base_path + "/Steps")
            if "NSteps" in steps.attrs:
                values = values[: int(steps.attrs["NSteps"])]
            self._step_times = values
            self._times = np.unique(values)
        return self._times

    def check_level_offsets(self, ilvl):
        """the per step offsets of a level must move forward and stay inside its datasets,
        else the file is corrupt, eg 32 bit offsets saturated on large files"""
        if ilvl in self._checked_levels:
            return
        self.times()
        nsteps = len(self._step_times)
        steps = self.file["VTKHDF"]["Steps"][f"Level{ilvl}"]
        num_boxes = np.asarray(steps["NumberOfAMRBox"][:nsteps], dtype=np.int64)
        box_offsets = np.asarray(steps["AMRBoxOffset"][:nsteps], dtype=np.int64)
        data_offsets = np.asarray(
            steps["PointDataOffset"]["data"][:nsteps], dtype=np.int64
        )
        num_box_rows = len(self.file["VTKHDF"][f"Level{ilvl}"]["AMRBox"])
        num_data_rows = len(self._data_set(ilvl))

        boxes_ok = np.all(box_offsets[1:] >= box_offsets[:-1] + num_boxes[:-1])
        boxes_ok = boxes_ok and np.all(box_offsets + num_boxes <= num_box_rows)
        data_ok = np.all(data_offsets[1:] >= data_offsets[:-1])
        data_ok = data_ok and np.all(
            data_offsets[1:][num_boxes[:-1] > 0] > data_offsets[:-1][num_boxes[:-1] > 0]
        )
        data_ok = data_ok and np.all(data_offsets <= num_data_rows)
        if not (boxes_ok and data_ok):
            raise ValueError(
                f"{self.filepath}: corrupt Steps/Level{ilvl} offsets, eg 32 bit offsets "
                "saturated in a file written before 64 bit offsets"
            )
        self._checked_levels.add(ilvl)

    def patches(self, ilvl, time_idx):
        """the boxes of a level at a step and the first data row of each patch, read once.
        Patch data is contiguous per level: all primal nodes, written 4 (1d) or 2 (2d) times
        to fill 3d cells"""
        key = (ilvl, time_idx)
        if key not in self._patches:
            self.check_level_offsets(ilvl)
            steps = self.file["VTKHDF"]["Steps"][f"Level{ilvl}"]
            num_boxes = int(steps["NumberOfAMRBox"][time_idx])
            box_offset = int(steps["AMRBoxOffset"][time_idx])
            data_offset = int(steps["PointDataOffset"]["data"][time_idx])
            rows = np.asarray(
                self.file["VTKHDF"][f"Level{ilvl}"]["AMRBox"][
                    box_offset : box_offset + num_boxes
                ]
            )
            dim = self.dimension
            lower, upper = rows[:, 0 : 2 * dim : 2], rows[:, 1 : 2 * dim : 2]
            sizes = np.prod(upper - lower + 2, axis=1) * [4, 2, 1][dim - 1]
            offsets = data_offset + np.concatenate(([0], np.cumsum(sizes)[:-1]))
            boxes = [Box(lo, up) for lo, up in zip(lower, upper)]
            self._patches[key] = (boxes, offsets.astype(np.int64))
        return self._patches[key]

    def _get_path_from(self, string):
        return get_path_from(self.file, string)


def get_all_available_quantities_from_h5(filepath, time=0, exclude=None, hier=None):
    if exclude is None:
        exclude = ["tags"]
    time = format_timestamp(time)
    path = Path(filepath)
    for h5 in path.glob("*.vtkhdf"):
        if h5.parent == path and h5.stem not in exclude:
            hier = hierarchy_fromvtkhdf(str(h5), time, hier)
    return hier


def make_layout(box, origin, cell_width, interp_order):
    ghosts_nbr = [0] * box.ndim
    return GridLayout(
        box, origin, cell_width, interp_order=interp_order, ghosts_nbr=ghosts_nbr
    )


def is_pop_fluid_file(basename):
    return (is_particle_file(basename) is False) and "pop" in basename


def is_particle_file(filename):
    return False


def pop_name(basename):
    return basename.split("_")[2]


def qty_from_basename(basename):
    if basename.startswith("ions_pop_"):
        return "_".join(basename.split("_")[3:])
    if basename.startswith("ions_"):
        return basename[len("ions_") :]
    return basename


def add_to_patchdata(vtk_file, ilvl, data_offset, patch_datas, basename, layout):
    """
    adds data in the h5_patch_grp in the given PatchData dict
    returns True if valid h5 patch found
    """

    if is_particle_file(basename):
        raise RuntimeError("Particle diagnostics are not supported under vtkhdf")

    dataset = vtk_file._data_set(ilvl)
    ncols = dataset.shape[1]  # decides if scalar/tensor

    if ncols not in _suffixes_per_ncols:
        raise RuntimeError(f"unsupported dataset column count {ncols} in {basename}")

    qty = qty_from_basename(basename)

    for cmp_id, suffix in enumerate(_suffixes_per_ncols[ncols]):
        ds_name = qty + suffix
        if ds_name not in field_qties:
            raise RuntimeError(
                "invalid dataset name : {} is not in {}".format(ds_name, field_qties)
            )
        pdata_name = field_qties[ds_name]
        pdata = FieldData(
            layout,
            pdata_name,
            dataset,
            cmp_id,
            data_offset,
        )
        if is_pop_fluid_file(basename):
            pdata_name = pop_name(basename) + "_" + pdata_name
        if pdata_name in patch_datas:
            raise ValueError("error - {} already in patchdata".format(ds_name))
        patch_datas[pdata_name] = pdata

    return True  # valid patch assumed


def patch_has_datasets(h5_patch_grp):
    return len(h5_patch_grp.keys()) > 0


def h5_filename_from(diagInfo):
    # diagInfo.quantity starts with a / , hence   [1:]
    return (diagInfo.quantity + ".vtkhdf").replace("/", "_")[1:]


def get_times_from_h5(filepath, as_float=True):
    times = VtkFile(filepath).times()  # written steps only, once each
    return times.astype(float) if as_float else list(times)


def create_from_all_times(time, hier):
    return time is None and hier is None


def create_from_times(times, hier):
    return times is not None and hier is None


def load_all_times(time, hier):
    return time is None and hier is not None


def load_one_time(time, hier):
    return time is not None and hier is not None


def patch_levels_from_h5(vtk_file, time, selection_box=None):
    """
    creates a dictionary of PatchLevels from a given time in a h5 file
    {ilvl: PatchLevel}
    """

    interp_order = vtk_file.interp_order
    basename = Path(vtk_file.file.filename).stem

    patch_levels = {}
    h5_patch = ""  # todo
    time_idx = vtk_file.time_idx(time)
    for ilvl in range(vtk_file.max_num_levels()):
        lvl_cell_width = vtk_file.level_spacing(ilvl)
        boxes, data_offsets = vtk_file.patches(ilvl, time_idx)

        patches = []

        for patch_box, data_offset in zip(boxes, data_offsets):
            origin = patch_box.lower * vtk_file.level_spacing(ilvl)

            intersect = None
            if selection_box is not None:
                pos_upper = [
                    orig + shape * dl
                    for orig, shape, dl in zip(origin, patch_box.shape, lvl_cell_width)
                ]
                pos_patch_box = Box(origin, pos_upper)
                intersect = selection_box * pos_patch_box

            if intersect is not None or selection_box is None:
                patch_datas = {}
                layout = make_layout(patch_box, origin, lvl_cell_width, interp_order)

                # currently, there is always data for patch
                add_to_patchdata(
                    vtk_file, ilvl, data_offset, patch_datas, basename, layout
                )
                patches.append(Patch(patch_datas, h5_patch, box=patch_box))

        if len(patches):
            patch_levels[ilvl] = PatchLevel(ilvl, patches)
    return patch_levels


def add_time_from_h5(hier, filepath, time, **kwargs):
    # add times to 'hier'
    # we may have a different selection box for that time as for already existing times
    # but we need to keep them, per time

    vtk_file = VtkFile(filepath)
    selection_box = kwargs.get("selection_box", None)
    if hier.has_time(time):
        add_data_from_h5(hier, filepath, time)

    else:
        patch_levels = patch_levels_from_h5(vtk_file, time, selection_box=selection_box)
        hier.add_time(time, patch_levels, vtk_file.file, selection_box=selection_box)

    return hier


def add_data_from_h5(hier, filepath, time):
    """
    adds new PatchDatas to an existing PatchHierarchy for an existing time

    Data will be added from the given filepath.
    Data will be extracted from the selection box of the hierarchy at that time.
    """
    if not hier.has_time(time):
        raise ValueError("time does not exist in hierarchy")

    vtk_file = VtkFile(filepath)

    # force using the hierarchy selection box at that time if existing
    selection_box = hier.selection_box.get(format_timestamp(time))
    patch_levels = patch_levels_from_h5(vtk_file, time, selection_box=selection_box)

    for ilvl, lvl in hier.levels(time).items():
        for ip, patch in enumerate(lvl.patches):
            patch.patch_datas.update(patch_levels[ilvl].patches[ip].patch_datas)

    return hier


def new_from_h5(filepath, times, **kwargs):
    selection_box = kwargs.get("selection_box", [None] * len(times))
    if phut.none_iterable(selection_box) and phut.all_iterables(times):
        selection_box = [selection_box] * len(times)

    time_hier = {}

    vtk_file = VtkFile(filepath)
    for it, time in enumerate(times):
        if isinstance(time, float):
            time = f"{time:.10f}"
        patch_levels = patch_levels_from_h5(
            vtk_file, time, selection_box=selection_box[it]
        )
        time_hier[time] = patch_levels

    # in hierarchy, we need to keep selection_box now
    # because we want that operations involving several hierarchies will need to check
    # that each time has the same patch layout.
    hier = PatchHierarchy(
        time_hier,
        vtk_file.domain_box,
        vtk_file.file,
        selection_box=selection_box,
    )

    return hier


def hierarchy_fromvtkhdf(h5_filename, time=None, hier=None, silent=True, **kwargs):
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
        for t in VtkFile(h5_filename).times():
            add_time_from_h5(hier, h5_filename, t, **kwargs)
        return hier

    raise ValueError("unreachable: all (time, hier) combinations handled above")
