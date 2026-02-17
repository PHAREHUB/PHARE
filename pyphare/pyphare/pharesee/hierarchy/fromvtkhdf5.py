import numpy as np
from pathlib import Path

from .for_vtk.patch import Patch
from .for_vtk.patchlevel import PatchLevel
from .for_vtk.patchdata import VtkFieldData as FieldData
from .for_vtk.hierarchy import PatchHierarchy
from .hierarchy import format_timestamp
from ...core.box import Box
from ...core.phare_utilities import (
    refinement_ratio,
    none_iterable,
    all_iterables,
)
from ...core.gridlayout import GridLayout
from .hierarchy_utils import field_qties

from pyphare.core import phare_utilities as phut
from pyphare.core.phare_utilities import listify

base_path = "VTKHDF"
level_base_path = base_path + "/Level"
step_level_path = base_path + "/Steps/Level"
h5_time_ds_path = "/VTKHDF/Steps/Values"

_qty_per_filename = {
    "EM_E": "EM_E",
    "EM_B": "EM_B",
    "ions_bulkVelocity": "bulkVelocity",
}

_vec_fields = {
    "EM_E",
    "EM_B",
    "bulkVelocity",
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
        self._domain_box = None
        self._level_spacing = {}

    def level_spacing(self, ilvl=0):
        if ilvl in self._level_spacing:
            return self._level_spacing[ilvl]
        level_group = self._get_path_from(level_base_path + str(ilvl))
        self._level_spacing[ilvl] = level_group.attrs["Spacing"][:]
        return self._level_spacing[ilvl]

    def level_boxes(self, ilvl=0, time=0):
        dim = self.dimension
        time_idx = self.time_idx(time)
        num_boxes = self._get_path_from(
            step_level_path + str(ilvl) + "/NumberOfAMRBox"
        )[time_idx]
        box_offset = self._get_path_from(step_level_path + str(ilvl) + "/AMRBoxOffset")[
            time_idx
        ]
        amr_box_ds = self._get_path_from(level_base_path + str(ilvl) + "/AMRBox")
        boxes = [0] * num_boxes

        for bi, boff in enumerate(range(box_offset, box_offset + num_boxes)):
            box_vals = amr_box_ds[boff]
            lo = np.zeros(dim)
            up = np.zeros(dim)
            for i, di in enumerate(range(0, dim * 2, 2)):
                lo[i] = box_vals[di + 0]
                up[i] = box_vals[di + 1]
            boxes[bi] = Box(lo, up)
        return boxes

    def max_num_levels(self):
        n = 0
        for i in range(11):
            if f"Level{i}" not in self.file["VTKHDF"]:
                break
            n += 1
        return n

    def time_idx(self, time):
        times = self.times()
        ret = np.where(np.isclose(times, float(time), 1e-10))
        return ret[0][0]

    def _amr_box_offset(self, ilvl, time_idx):
        return self.file["VTKHDF"]["Steps"][f"Level{ilvl}"]["AMRBoxOffset"][time_idx]

    def _num_boxes(self, ilvl, time_idx):
        return self.file["VTKHDF"]["Steps"][f"Level{ilvl}"]["NumberOfAMRBox"][time_idx]

    def _data_offset(self, ilvl, time_idx):
        return self.file["VTKHDF"]["Steps"][f"Level{ilvl}"]["PointDataOffset"]["data"][
            time_idx
        ]

    def _data_set(self, ilvl):
        return self.file["VTKHDF"][f"Level{ilvl}"]["PointData"]["data"]

    @property
    def dimension(self):
        return self.file.attrs["dimension"]

    @property
    def interp_order(self):
        return self.file.attrs["interpOrder"]

    @property
    def domain_box(self):
        if self._domain_box is None:
            dbox_attr = self.file.attrs["domain_box"]
            self._domain_box = Box([0] * len(dbox_attr), dbox_attr)
        return self._domain_box

    def has_time(self, time):
        return self.time_idx(time) is not None

    def times(self, reset=False):
        if reset:
            self._times = None
        if self._times is None:
            self._times = self._get_path_from(h5_time_ds_path)[:]
        return self._times

    def _get_path_from(self, string):
        return get_path_from(self.file, string)


class VtkPatchLevelInfo:
    def __init__(self, vtk_file, ilvl, time_idx):
        self.vtk_file = vtk_file
        self.ilvl = ilvl
        self.time_idx = time_idx
        self.num_boxes = vtk_file._num_boxes(ilvl, self.time_idx)
        self.data_offset = vtk_file._data_offset(ilvl, self.time_idx)
        self.box_offset = vtk_file._amr_box_offset(ilvl, self.time_idx)
        self.rolling_data_offset = self.data_offset

    def boxes(self):
        dim = self.vtk_file.dimension
        amr_box_ds = self.vtk_file.file["VTKHDF"][f"Level{self.ilvl}"]["AMRBox"]
        box_vals_list = amr_box_ds[self.box_offset : self.box_offset + self.num_boxes]
        boxes = [0] * self.num_boxes
        for bi in range(self.num_boxes):
            lo = np.zeros(dim)
            up = np.zeros(dim)
            box_vals = box_vals_list[bi]
            for i, di in enumerate(range(0, dim * 2, 2)):
                lo[i] = box_vals[di + 0]
                up[i] = box_vals[di + 1]
            boxes[bi] = Box(lo, up)
        return boxes

    def __deepcopy__(self, memo):
        no_copy_keys = ["vtk_file"]  # do not copy these things
        cpy = phut.deep_copy(self, memo, no_copy_keys)
        cpy.vtk_file = self.vtk_file
        return cpy


def get_all_available_quantities_from_h5(filepath, time=0, exclude=["tags"], hier=None):
    time = format_timestamp(time)
    path = Path(filepath)
    for h5 in path.glob("*.vtkhdf"):
        if h5.parent == path and h5.stem not in exclude:
            hier = hierarchy_fromvtkhdf(str(h5), time, hier)
    return hier


def make_layout(box, origin, cell_width, interp_order):
    return GridLayout(box, origin, cell_width, interp_order=interp_order)


def is_pop_fluid_file(basename):
    return (is_particle_file(basename) is False) and "pop" in basename


def is_particle_file(filename):
    return False


def pop_name(basename):
    return Path(".vtkhdf").stem.split("_")[2]


def add_to_patchdata(vtk_file, lvl_info, patch_idx, patch_datas, basename, layout):
    """
    adds data in the h5_patch_grp in the given PatchData dict
    returns True if valid h5 patch found
    """

    X_TIMES = [4, 2, 1][vtk_file.dimension - 1]  # data in file is sometimes duplicated

    if is_particle_file(basename):
        raise RuntimeError("Particle diagnostics are not supported under vtkhdf")

    def _do_field(cmp_id, qty):
        if qty not in field_qties:
            raise RuntimeError(
                "invalid dataset name : {} is not in {}".format(qty, field_qties)
            )
        pdata_name = field_qties[qty]
        pdata = FieldData(
            layout,
            pdata_name,
            vtk_file._data_set(lvl_info.ilvl),
            cmp_id,
            lvl_info.rolling_data_offset,
        )

        if is_pop_fluid_file(basename):
            pdata_name = pop_name(basename) + "_" + pdata_name
        if pdata_name in patch_datas:
            raise ValueError("error - {} already in patchdata".format(qty))
        patch_datas[pdata_name] = pdata
        return pdata.field_box.size()  # constant as all primal

    data_size = 0
    qty = _qty_per_filename[basename]
    if qty in _vec_fields:
        for cmp_id, cmp in enumerate(["x", "y", "z"]):
            data_size = _do_field(cmp_id, f"{qty}_{cmp}")
    else:  # assume scalar field
        data_size = _do_field(0, qty)

    lvl_info.rolling_data_offset += data_size * X_TIMES
    return True  # valid patch assumed


def patch_has_datasets(h5_patch_grp):
    return len(h5_patch_grp.keys()) > 0


def h5_filename_from(diagInfo):
    # diagInfo.quantity starts with a / , hence   [1:]
    return (diagInfo.quantity + ".vtkhdf").replace("/", "_")[1:]


def get_times_from_h5(filepath, as_float=True):
    import h5py  # see doc/conventions.md section 2.1.1

    with h5py.File(filepath, "r") as f:
        ds = get_path_from(f, h5_time_ds_path)
        times = list(ds[:])
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
        vtk_patch_level_info = VtkPatchLevelInfo(vtk_file, ilvl, time_idx)

        patches = []

        for ipatch, patch_box in enumerate(vtk_patch_level_info.boxes()):
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
                    vtk_file,
                    vtk_patch_level_info,
                    ipatch,
                    patch_datas,
                    basename,
                    layout,
                )
                patches.append(Patch(patch_datas, h5_patch, layout=layout))

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
    if hier.selection_box is not None:
        selection_box = hier.selection_box[time]
    else:
        selection_box = None
    patch_levels = patch_levels_from_h5(vtk_file, time, selection_box=selection_box)

    for ilvl, lvl in hier.levels(time).items():
        for ip, patch in enumerate(lvl.patches):
            patch.patch_datas.update(patch_levels[ilvl].patches[ip].patch_datas)

    return hier


def new_from_h5(filepath, times, **kwargs):
    selection_box = kwargs.get("selection_box", [None] * len(times))
    if none_iterable(selection_box) and all_iterables(times):
        selection_box = [selection_box] * len(times)

    patch_levels_per_time = []

    vtk_file = VtkFile(filepath)
    for it, time in enumerate(times):
        if isinstance(time, float):
            time = f"{time:.10f}"
        patch_levels = patch_levels_from_h5(
            vtk_file, time, selection_box=selection_box[it]
        )
        patch_levels_per_time.append(patch_levels)

    # in hierarchy, we need to keep selection_box now
    # because we want that operations involving several hierarchies will need to check
    # that each time has the same patch layout.
    hier = PatchHierarchy(
        patch_levels_per_time,
        vtk_file.domain_box,
        refinement_ratio,
        times,
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
        for t in get_times_from_h5(h5_filename):
            add_time_from_h5(hier, h5_filename, t, **kwargs)
        return hier

    assert False
