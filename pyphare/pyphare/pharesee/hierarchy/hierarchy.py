#
#
#


import numpy as np

from .patch import Patch
from ...core import box as boxm
from .patchlevel import PatchLevel
from ...core import phare_utilities as phut


def format_timestamp(timestamp):
    if isinstance(timestamp, str):
        return timestamp
    return "{:.10f}".format(timestamp)


class PatchHierarchy(object):
    """is a collection of patch levels"""

    def __init__(
        self,
        patch_levels,
        domain_box,
        refinement_ratio=2,
        times=[0.0],
        data_files=None,
        selection_box=None,
        ephemerals=None,
    ):
        if not isinstance(times, (tuple, list)):
            times = phut.listify(times)

        if not isinstance(patch_levels, (tuple, list)):
            patch_levels = phut.listify(patch_levels)

        self.selection_box = selection_box
        if self.selection_box is not None:
            if not isinstance(self.selection_box, (tuple, list)):
                self.selection_box = phut.listify(self.selection_box)
            self.selection_box = {
                format_timestamp(t): box for t, box in zip(times, self.selection_box)
            }
            assert len(times) == len(self.selection_box)

        assert len(times) == len(patch_levels)

        self.patch_levels = patch_levels
        self.ndim = len(domain_box.lower)
        self.time_hier = {}
        self.time_hier.update(
            {format_timestamp(t): pl for t, pl in zip(times, patch_levels)}
        )

        self.domain_box = domain_box
        self.refinement_ratio = refinement_ratio

        self._sim = None

        self.data_files = {}
        if data_files is not None and isinstance(data_files, dict):
            self.data_files = data_files
        elif data_files is not None:
            if hasattr(self, "data_files"):
                self.data_files.update({data_files.filename: data_files})
            else:
                self.data_files = {data_files.filename: data_files}

        self.ephemerals = ephemerals
        self.update()

    def __deepcopy__(self, memo):
        no_copy_keys = ["data_files"]  # do not copy these things
        return phut.deep_copy(self, memo, no_copy_keys)

    def __getitem__(self, qty):
        return self.__dict__[qty]

    def update(self):
        if len(self.quantities()) > 1:
            for qty in self.quantities():
                for time, levels in self.time_hier.items():
                    new_lvls = {}
                    for ilvl, level in levels.items():
                        patches = []
                        for patch in level.patches:
                            patches += [Patch({qty: patch.patch_datas[qty]}, patch.id)]
                        new_lvls[ilvl] = PatchLevel(ilvl, patches)
                    if qty not in self.__dict__:
                        self.__dict__[qty] = PatchHierarchy(
                            new_lvls,
                            self.domain_box,
                            selection_box=self.domain_box,
                            times=time,
                            data_files=self.data_files,
                        )
                    else:
                        self.__dict__[qty].time_hier[time] = new_lvls

    def nbytes(self):
        n = 0
        for t in self.times():
            for lvl in self.levels(t).values():
                for p in lvl.patches:
                    for pd in p.patch_datas.values():
                        n += pd.dataset.nbytes
        return n

    def nbrPatches(self):
        n = 0
        for t in self.times():
            for lvl in self.levels(t).values():
                n += len(lvl.patches)
        return n

    @property
    def sim(self):
        if self._sim:
            return self._sim

        # data_files has a key/value per h5 filename.
        # but the "serialized_simulation" in "py_attrs" should be the same for all files
        # used by the hierarchy. So we just take the first one.
        first_file = list(self.data_files.values())[0]
        if "py_attrs" not in first_file.keys():
            raise ValueError("Simulation is not available for deserialization")

        from ...pharein.simulation import deserialize

        try:
            self._sim = deserialize(
                first_file["py_attrs"].attrs["serialized_simulation"]
            )
        except Exception as e:
            raise RuntimeError(f"Failed to deserialize simulation from data file : {e}")
        return self._sim

    def __call__(self, qty=None, **kwargs):
        # take slice/slab of 1/2d array from 2/3d array
        def cuts(c, coord):
            return c > coord.min() and c < coord.max()

        class Extractor:
            def __init__(self):
                self.exclusions = []

            def extract(self, coord, data):
                mask = coord == coord
                for exclusion in self.exclusions:
                    idx = np.where(
                        (coord > exclusion[0] - 1e-6) & (coord < exclusion[1] + 1e-6)
                    )[0]
                    mask[idx] = False

                self.exclusions += [(coord.min(), coord.max())]
                return coord[mask], data[mask]

        def domain_coords(patch, qty):
            pd = patch.patch_datas[qty]
            nbrGhosts = pd.ghosts_nbr[0]
            return pd.x[nbrGhosts:-nbrGhosts], pd.y[nbrGhosts:-nbrGhosts]

        if len(kwargs) < 1 or len(kwargs) > 3:
            raise ValueError("Error - must provide coordinates")
        if qty is None:
            if len(self.quantities()) == 1:
                qty = self.quantities()[0]
            else:
                raise ValueError(
                    "The PatchHierarchy has several quantities but none is specified"
                )

        if "x" in kwargs:
            c = kwargs["x"]
            slice_dim = 1
            cst_dim = 0
        else:
            c = kwargs["y"]
            slice_dim = 0
            cst_dim = 1

        extractor = Extractor()
        datas = []
        coords = []
        ilvls = list(self.levels().keys())[::-1]

        for ilvl in ilvls:
            lvl = self.patch_levels[ilvl]
            for patch in lvl.patches:
                slice_coord = domain_coords(patch, qty)[slice_dim]
                cst_coord = domain_coords(patch, qty)[cst_dim]

                if cuts(c, cst_coord):
                    data = patch(qty, **kwargs)
                    coord_keep, data_keep = extractor.extract(slice_coord, data)
                    datas += [data_keep]
                    coords += [coord_keep]

        cut = np.concatenate(datas)
        coords = np.concatenate(coords)
        ic = np.argsort(coords)
        coords = coords[ic]
        cut = cut[ic]
        return coords, cut

    def _default_time(self):
        return self.times()[0]

    def finest_level(self, time=None):
        if time is None:
            time = self._default_time()
        return max(list(self.levels(time=time).keys()))

    def levels(self, time=None):
        if time is None:
            time = self._default_time()
        return self.time_hier[format_timestamp(time)]

    def level(self, level_number, time=None):
        return self.levels(time)[level_number]

    def levelNbr(self, time=None):
        if time is None:
            time = self._default_time()
        return len(self.levels(time).items())

    def levelNbrs(self, time=None):
        if time is None:
            time = self._default_time()
        return list(self.levels(time).keys())

    def add_time(self, time, patch_level, h5file, selection_box=None):
        formated_time = format_timestamp(time)

        self.time_hier[format_timestamp(time)] = patch_level
        if selection_box is not None:
            self.selection_box[formated_time] = selection_box

        self.data_files[h5file.filename] = h5file
        self.update()

    def is_homogeneous(self):
        """
        return True if all patches of all levels at all times
        have the same patch data quantities
        """
        qties = self._quantities()
        it_is = True
        for time, levels in self.time_hier.items():
            for ilvl, lvl in levels.items():
                for patch in lvl.patches:
                    pdnames = list(patch.patch_datas.keys())
                    if len(pdnames):  # do not compare empty patches
                        it_is &= qties == pdnames
        return it_is

    def _quantities(self):
        # we return the name of the patchdatas of the first level that has
        # patches with data. checking that patchdatas are not {} is important
        # since some patches might be empty (e.g. level ghost patchdatas on L0)
        for ilvl, lvl in self.levels().items():
            if len(lvl.patches) > 0 and len(lvl.patches[0].patch_datas):
                return list(self.level(ilvl).patches[0].patch_datas.keys())
        return []

    def quantities(self):
        if not self.is_homogeneous():
            print("WARNING - hierarchy is not homogeneous")
        return self._quantities()

    def global_min(self, qty, **kwargs):
        time = kwargs.get("time", self._default_time())
        first = True
        for ilvl, lvl in self.levels(time).items():
            for patch in lvl.patches:
                pd = patch.patch_datas[qty]
                if first:
                    m = np.nanmin(pd.dataset[:])
                    first = False
                else:
                    data_and_min = np.concatenate(([m], pd.dataset[:].flatten()))
                    m = np.nanmin(data_and_min)

        return m

    def global_max(self, qty, **kwargs):
        time = kwargs.get("time", self._default_time())
        first = True
        for _, lvl in self.levels(time).items():
            for patch in lvl.patches:
                pd = patch.patch_datas[qty]
                if first:
                    m = np.nanmax(pd.dataset[:])
                    first = False
                else:
                    data_and_max = np.concatenate(([m], pd.dataset[:].flatten()))
                    m = np.nanmax(data_and_max)

        return m

    def refined_domain_box(self, level_number):
        """
        returns the domain box refined for a given level number
        """
        assert level_number >= 0
        return boxm.refine(self.domain_box, self.refinement_ratio**level_number)

    def level_domain_box(self, level_number):
        if level_number == 0:
            return self.domain_box
        return self.refined_domain_box(level_number)

    def __str__(self):
        s = "Hierarchy: \n"
        for t, patch_levels in self.time_hier.items():
            s = s + "Time {}\n".format(t)
            for ilvl, lvl in patch_levels.items():
                s = s + "Level {}\n".format(ilvl)
                for ip, patch in enumerate(lvl.patches):
                    for qty_name, pd in patch.patch_datas.items():
                        pdstr = "    P{ip} {type} {pdname} box is {box} and ghost box is {gbox}"
                        s = s + pdstr.format(
                            ip=ip,
                            type=type(pd.dataset),
                            pdname=qty_name,
                            box=patch.box,
                            gbox=pd.ghost_box,
                        )
                        s = s + "\n"
        return s

    def has_time(self, time):
        return format_timestamp(time) in self.time_hier

    def has_file(self, filename):
        return filename in self.data_files

    def times(self):
        # return np.sort(np.asarray(list(self.time_hier.keys()), dtype=np.float32))
        return np.sort(np.asarray(list(self.time_hier.keys())))

    def plot_patches(self, save=False):
        from .plotting.plot_fields import plot_patches

        return plot_patches(self, save=save)

    def box_to_Rectangle(self, box):
        from .plotting import box_to_Rectangle

        return box_to_Rectangle(box)

    def plot_2d_patches(self, ilvl, collections, **kwargs):
        from .plotting.plot_fields import plot_2d_patches

        return plot_2d_patches(self, ilvl, collections, **kwargs)

    def plot1d(self, **kwargs):
        from .plotting.plot_fields import plot1d

        return plot1d(self, **kwargs)

    def plot2d(self, **kwargs):
        from .plotting.plot_fields import plot2d

        return plot2d(self, **kwargs)

    def plot(self, **kwargs):
        from .plotting.plot_fields import plot

        return plot(self, **kwargs)

    def dist_plot(self, **kwargs):
        from .plotting.plot_particles import hierarchy_dist_plot

        return hierarchy_dist_plot(self, **kwargs)
