import numpy as np
import matplotlib.pyplot as plt

from .patch import Patch
from .patchlevel import PatchLevel
from ...core.box import Box
from ...core import box as boxm
from ...core.phare_utilities import listify
from ...core.phare_utilities import deep_copy
from ...core.phare_utilities import refinement_ratio


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
        **kwargs,
    ):
        if not isinstance(times, (tuple, list)):
            times = listify(times)

        if not isinstance(patch_levels, (tuple, list)):
            patch_levels = listify(patch_levels)

        self.selection_box = kwargs.get("selection_box", None)
        if self.selection_box is not None:
            if not isinstance(self.selection_box, (tuple, list)):
                self.selection_box = listify(self.selection_box)
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

        if data_files is not None and isinstance(data_files, dict):
            self.data_files = data_files
        elif data_files is not None:
            if hasattr(self, "data_files"):
                self.data_files.update({data_files.filename: data_files})
            else:
                self.data_files = {data_files.filename: data_files}
        else:
            self.data_files = {}

        self.update()

    def __deepcopy__(self, memo):
        no_copy_keys = ["data_files"]  # do not copy these things
        return deep_copy(self, memo, no_copy_keys)

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
        fig, ax = plt.subplots(figsize=(10, 3))
        for ilvl, lvl in self.levels(0.0).items():
            lvl_offset = ilvl * 0.1
            for patch in lvl.patches:
                dx = patch.dl[0]
                x0 = patch.box.lower * dx
                x1 = patch.box.upper * dx
                xcells = np.arange(x0, x1 + dx, dx)
                y = lvl_offset + np.zeros_like(xcells)
                ax.plot(xcells, y, marker=".")

        if save:
            fig.savefig("hierarchy.png")

    def box_to_Rectangle(self, box):
        from matplotlib.patches import Rectangle

        return Rectangle(box.lower, *box.shape)

    def plot_2d_patches(self, ilvl, collections, **kwargs):
        if isinstance(collections, list) and all(
            [isinstance(el, Box) for el in collections]
        ):
            collections = [{"boxes": collections}]

        from matplotlib.collections import PatchCollection

        level_domain_box = self.level_domain_box(ilvl)
        mi, ma = level_domain_box.lower.min(), level_domain_box.upper.max()

        fig, ax = kwargs.get("subplot", plt.subplots(figsize=(6, 6)))

        for collection in collections:
            facecolor = collection.get("facecolor", "none")
            edgecolor = collection.get("edgecolor", "purple")
            alpha = collection.get("alpha", 1)
            rects = [self.box_to_Rectangle(box) for box in collection["boxes"]]

            ax.add_collection(
                PatchCollection(
                    rects, facecolor=facecolor, alpha=alpha, edgecolor=edgecolor
                )
            )

        if "title" in kwargs:
            from textwrap import wrap

            xfigsize = int(fig.get_size_inches()[0] * 10)  # 10 characters per inch
            ax.set_title("\n".join(wrap(kwargs["title"], xfigsize)))

        major_ticks = np.arange(mi - 5, ma + 5 + 5, 5)
        ax.set_xticks(major_ticks)
        ax.set_yticks(major_ticks)

        minor_ticks = np.arange(mi - 5, ma + 5 + 5, 1)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(minor_ticks, minor=True)

        ax.grid(which="both")

        return fig

    def plot1d(self, **kwargs):
        """
        plot
        """
        usr_lvls = kwargs.get("levels", (0,))
        qty = kwargs.get("qty", None)
        time = kwargs.get("time", self.times()[0])

        if "ax" not in kwargs:
            fig, ax = plt.subplots()
        else:
            ax = kwargs["ax"]
            fig = ax.figure
        for lvl_nbr, level in self.levels(time).items():
            if lvl_nbr not in usr_lvls:
                continue
            for ip, patch in enumerate(level.patches):
                pdata_nbr = len(patch.patch_datas)
                pdata_names = list(patch.patch_datas.keys())
                if qty is None and pdata_nbr != 1:
                    multiple = "multiple quantities in patch, "
                    err = (
                        multiple
                        + "please specify a quantity in "
                        + " ".join(pdata_names)
                    )
                    raise ValueError(err)
                if qty is None:
                    qty = pdata_names[0]

                layout = patch.patch_datas[qty].layout
                nbrGhosts = layout.nbrGhostFor(qty)
                val = patch.patch_datas[qty][patch.box]
                x = patch.patch_datas[qty].x[nbrGhosts[0] : -nbrGhosts[0]]
                label = "L{level}P{patch}".format(level=lvl_nbr, patch=ip)
                marker = kwargs.get("marker", "")
                ls = kwargs.get("ls", "--")
                color = kwargs.get("color", "k")
                ax.plot(x, val, label=label, marker=marker, ls=ls, color=color)

        ax.set_title(kwargs.get("title", ""))
        ax.set_xlabel(kwargs.get("xlabel", "x"))
        ax.set_ylabel(kwargs.get("ylabel", qty))
        if "xlim" in kwargs:
            ax.set_xlim(kwargs["xlim"])
        if "ylim" in kwargs:
            ax.set_ylim(kwargs["ylim"])

        if kwargs.get("legend", None) is not None:
            ax.legend()

        if "filename" in kwargs:
            fig.savefig(kwargs["filename"])

    def plot2d(self, **kwargs):
        from matplotlib.patches import Rectangle
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        time = kwargs.get("time", self._default_time())
        usr_lvls = kwargs.get("levels", self.levelNbrs(time))
        default_qty = None
        if len(self.quantities()) == 1:
            default_qty = self.quantities()[0]
        qty = kwargs.get("qty", default_qty)

        if "ax" not in kwargs:
            fig, ax = plt.subplots()
        else:
            ax = kwargs["ax"]
            fig = ax.figure

        glob_min = self.global_min(qty)
        glob_max = self.global_max(qty)
        # assumes max 5 levels...
        patchcolors = {ilvl: "k" for ilvl in usr_lvls}
        patchcolors = kwargs.get("patchcolors", patchcolors)
        if not isinstance(patchcolors, dict):
            patchcolors = dict(zip(usr_lvls, patchcolors))

        linewidths = {ilvl: 1 for ilvl in usr_lvls}
        linewidths = kwargs.get("lw", linewidths)
        if not isinstance(linewidths, dict):
            linewidths = dict(zip(usr_lvls, linewidths))
        linestyles = {ilvl: "-" for ilvl in usr_lvls}
        linestyles = kwargs.get("ls", linestyles)
        if not isinstance(linestyles, dict):
            linestyles = dict(zip(usr_lvls, linestyles))

        for lvl_nbr, lvl in self.levels(time).items():
            if lvl_nbr not in usr_lvls:
                continue
            for patch in self.level(lvl_nbr, time).patches:
                pdat = patch.patch_datas[qty]
                data = pdat.dataset[:]
                nbrGhosts = pdat.ghosts_nbr
                x = pdat.x
                y = pdat.y

                # if nbrGhosts is 0, we cannot do array[0,-0]
                if np.all(nbrGhosts == np.zeros_like(nbrGhosts)):
                    x = np.copy(x)
                    y = np.copy(y)
                else:
                    data = pdat[patch.box]
                    x = np.copy(x[nbrGhosts[0] : -nbrGhosts[0]])
                    y = np.copy(y[nbrGhosts[1] : -nbrGhosts[1]])
                dx, dy = pdat.layout.dl
                x -= dx * 0.5
                y -= dy * 0.5
                x = np.append(x, x[-1] + dx)
                y = np.append(y, y[-1] + dy)
                im = ax.pcolormesh(
                    x,
                    y,
                    data.T,
                    cmap=kwargs.get("cmap", "Spectral_r"),
                    vmin=kwargs.get("vmin", glob_min - 1e-6),
                    vmax=kwargs.get("vmax", glob_max + 1e-6),
                )

                if kwargs.get("plot_patches", False) is True:
                    r = Rectangle(
                        (patch.box.lower[0] * dx, patch.box.lower[1] * dy),
                        patch.box.shape[0] * dx,
                        patch.box.shape[1] * dy,
                        fc="none",
                        ec=patchcolors[lvl_nbr],
                        alpha=0.4,
                        lw=linewidths[lvl_nbr],
                        ls=linestyles[lvl_nbr],
                    )
                    ax.add_patch(r)

        ax.set_aspect(kwargs.get("aspect", "equal"))
        ax.set_title(kwargs.get("title", ""))
        ax.set_xlabel(kwargs.get("xlabel", "x"))
        ax.set_ylabel(kwargs.get("ylabel", "y"))
        if "xlim" in kwargs:
            ax.set_xlim(kwargs["xlim"])
        if "ylim" in kwargs:
            ax.set_ylim(kwargs["ylim"])

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        plt.colorbar(im, ax=ax, cax=cax)

        if kwargs.get("legend", None) is not None:
            ax.legend()

        if "filename" in kwargs:
            fig.savefig(kwargs["filename"], dpi=kwargs.get("dpi", 200))

        return fig, ax

    def plot(self, **kwargs):
        if self.ndim == 1:
            return self.plot1d(**kwargs)
        elif self.ndim == 2:
            return self.plot2d(**kwargs)

    def dist_plot(self, **kwargs):
        """
        plot phase space of a particle hierarchy
        """
        import copy

        from ..plotting import dist_plot as dp

        usr_lvls = kwargs.get("levels", (0,))
        finest = kwargs.get("finest", False)
        pops = kwargs.get("pop", [])
        time = kwargs.get("time", self.times()[0])
        axis = kwargs.get("axis", ("Vx", "Vy"))
        all_pops = list(self.level(0, time).patches[0].patch_datas.keys())

        vmin = kwargs.get("vmin", -2)
        vmax = kwargs.get("vmax", 2)
        dv = kwargs.get("dv", 0.05)
        vbins = vmin + dv * np.arange(int((vmax - vmin) / dv))

        if finest:
            final = finest_part_data(self)
            if axis[0] == "x":
                xbins = amr_grid(self, time)
                bins = (xbins, vbins)
            else:
                bins = (vbins, vbins)
            kwargs["bins"] = bins

        else:
            final = {pop: None for pop in all_pops}
            for lvl_nbr, level in self.levels(time).items():
                if lvl_nbr not in usr_lvls:
                    continue
                for ip, patch in enumerate(level.patches):
                    if len(pops) == 0:
                        pops = list(patch.patch_datas.keys())

                    for pop in pops:
                        tmp = copy.copy(patch.patch_datas[pop].dataset)

                        if final[pop] is None:
                            final[pop] = tmp
                        else:
                            final[pop].add(tmp)

        # select particles
        if "select" in kwargs:
            for pop, particles in final.items():
                final[pop] = kwargs["select"](particles)

        return final, dp(final, **kwargs)


def finest_part_data(hierarchy, time=None):
    """
    returns a dict {popname : Particles}
    Particles contained in the dict are those from
    the finest patches available at a given location
    """
    from copy import deepcopy

    from ..particles import remove

    # we are going to return a dict {popname : Particles}
    # we prepare it with population names
    aPatch = hierarchy.level(0, time=time).patches[0]
    particles = {popname: None for popname in aPatch.patch_datas.keys()}

    # our strategy is to explore the hierarchy from the finest
    # level to the coarsest. at Each level we keep only particles
    # that are in cells that are not overlaped by finer boxes

    # this dict keeps boxes for patches at each level
    # each level will thus need this dict to see next finer boxes
    lvlPatchBoxes = {ilvl: [] for ilvl in range(hierarchy.finest_level(time) + 1)}

    for ilvl in range(hierarchy.finest_level(time) + 1)[::-1]:
        plvl = hierarchy.level(ilvl, time=time)
        for ip, patch in enumerate(plvl.patches):
            lvlPatchBoxes[ilvl].append(patch.box)
            for popname, pdata in patch.patch_datas.items():
                # if we're at the finest level
                # we need to keep all particles
                if ilvl == hierarchy.finest_level(time):
                    if particles[popname] is None:
                        particles[popname] = deepcopy(pdata.dataset)
                    else:
                        particles[popname].add(deepcopy(pdata.dataset))

                # if there is a finer level
                # we need to keep only those of the current patch
                # that are not in cells overlaped by finer patch boxes
                else:
                    icells = pdata.dataset.iCells
                    parts = deepcopy(pdata.dataset)
                    create = True
                    for finerBox in lvlPatchBoxes[ilvl + 1]:
                        coarseFinerBox = boxm.coarsen(finerBox, refinement_ratio)
                        within = np.where(
                            (icells >= coarseFinerBox.lower[0])
                            & (icells <= coarseFinerBox.upper[0])
                        )[0]
                        if create:
                            toRemove = within
                            create = False
                        else:
                            toRemove = np.concatenate((toRemove, within))

                    if toRemove.size != 0:
                        parts = remove(parts, toRemove)
                    if parts is not None:
                        particles[popname].add(parts)
    return particles


def amr_grid(hierarchy, time):
    """returns a non-uniform contiguous primal grid
    associated to the given hierarchy
    """
    lvlPatchBoxes = {ilvl: [] for ilvl in range(hierarchy.finest_level() + 1)}
    finalCells = {ilvl: None for ilvl in range(hierarchy.finest_level() + 1)}
    lvl = hierarchy.levels(time)

    for ilvl in range(hierarchy.finest_level(time) + 1)[::-1]:
        sorted_patches = sorted(lvl[ilvl].patches, key=lambda p: p.layout.box.lower[0])

        for ip, patch in enumerate(sorted_patches):
            box = patch.layout.box
            lvlPatchBoxes[ilvl].append(box)

            # we create a list of all cells in the current patch
            # remember that if the box upper cell is, say = 40,
            # it means that the upper node is the lower node of cell 41
            # so to get all primal nodes of a patch we need to include
            # one past the upper cell.
            # this said we do not want to include that last primal nodes
            # all the time because that would be a duplicate with the lower
            # node of the next patch. We only want to add it for the LAST
            # (because sorted) patch. We also do not want to do it on levels
            # other than L0 because the last primal node of the last patch
            # of L_i is the first primal node of a L_{i-1} node, so including it
            # would also mean adding a duplicate.
            last = 1 if ilvl == 0 and ip == len(sorted_patches) - 1 else 0
            cells = np.arange(box.lower[0], box.upper[0] + 1 + last)

            # finest level has no next finer so we take all cells
            if ilvl == hierarchy.finest_level(time):
                if finalCells[ilvl] is None:
                    finalCells[ilvl] = cells
                else:
                    finalCells[ilvl] = np.concatenate((finalCells[ilvl], cells))

            else:
                # on other levels
                # we take only grids not overlaped by next finer
                coarsenedNextFinerBoxes = [
                    boxm.coarsen(b, refinement_ratio) for b in lvlPatchBoxes[ilvl + 1]
                ]
                for coarseBox in coarsenedNextFinerBoxes:
                    ccells = np.arange(coarseBox.lower[0], coarseBox.upper[0] + 1)
                    inter, icells, iccells = np.intersect1d(
                        cells, ccells, return_indices=True
                    )
                    cells = np.delete(cells, icells)
                if len(cells):
                    if finalCells[ilvl] is None:
                        finalCells[ilvl] = cells
                    else:
                        finalCells[ilvl] = np.unique(
                            np.concatenate((finalCells[ilvl], cells))
                        )

    # now we have all cells for each level we
    # just need to compute the primal coordinates
    # and concatenate in a single array
    for ilvl in range(hierarchy.finest_level() + 1):
        if ilvl == 0:
            x = finalCells[ilvl] * hierarchy.level(ilvl).patches[0].layout.dl[0]
        else:
            xx = finalCells[ilvl] * hierarchy.level(ilvl).patches[0].layout.dl[0]
            x = np.concatenate((x, xx))

    return np.sort(x)
