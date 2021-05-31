

import os
import numpy as np

from .particles import Particles

from ..core import box as boxm
from ..core.box import Box
from ..core.gridlayout import GridLayout
import matplotlib.pyplot as plt
from ..core.phare_utilities import np_array_ify, is_scalar, listify, refinement_ratio


class PatchData:
    """
    base class for FieldData and ParticleData
    this class just factors common geometrical properties
    """
    def __init__(self, layout, quantity):
        """
        :param layout: a GridLayout representing the domain on which the data
             is defined
        :param quantity: ['field', 'particle']
        """
        self.quantity = quantity
        self.box      = layout.box
        self.origin   = layout.origin
        self.layout   = layout




class FieldData(PatchData):
    """
    Concrete type of PatchData representing a physical quantity
    defined on a grid.
    """

    @property
    def x(self):
        if self._x is None:
            self._x = self.layout.yeeCoordsFor(self.field_name, "x")
        return self._x

    @property
    def y(self):
        if self._y is None:
            self._y = self.layout.yeeCoordsFor(self.field_name, "y")
        return self._y

    @property
    def z(self):
        if self._z is None:
            self._z = self.layout.yeeCoordsFor(self.field_name, "z")
        return self._z

    def primal_directions(self):
        return self.size - self.ghost_box.shape

    def __str__(self):
        return "FieldData: (box=({}, {}), key={})".format(self.layout.box, self.layout.box.shape, self.field_name)
    def __repr__(self):
        return self.__str__()

    def select(self, box):
        """
          return view of internal data based on overlap of input box
             returns a view +1 in size in primal directions
        """
        assert isinstance(box, Box) and box.ndim == self.box.ndim

        gbox = self.ghost_box.copy()
        gbox.upper += self.primal_directions()

        box = box.copy()
        box.upper += self.primal_directions()

        overlap = box * gbox
        if overlap is not None:
            lower = self.layout.AMRIndexToLocal(dim=box.ndim - 1, index=overlap.lower)
            upper  = self.layout.AMRIndexToLocal(dim=box.ndim - 1, index=overlap.upper)
            if box.ndim == 1:
                return self.dataset[lower[0] : upper[0] + 1]
            if box.ndim == 2:
                return self.dataset[lower[0]:upper[0] + 1 , lower[1] : upper[1] + 1]
        return np.array([])


    def __getitem__(self, box):
        return self.select(box)


    def __init__(self, layout, field_name, data, **kwargs):
        """
        :param layout: A GridLayout representing the domain on which data is defined
        :param field_name: the name of the field (e.g. "Bx")
        :param data: the dataset from which data can be accessed
        """
        super().__init__(layout, 'field')
        self._x = None
        self._y = None
        self._z = None

        self.layout = layout
        self.field_name = field_name
        self.name = field_name
        self.dl = np.asarray(layout.dl)
        self.ndim = layout.box.ndim
        self.ghosts_nbr = np.zeros(self.ndim, dtype=int)

        if field_name in layout.centering["X"]:
            directions = ["X", "Y", "Z"][:layout.box.ndim] # drop unused directions
            centerings = [layout.qtyCentering(field_name, direction) for direction in directions]
        elif "centering" in kwargs:
            if isinstance(kwargs["centering"], list):
                centerings = kwargs["centering"]
                assert len(centerings) == self.ndim
            else:
                if self.ndim != 1:
                    raise ValueError("FieldData invalid dimenion for centering argument, expected list for dim > 1")
                centerings = [kwargs["centering"]]
        else:
            raise ValueError("centering not specified and cannot be inferred from field name")

        for i, centering in enumerate(centerings):
            self.ghosts_nbr[i] = layout.nbrGhosts(layout.interp_order, centering)

        self.ghost_box = boxm.grow(layout.box, self.ghosts_nbr)

        self.size = np.copy(self.ghost_box.shape)
        self.offset = np.zeros(self.ndim)

        for i, centering in enumerate(centerings):
            if centering == "primal":
                self.size[i] = self.ghost_box.shape[i] + 1
            else:
                self.size[i] = self.ghost_box.shape[i]
                self.offset[i] = 0.5*self.dl[i]

        self.dataset = data



class ParticleData(PatchData):
    """
    Concrete type of PatchData representing particles in a region
    """
    def __init__(self, layout, data, pop_name):
        """
        :param layout: A GridLayout object representing the domain in which particles are
        :param data: dataset containing particles
        """
        super().__init__(layout, 'particles')
        self.dataset = data
        self.pop_name = pop_name
        self.name = pop_name
        self.ndim = layout.box.ndim

        self.pop_name = pop_name
        if layout.interp_order == 1:
            self.ghosts_nbr = np.array([1] * layout.box.ndim)
        elif layout.interp_order == 2 or layout.interp_order == 3:
            self.ghosts_nbr = np.array([2] * layout.box.ndim)
        else:
            raise RuntimeError("invalid interpolation order {}".format(layout.interp_order))

        self.ghost_box = boxm.grow(layout.box, self.ghosts_nbr)
        assert (self.box.lower == self.ghost_box.lower + self.ghosts_nbr).all()


class Patch:
    """
    A patch represents a hyper-rectangular region of space
    """

    def __init__(self, patch_datas):
        """
        :param patch_datas: a list of PatchData objects
        these are assumed to "belong" to the Patch so to
        share the same origin, mesh size and box.
        """
        pdata0 = list(patch_datas.values())[0] #0 represents all others
        self.layout = pdata0.layout
        self.box = pdata0.layout.box
        self.origin = pdata0.layout.origin
        self.dl = pdata0.layout.dl
        self.patch_datas = patch_datas




class PatchLevel:
    """is a collection of patches """

    def __init__(self, lvl_nbr, patches):
        self.level_number = lvl_nbr
        self.patches = patches

    def __iter__(self):
        return self.patches.__iter__()

    def level_range(self):
        name = list(self.patches[0].patch_datas.keys())[0]
        return min([patch.patch_datas[name].x.min() for patch in self.patches]),\
                max([patch.patch_datas[name].x.max() for patch in self.patches])



def are_adjacent(lower, upper, atol=1e-6):
    return np.abs(upper[0]-lower[-1]) < atol


def overlap_mask(x, level, qty):
    """
    returns a mask for x where x is overlaped by
    the qty patch datas on the specified level
    """

    is_overlaped = x != x #lgtm [py/comparison-of-identical-expressions]

    # the strategy is to exclude from this patch data
    # all nodes and values that are overlaped by next finer values
    # the is_overlaped array will serve as a mask

    for patch in level.patches:
        pdata = patch.patch_datas[qty]
        ghosts_nbr = pdata.ghosts_nbr[0]
        finex = pdata.x[ghosts_nbr:-ghosts_nbr]
        fine_dx = finex[1]-finex[0]
        local_dx = x[1]-x[0]
        if fine_dx < local_dx:
            # the fine patch min and max are widened
            # a bit so we include border nodes as overlaped
            xmin, xmax = finex.min() - fine_dx/2, finex.max() + fine_dx/2
            overlaped_idx = np.where( (x > xmin) & (x < xmax))[0]
            is_overlaped[overlaped_idx] = True
        else:
            raise ValueError("level needs to have finer grid resolution than that of x")

    return is_overlaped






def finest_field(hierarchy, qty, time=None):
    """returns a non-uniform contiguous (value,grid)
       associated to the given hierarchy
       the quantity "qty" must be a field
    """
    lvl = hierarchy.levels(time)

    for ilvl in range(hierarchy.finest_level(time)+1)[::-1]:
        sorted_patches = sorted(lvl[ilvl].patches,
                                key= lambda p:p.box.lower[0])

        for ip, patch in enumerate(sorted_patches):

            pdata = patch.patch_datas[qty]
            ghosts_nbr = pdata.ghosts_nbr[0]
            patch_values = pdata.dataset[ghosts_nbr:-ghosts_nbr]
            x = pdata.x[ghosts_nbr:-ghosts_nbr]

            if ilvl == hierarchy.finest_level(time):
                if ip == 0:
                    final = patch_values
                    final_x = x
                    last_x = x
                else:
                    if are_adjacent(last_x, x):
                        first_node = 1
                    else:
                        first_node = 0

                    final = np.concatenate((final, patch_values[first_node:]))
                    final_x = np.concatenate((final_x, x[first_node:]))
                    last_x = x[first_node:]

            else:
                is_overlaped = overlap_mask(x, hierarchy.level(ilvl+1, time), qty)
                unique_values = patch_values[~is_overlaped]
                unique_x = x[~is_overlaped]

                if unique_x.size:
                    if are_adjacent(last_x, unique_x):
                        first_node = 1
                    else:
                        first_node = 0

                    last_x = unique_x[first_node:]

                    tmp_x = np.concatenate((final_x, unique_x[first_node:]))
                    final_x = np.sort(tmp_x)

                    # ok now we need to insert unique_values in
                    # final. unique_values should be inserted in final
                    # in the same order unique_x should be in final_x
                    # i.e. with the indexing that sorts tmp_x

                    order = np.argsort(tmp_x)
                    tmp = np.concatenate((final, unique_values[first_node:-1]))
                    final = np.concatenate((final, unique_values[first_node:]))[order]

    return final, final_x




def finest_part_data(hierarchy, time=None):
    """
    returns a dict {popname : Particles}
    Particles contained in the dict are those from
    the finest patches available at a given location
    """
    from .particles import remove
    from copy import deepcopy

    # we are going to return a dict {popname : Particles}
    # we prepare it with population names
    aPatch = hierarchy.level(0, time=time).patches[0]
    particles = {popname:None for popname in aPatch.patch_datas.keys()}

    # our strategy is to explore the hierarchy from the finest
    # level to the coarsest. at Each level we keep only particles
    # that are in cells that are not overlaped by finer boxes

    # this dict keeps boxes for patches at each level
    # each level will thus need this dict to see next finer boxes
    lvlPatchBoxes = {ilvl:[] for ilvl in range(hierarchy.finest_level(time)+1)}

    for ilvl in range(hierarchy.finest_level(time)+1)[::-1]:
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
                    for finerBox in lvlPatchBoxes[ilvl+1]:
                        coarseFinerBox = boxm.coarsen(finerBox, refinement_ratio)
                        within = np.where((icells >= coarseFinerBox.lower[0]) &\
                                                 (icells <= coarseFinerBox.upper[0]))[0]
                        if create:
                            toRemove = within
                            create=False
                        else:
                            toRemove = np.concatenate((toRemove, within))

                    parts = remove(parts, toRemove)
                    if parts is not None:
                        particles[popname].add(parts)
    return particles

class PatchHierarchy:
    """is a collection of patch levels """


    def __init__(self, patch_levels, domain_box, refinement_ratio=2, time=0., data_files=None):
        self.patch_levels = patch_levels
        self.ndim = len(domain_box.lower)
        self.time_hier = {}
        self.time_hier.update({self.format_timestamp(time):patch_levels})

        self.domain_box = domain_box
        self.refinement_ratio = refinement_ratio

        self.data_files = {}

        if data_files is not None:
            self.data_files.update(data_files)

    def _default_time(self):
        return self.times()[0]

    def finest_level(self, time=None):
        if time is None:
            time  = self._default_time()
        return max(list(self.levels(time=time).keys()))


    def levels(self, time=None):
        if time is None:
            time = self._default_time()
        return self.time_hier[self.format_timestamp(time)]


    def level(self, level_number, time=None):
        return self.levels(time)[level_number]


    def levelNbr(self, time):
        return len(self.levels(time).items())

    def levelNbrs(self, time):
        return list(self.levels(time).keys())

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
                    it_is &= qties == list(patch.patch_datas.keys())
        return it_is

    def _quantities(self):
        return list(self.level(0).patches[0].patch_datas.keys())

    def quantities(self):
        if not self.is_homogeneous():
            raise RuntimeError("Error - hierarchy is not homogeneous")
        return self._quantities()

    def refined_domain_box(self, level_number):
        """
        returns the domain box refined for a given level number
        """
        assert (level_number >= 0)
        return boxm.refine(self.domain_box, self.refinement_ratio ** level_number)


    def format_timestamp(self, timestamp):
        if isinstance(timestamp, str):
            return timestamp
        return "{:.10f}".format(timestamp)


    def level_domain_box(self, level_number):
        if level_number == 0:
            return self.domain_box
        return self.refined_domain_box(level_number)


    def __str__(self):
        s = "Hierarchy: \n"
        for t, patch_levels in self.time_hier.items():
            for ilvl, lvl in patch_levels.items():
                s = s + "Level {}\n".format(ilvl)
                for ip, patch in enumerate(lvl.patches):
                    for qty_name, pd in patch.patch_datas.items():
                        pdstr = "    P{ip} {type} {pdname} box is {box} and ghost box is {gbox}"
                        s = s + pdstr.format(ip=ip, type=type(pd.dataset), pdname=qty_name,
                                             box=patch.box, gbox=pd.ghost_box)
                        s = s + "\n"
        return s

    def times(self):
        return np.sort(np.asarray(list(self.time_hier.keys())))


    def plot_patches(self, save=False):
        fig, ax = plt.subplots(figsize=(10, 3))
        for ilvl, lvl in self.levels(0.).items():
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
        if isinstance(collections, list) and all([isinstance(el, Box) for el in collections]):
            collections = [{"boxes" : collections}]

        from matplotlib.collections import PatchCollection

        level_domain_box = self.level_domain_box(ilvl)
        mi, ma = level_domain_box.lower.min(), level_domain_box.upper.max()

        fig, ax = kwargs.get("subplot", plt.subplots(figsize=(6, 6)))

        for collection in collections:
            facecolor = collection.get("facecolor", "none")
            edgecolor = collection.get("edgecolor", 'purple')
            alpha = collection.get("alpha", 1)
            rects = [self.box_to_Rectangle(box) for box in collection["boxes"]]

            ax.add_collection(PatchCollection(rects,
                              facecolor=facecolor,
                              alpha=alpha,
                              edgecolor=edgecolor))

        if "title" in kwargs:
            from textwrap import wrap
            xfigsize = int(fig.get_size_inches()[0] * 10) # 10 characters per inch
            ax.set_title("\n".join(wrap(kwargs["title"], xfigsize)))

        major_ticks = np.arange(mi - 5, ma + 5 + 5, 5)
        ax.set_xticks(major_ticks)
        ax.set_yticks(major_ticks)

        minor_ticks = np.arange(mi - 5, ma + 5 + 5, 1)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(minor_ticks, minor=True)

        ax.grid(which='both')

        return fig



    def plot(self, **kwargs):
        """
        plot
        """
        usr_lvls = kwargs.get("levels",(0,))
        qty = kwargs.get("qty",None)
        time = kwargs.get("time", self.times()[0])

        if "ax" not in kwargs:
            fig, ax = plt.subplots()
        else:
            ax = kwargs["ax"]
            fig = ax.figure
        for lvl_nbr,level in self.levels(time).items():
            if lvl_nbr not in usr_lvls:
                continue
            for ip, patch in enumerate(level.patches):
                pdata_nbr = len(patch.patch_datas)
                if qty is None and pdata_nbr != 1:
                    multiple = "multiple quantities in patch, "
                    err = multiple + "please specify a quantity in  "+ pdata_names
                    raise ValueError(err)
                if qty is None:
                    qty = list(patch.patch_datas.keys())[0]

                val = patch.patch_datas[qty].dataset[:]
                x   = patch.patch_datas[qty].x
                label = "L{level}P{patch}".format(level=lvl_nbr,patch=ip)
                marker=kwargs.get("marker", "")
                ls = kwargs.get("ls","--")
                ax.plot(x, val, label=label, marker=marker, ls=ls)

        ax.set_title(kwargs.get("title",""))
        ax.set_xlabel(kwargs.get("xlabel","x"))
        ax.set_ylabel(kwargs.get("ylabel", qty))
        if "xlim" in kwargs:
            ax.set_xlim(kwargs["xlim"])
        if "ylim" in kwargs:
            ax.set_ylim(kwargs["ylim"])

        if kwargs.get("legend", None) is not None:
            ax.legend()

        if "filename" in kwargs:
            fig.savefig(kwargs["filename"])


    def dist_plot(self, **kwargs):
        """
        plot phase space of a particle hierarchy
        """
        import copy
        from .plotting import dist_plot as dp
        usr_lvls = kwargs.get("levels",(0,))
        finest = kwargs.get("finest", False)
        pops = kwargs.get("pop",[])
        time = kwargs.get("time", self.times()[0])
        axis = kwargs.get("axis", ("Vx", "Vy"))
        all_pops = list(self.level(0,time).patches[0].patch_datas.keys())

        vmin = kwargs.get("vmin", -2)
        vmax = kwargs.get("vmax", 2)
        dv = kwargs.get("dv", 0.05)
        vbins = vmin + dv*np.arange(int((vmax-vmin)/dv))

        if finest:
            final = finest_part_data(self)
            if axis[0] == "x":
                xbins = amr_grid(self, time)
                bins = (xbins, vbins)
            else:
                bins=(vbins, vbins)
            kwargs["bins"] = bins

        else:
            final = {pop:None for pop in all_pops}
            for lvl_nbr,level in self.levels(time).items():
                if lvl_nbr not in usr_lvls:
                    continue
                for ip, patch in enumerate(level.patches):
                    if len(pops)==0:
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





def amr_grid(hierarchy, time):
    """returns a non-uniform contiguous primal grid
       associated to the given hierarchy
    """
    lvlPatchBoxes = {ilvl:[] for ilvl in range(hierarchy.finest_level()+1)}
    finalCells = {ilvl:None for ilvl in range(hierarchy.finest_level()+1)}
    lvl = hierarchy.levels(time)

    for ilvl in range(hierarchy.finest_level(time)+1)[::-1]:

        sorted_patches = sorted(lvl[ilvl].patches,
                                key= lambda p:p.layout.box.lower[0])

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
            last = 1 if ilvl==0 and ip == len(sorted_patches)-1 else 0
            cells = np.arange(box.lower[0], box.upper[0]+1 + last)

            # finest level has no next finer so we take all cells
            if ilvl == hierarchy.finest_level(time):
                if finalCells[ilvl] is None:
                    finalCells[ilvl] = cells
                else:
                    finalCells[ilvl] = np.concatenate((finalCells[ilvl],cells))

            else:
                # on other levels
                # we take only grids not overlaped by next finer
                coarsenedNextFinerBoxes = [boxm.coarsen(b, refinement_ratio) for b in lvlPatchBoxes[ilvl+1]]
                for coarseBox in coarsenedNextFinerBoxes:
                    ccells = np.arange(coarseBox.lower[0], coarseBox.upper[0]+1)
                    inter,icells, iccells = np.intersect1d(cells, ccells, return_indices=True)
                    cells = np.delete(cells, icells)
                if len(cells):
                    if finalCells[ilvl] is None:
                        finalCells[ilvl] = cells
                    else:
                        finalCells[ilvl] = np.unique(np.concatenate((finalCells[ilvl],cells)))

    # now we have all cells for each level we
    # just need to compute the primal coordinates
    # and concatenate in a single array
    for ilvl in range(hierarchy.finest_level()+1):
        if ilvl == 0:
            x = finalCells[ilvl]*hierarchy.level(ilvl).patches[0].layout.dl[0]
        else:
            xx = finalCells[ilvl]*hierarchy.level(ilvl).patches[0].layout.dl[0]
            x = np.concatenate((x, xx))

    return np.sort(x)



def is_root_lvl(patch_level):
    return patch_level.level_number == 0








field_qties = {"EM_B_x": "Bx",
               "EM_B_y": "By",
               "EM_B_z": "Bz",
               "EM_E_x": "Ex",
               "EM_E_y": "Ey",
               "EM_E_z": "Ez",
               "flux_x": "Fx",
               "flux_y": "Fy",
               "flux_z": "Fz",
               "bulkVelocity_x": "Vx",
               "bulkVelocity_y": "Vy",
               "bulkVelocity_z": "Vz",
               "density": "rho"}


particle_files_patterns = ("domain", "patchGhost", "levelGhost")


def is_particle_file(filename):
    return any([pattern in filename for pattern in particle_files_patterns])



def particle_dataset_name(basename):
    """
    return "alpha_domain" from ions_pop_alpha_domain.h5
    """
    popname = basename.strip(".h5").split("_")[-2]
    part_type = basename.strip(".h5").split("_")[-1]
    dataset_name = popname + "_" + part_type

    return dataset_name


def is_pop_fluid_file(basename):
    return (is_particle_file(basename) is False) and "pop" in basename




def make_layout(h5_patch_grp, cell_width, interp_order):
    origin = h5_patch_grp.attrs['origin']
    upper = h5_patch_grp.attrs['upper']
    lower = h5_patch_grp.attrs['lower']
    return GridLayout(Box(lower, upper), origin, cell_width, interp_order=interp_order)



def create_from_all_times(time, hier):
    return time is None and hier is None


def create_from_one_time(time, hier):
    return time is not None and hier is None


def load_all_times(time, hier):
    return time is None and hier is not None


def load_one_time(time, hier):
    return time is not None and hier is not None


def compute_hier_from(h, compute):
    """
    returns a hierarchy resulting from calling 'compute'
    on each patch of the given hierarchy 'h'

    compute is a function taking a Patch and returning
    a list of dicts with the following keys:

        "data": ndarray containing the data
        "name": str, name of the data living on that patch, must be unique
        "centering": str, ["dual", "primal"]

     caveat: routine only works in 1D so far.
    """
    patch_levels = {}
    for ilvl, lvl in h.patch_levels.items():
        patches = {}
        for ip, patch in enumerate(lvl.patches):
            new_patch_datas = {}
            layout = patch.layout
            datas = compute(patch)
            for data in datas:
                pd = FieldData(layout, data["name"],
                               data["data"],
                               centering=data["centering"])
                new_patch_datas[data["name"]] = pd
            if ilvl not in patches:
                patches[ilvl] = []
            patches[ilvl].append(Patch(new_patch_datas))

        patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

    return PatchHierarchy(patch_levels, h.domain_box, refinement_ratio)




def pop_name(basename):
    return basename.strip(".h5").split("_")[-2]


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
            dl[:,i] = layout.dl[i]

        particles = Particles(icells=h5_patch_grp["iCell"],
                              deltas=h5_patch_grp["delta"],
                              v=v,
                              weights=h5_patch_grp["weight"],
                              charges=h5_patch_grp["charge"],
                              dl=dl)

        pdname = particle_dataset_name(basename)
        if pdname in patch_datas:
            raise ValueError("error - {} already in patchdata".format(pdname))

        patch_datas[pdname] = ParticleData(layout, particles, pop_name(basename))

    else:
        for dataset_name in h5_patch_grp.keys():

            dataset = h5_patch_grp[dataset_name]

            if dataset_name not in field_qties:
                raise RuntimeError(
                    "invalid dataset name : {} is not in {}".format(dataset_name, field_qties))

            pdata = FieldData(layout, field_qties[dataset_name], dataset)

            pdata_name = field_qties[dataset_name]

            if is_pop_fluid_file(basename):
                pdata_name = pop_name(basename) + "_" + pdata_name


            if dataset_name in patch_datas:
                raise ValueError("error - {} already in patchdata".format(dataset_name))

            patch_datas[pdata_name] = pdata

    return True # valid patch assumed



def patch_has_datasets(h5_patch_grp):
    return len(h5_patch_grp.keys())>0



h5_time_grp_key = "t"

def hierarchy_fromh5(h5_filename, time, hier, silent=True):
    import h5py
    data_file = h5py.File(h5_filename, "r")
    basename = os.path.basename(h5_filename)

    root_cell_width = np.asarray(data_file.attrs["cell_width"])
    interp = data_file.attrs["interpOrder"]
    domain_box = Box([0] * len(data_file.attrs["domain_box"]), data_file.attrs["domain_box"])

    if create_from_all_times(time, hier):
        # first create from first time
        # then add all other times
        if not silent:
            print("creating hierarchy from all times in file")
        times = list(data_file[h5_time_grp_key].keys())
        hier = hierarchy_fromh5(h5_filename, time=times[0], hier=hier)
        if len(times) > 1:
            for t in times[1:]:
                hierarchy_fromh5(h5_filename, t, hier)
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
            lvl_cell_width = root_cell_width / refinement_ratio ** ilvl
            patches = {}

            for pkey in h5_patch_lvl_grp.keys():

                h5_patch_grp = h5_patch_lvl_grp[pkey]

                if patch_has_datasets(h5_patch_grp):
                    patch_datas = {}
                    layout = make_layout(h5_patch_grp, lvl_cell_width, interp)
                    add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)

                    if ilvl not in patches:
                        patches[ilvl] = []

                    patches[ilvl].append(Patch(patch_datas))

                    patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

        diag_hier = PatchHierarchy(patch_levels, domain_box, refinement_ratio, t, data_file)

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
                lvl_cell_width = root_cell_width / refinement_ratio ** ilvl

                for ipatch, pkey in enumerate(h5_time_grp[plvl_key].keys()):
                    h5_patch_grp = h5_time_grp[plvl_key][pkey]

                    if patch_has_datasets(h5_patch_grp):
                        hier_patch = patch_levels[ilvl].patches[ipatch]
                        origin = h5_time_grp[plvl_key][pkey].attrs['origin']
                        upper = h5_time_grp[plvl_key][pkey].attrs['upper']
                        lower = h5_time_grp[plvl_key][pkey].attrs['lower']
                        file_patch_box = Box(lower, upper)

                        assert file_patch_box == hier_patch.box
                        assert (abs(origin - hier_patch.origin) < 1e-6).all()
                        assert (abs(lvl_cell_width - hier_patch.dl) < 1e-6).all()

                        layout = make_layout(h5_patch_grp, lvl_cell_width, interp)
                        add_to_patchdata(hier_patch.patch_datas, h5_patch_grp, basename, layout)

            return hier

        if not silent:
            print("adding data to new time")
        # time does not exist in the hierarchy
        # we have to create a brand new set of patchLevels
        # containing patches, and load data in their patchdatas

        patch_levels = {}

        for plvl_key in h5_time_grp.keys():
            ilvl = int(plvl_key[2:])

            lvl_cell_width = root_cell_width / refinement_ratio ** ilvl
            lvl_patches = []

            for ipatch, pkey in enumerate(h5_time_grp[plvl_key].keys()):
                h5_patch_grp = h5_time_grp[plvl_key][pkey]

                if patch_has_datasets(h5_patch_grp):
                    layout = make_layout(h5_patch_grp, lvl_cell_width, interp)
                    patch_datas = {}
                    add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)
                    lvl_patches.append(Patch(patch_datas))

            patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

        hier.time_hier[t] = patch_levels
        return hier

    if load_all_times(time, hier):
        if not silent:
            print("loading all times in existing hier")
        for time in data_file[h5_time_grp_key].keys():
            hier = hierarchy_fromh5(h5_filename, time, hier)

        return hier


def quantidic(ilvl, wrangler):
    pl = wrangler.getPatchLevel(ilvl)

    return  {"density": pl.getDensity,
             "bulkVelocity_x":pl.getVix,
             "bulkVelocity_y": pl.getViy,
             "bulkVelocity_z":pl.getViz,
             "EM_B_x":pl.getBx,
             "EM_B_y": pl.getBy,
             "EM_B_z": pl.getBz,
             "EM_E_x": pl.getEx,
             "EM_E_y": pl.getEy,
             "EM_E_z": pl.getEz,
             "flux_x": pl.getFx,
             "flux_y": pl.getFy,
             "flux_z": pl.getFz,
             "particles":pl.getParticles}



def isFieldQty(qty):
    return qty in ("density",
                   "bulkVelocity_x",
                   "bulkVelocity_y",
                   "bulkVelocity_z",
                   "EM_B_x",
                   "EM_B_y",
                   "EM_B_z",
                   "EM_E_x",
                   "EM_E_y",
                   "EM_E_z",
                   "flux_x", "flux_y", "flux_z")



def hierarchy_from_sim(simulator, qty, pop=""):
    dw = simulator.data_wrangler()
    nbr_levels = dw.getNumberOfLevels()
    patch_levels = {}

    root_cell_width = simulator.cell_width()
    domain_box = Box([0] * len(root_cell_width) , simulator.domain_box())
    assert len(domain_box.ndim) == len(simulator.domain_box().ndim)

    for ilvl in range(nbr_levels):

        lvl_cell_width = root_cell_width / refinement_ratio ** ilvl

        patches = {ilvl : [] for ilvl in range(nbr_levels)}
        getters = quantidic(ilvl, dw)


        if isFieldQty(qty):
            wpatches = getters[qty]()
            for patch in wpatches:
                patch_datas = {}
                lower = patch.lower
                upper = patch.upper
                origin = patch.origin
                layout = GridLayout(Box(lower, upper), origin, lvl_cell_width, interp_order = simulator.interporder())
                pdata = FieldData(layout, field_qties[qty], patch.data)
                patch_datas[qty] = pdata
                patches[ilvl].append(Patch(patch_datas))

        elif qty == "particles":

            if pop=="":
                raise ValueError("must specify pop argument for particles")
            # here the getter returns a dict like this
            # {'protons': {'patchGhost': [<pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78970>,
            #<pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78f70>],
            # 'domain': [<pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78d70>,
            # <pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78770>]}}

            # domain particles are assumed to always be here
            # but patchGhost and levelGhost may not be, depending on the level

            populationdict = getters[qty](pop)[pop]


            dom_dw_patches = populationdict["domain"]
            for patch in dom_dw_patches:
                patch_datas= {}

                lower = patch.lower
                upper = patch.upper
                origin = patch.origin
                layout = GridLayout(Box(lower, upper), origin, lvl_cell_width, interp_order=simulator.interp_order())
                v = np.asarray(patch.data.v).reshape(int(len(patch.data.v) / 3), 3)

                domain_particles = Particles(icells = np.asarray(patch.data.iCell),
                                             deltas = np.asarray(patch.data.delta),
                                             v      = v,
                                             weights = np.asarray(patch.data.weight),
                                             charges = np.asarray(patch.data.charge))

                patch_datas[pop +"_particles"] = ParticleData(layout, domain_particles, pop)
                patches[ilvl].append(Patch(patch_datas))



            # ok now let's add the patchGhost if present
            # note that patchGhost patches may not be the same list as the
            # domain patches... since not all patches may not have patchGhost while they do have
            # domain... while looping on the patchGhost items, we need to search in
            # the already created patches which one to which add the patchGhost particles


            for ghostParticles in ["patchGhost", "levelGhost"]:
                if ghostParticles in populationdict:
                    for dwpatch in populationdict[ghostParticles]:

                        v = np.asarray(dwpatch.data.v)
                        s = v.size
                        v = v[:].reshape(int(s / 3), 3)


                        patchGhost_part = Particles(icells=np.asarray(dwpatch.data.iCell),
                                                     deltas=np.asarray(dwpatch.data.delta),
                                                     v= v,
                                                     weights=np.asarray(dwpatch.data.weight),
                                                     charges=np.asarray(dwpatch.data.charge))

                        box = Box(dwpatch.lower, dwpatch.upper)

                        # now search which of the already created patches has the same box
                        # once found we add the new particles to the ones already present

                        patch = [p for p in patches[ilvl] if p.box == box][0]
                        patch.patch_datas[pop+"_particles"].dataset.add(patchGhost_part)






        else:
            raise ValueError("{} is not a valid quantity".format(qty))



        patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

    return PatchHierarchy(patch_levels, domain_box, time=simulator.currentTime())









def hierarchy_from(simulator=None, qty= None, pop = "", h5_filename=None, time=None, hier=None):
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
        return hierarchy_fromh5(h5_filename, time, hier)

    if simulator is not None and qty is not None:
        return hierarchy_from_sim(simulator, qty, pop=pop)

    raise ValueError("can't make hierarchy")



def merge_particles(hierarchy):

    for time, patch_levels in hierarchy.time_hier.items():
        for ilvl, plvl in patch_levels.items():
            for ip, patch in enumerate(plvl.patches):
                pdatas = patch.patch_datas
                domain_pdata = [(pdname,pd) for pdname, pd in pdatas.items() if "domain" in pdname][0]

                pghost_pdatas = [(pdname,pd) for pdname, pd in pdatas.items() if "patchGhost" in pdname]
                lghost_pdatas = [(pdname,pd) for pdname, pd in pdatas.items() if "levelGhost" in pdname]

                pghost_pdata = pghost_pdatas[0] if pghost_pdatas else None
                lghost_pdata = lghost_pdatas[0] if lghost_pdatas else None

                if pghost_pdata is not None:
                    domain_pdata[1].dataset.add(pghost_pdata[1].dataset)
                    del pdatas[pghost_pdata[0]]

                if lghost_pdata is not None:
                    domain_pdata[1].dataset.add(lghost_pdata[1].dataset)
                    del pdatas[lghost_pdata[0]]

                popname = domain_pdata[0].split('_')[0]
                pdatas[popname+"_particles"] = pdatas[domain_pdata[0]]
                del pdatas[domain_pdata[0]]


def h5_filename_from(diagInfo):
    # diagInfo.quantity starts with a / , hence   [1:]
    return (diagInfo.quantity + ".h5").replace('/', '_')[1:]

def get_times_from_h5(filepath):
    import h5py
    f = h5py.File(filepath, 'r')
    times = np.array(sorted([float(s) for s in list(f["t"].keys())]))
    f.close()
    return times
