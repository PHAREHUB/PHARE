
from ..core import box as boxm
from ..core.box import Box
import numpy as np
import os
from ..core.gridlayout import GridLayout
from .particles import Particles
from ..data.wrangler import DataWrangler
import matplotlib.pyplot as plt
import seaborn as sns
from ..core.phare_utilities import np_array_ify, is_scalar


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
    def _mesh_coords(self, i):
        return self.origin[i] - self.ghosts_nbr * self.dl[i] + np.arange(self.size[i]) * self.dl[i] + self.offset[i]

    @property
    def x(self):
        if self._x is None:
            self._x = self._mesh_coords(0)
        return self._x

    def __init__(self, layout, field_name, data, **kwargs):
        """
        :param layout: A GridLayout representing the domain on which data is defined
        :param field_name: the name of the field (e.g. "Bx")
        :param data: the dataset from which data can be accessed
        """
        super().__init__(layout, 'field')
        self._x = None

        self.layout = layout
        self.field_name = field_name
        self.dx = layout.dl[0] # dropped in 2d_py_init PR - use dl[0]
        self.dl = layout.dl
        self.dim = layout.box.dim()
        self.ghosts_nbr = np.zeros(self.dim)

        if field_name in layout.centering["X"]:
            directions = ["X", "Y", "Z"][:self.dim] # drop unused directions
            centerings = [layout.qtyCentering(field_name, direction) for direction in directions]
        elif "centering" in kwargs:
            if isinstance(kwargs["centering"], list):
                centerings = kwargs["centering"]
                assert len(centerings) == self.dim
            else:
                if self.dim != 1:
                    raise ValeuError("FieldData invalid dimenion for centering argument, expected list for dim > 1")
                centerings = [kwargs["centering"]]
        else:
            ValueError("centering not specified and cannot be inferred from field name")

        for i, centering in enumerate(centerings):
            self.ghosts_nbr[i] = layout.nbrGhosts(layout.interp_order, centering)

        self.ghost_box = boxm.grow(layout.box, self.ghosts_nbr)

        self.size = np.zeros(self.dim, dtype=int)
        self.offset = np.zeros(self.dim)

        for i, centering in enumerate(centerings):
            if centering == "primal":
                self.size[i] = self.ghost_box.shape()[i] + 1
            else:
                self.size[i] = self.ghost_box.shape()[i]
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
        if layout.interp_order == 1:
            self.ghosts_nbr = [1] * layout.box.dim()
        elif layout.interp_order == 2 or layout.interp_order == 3:
            self.ghosts_nbr = [2] * layout.box.dim()
        else:
            raise RuntimeError("invalid interpolation order {}".format(layout.interp_order))
        self.ghost_box = boxm.grow(layout.box, self.ghosts_nbr)



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
        self.dx = pdata0.layout.dl[0]
        self.patch_datas = patch_datas


class PatchLevel:
    """is a collection of patches """

    def __init__(self, lvl_nbr, patches):
        self.level_number = lvl_nbr
        self.patches = patches


    def level_range(self):
        name = list(self.patches[0].patch_datas.keys())[0]
        return min([patch.patch_datas[name].x.min() for patch in self.patches]),\
                max([patch.patch_datas[name].x.max() for patch in self.patches])


class PatchHierarchy:
    """is a collection of patch levels """


    def __init__(self, patch_levels, domain_box, refinement_ratio=2, time=0., data_files=None):
        self.patch_levels = patch_levels
        self.dim = len(domain_box.lower)
        self.time_hier = {}
        self.time_hier.update({time:patch_levels})

        self.domain_box = domain_box
        self.refinement_ratio = refinement_ratio

        self.data_files = {}

        if data_files is not None:
            self.data_files.update(data_files)


    def level(self, level_number, time=0.):
        return self.time_hier[time][level_number]



    def levels(self, time=0.):
        return self.time_hier[time]


    def refined_domain_box(self, level_number):
        """
        returns the domain box refined for a given level number
        """
        assert (level_number >= 0)
        return boxm.refine(self.domain_box, self.refinement_ratio ** level_number)

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

    def plot_patches(self):
        fig, ax = plt.subplots(figsize=(10, 3))
        for ilvl, lvl in self.time_hier[0.].items():
            lvl_offset = ilvl * 0.1
            for patch in lvl.patches:
                dx = patch.dx
                x0 = patch.box.lower * dx
                x1 = patch.box.upper * dx
                xcells = np.arange(x0, x1 + dx, dx)
                y = lvl_offset + np.zeros_like(xcells)
                ax.plot(xcells, y, marker=".")

        fig.savefig("hierarchy.png")


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
                    pdata_names = "("+",".join(["'{}'".format(l) for l in patch.patch_datas])+")"
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
        ax.legend()

        if "filename" in kwargs:
            fig.savefig(filename)


    def dist_plot(self, **kwargs):
        """
        plot
        """

        usr_lvls = kwargs.get("levels",(0,))
        qty = kwargs.get("qty",None)
        time = kwargs.get("time", self.times()[0])
        cpt = 0
        fig, ax = plt.subplots()
        for lvl_nbr,level in self.levels(time).items():
            if lvl_nbr not in usr_lvls:
                continue
            for ip, patch in enumerate(level.patches):

                pdata_nbr = len(patch.patch_datas)
                if qty is None and pdata_nbr != 1:
                    pdata_names = "("+",".join(["'{}'".format(l) for l in patch.patch_datas])+")"
                    multiple = "multiple quantities in patch, "
                    err = multiple + "please specify a quantity in  "+ pdata_names
                    raise ValueError(err)
                if qty is None:
                    qty = list(patch.patch_datas.keys())[0]

                tmp = patch.patch_datas[qty].dataset

                # select particles
                if "select" in kwargs:
                    selected = kwargs["select"](tmp)
                else:
                    selected = tmp
                if cpt == 0:
                    final = selected
                else:
                    final.add(selected)

                cpt += 1

        axis = kwargs.get("axis",("Vx","Vy"))

        vaxis = {"Vx":0, "Vy":1, "Vz":2}
        if axis[0] in vaxis:
            x = final.v[:,vaxis[axis[0]]]
        elif axis[0] == "x":
            x = final.x
        if axis[1] in vaxis:
            y = final.v[:,vaxis[axis[1]]]

        bins = kwargs.get("bins", (50,50))
        ax.hist2d(x, y,
                  bins=kwargs.get("bins", bins),
                 cmap='jet')

        if kwargs.get("kde",False) is True:
            sns.kdeplot(x, y, ax=ax, color="w")


        ax.set_title(kwargs.get("title",""))
        ax.set_xlabel(kwargs.get("xlabel", axis[0]))
        ax.set_ylabel(kwargs.get("ylabel", axis[1]))
        if "xlim" in kwargs:
            ax.set_xlim(kwargs["xlim"])
        if "ylim" in kwargs:
            ax.set_ylim(kwargs["ylim"])
        ax.legend()

        if "bulk" in kwargs:
            if kwargs["bulk"] is True:
                if axis[0] in vaxis:
                    ax.axvline(final.v[:,vaxis[axis[0]]].mean(), color="w",ls="--")
                if axis[1] in vaxis:
                    ax.axhline(final.v[:,vaxis[axis[1]]].mean(), color="w",ls="--")


        if "filename" in kwargs:
            fig.savefig(filename)








def is_root_lvl(patch_level):
    return patch_level.level_number == 0








field_qties = {"EM_B_x": "Bx",
               "EM_B_y": "By",
               "EM_B_z": "Bz",
               "EM_E_x": "Ex",
               "EM_E_y": "Ey",
               "EM_E_z": "Ez",
               "flux_x": "Vx",
               "flux_y": "Vy",
               "flux_z": "Vz",
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



def make_layout(h5_patch_grp, cell_width):
    nbrCells = h5_patch_grp.attrs['nbrCells']
    origin = float(h5_patch_grp.attrs['origin'])
    upper = int(h5_patch_grp.attrs['upper'])
    lower = int(h5_patch_grp.attrs['lower'])

    return GridLayout(Box(lower, upper), origin, cell_width)



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

    return PatchHierarchy(patch_levels, h.domain_box, 2)




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

        particles = Particles(icells=h5_patch_grp["iCell"],
                              deltas=h5_patch_grp["delta"],
                              v=v,
                              weights=h5_patch_grp["weight"],
                              charges=h5_patch_grp["charge"],
                              dl=np_array_ify(layout.dl))

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

            if is_pop_fluid_file(basename):
                dataset_name = pop_name(basename) + "_" + dataset_name


            if dataset_name in patch_datas:
                raise ValueError("error - {} already in patchdata".format(dataset_name))

            patch_datas[dataset_name] = pdata

    return True # valid patch assumed



def patch_has_datasets(h5_patch_grp):
    return len(h5_patch_grp.keys())>0


def hierarchy_fromh5(h5_filename, time, hier, silent=True):
    import h5py
    data_file = h5py.File(h5_filename, "r")
    basename = os.path.basename(h5_filename)
    root_cell_width = float(data_file.attrs["cell_width"])
    domain_box = Box(0, int(data_file.attrs["domain_box"]))

    if create_from_all_times(time, hier):
        # first create from first time
        # then add all other times
        if not silent:
            print("creating hierarchy from all times in file")
        times = list(data_file.keys())
        hier = hierarchy_fromh5(h5_filename, time=times[0], hier=hier)
        if len(times) > 1:
            for t in times[1:]:
                hierarchy_fromh5(h5_filename, t, hier)
        return hier

    if create_from_one_time(time, hier):
        if not silent:
            print("creating hierarchy from time {}".format(time))
        t = float(time.strip("t"))

        h5_time_grp = data_file[time]
        patch_levels = {}

        for plvl_key in h5_time_grp.keys():

            h5_patch_lvl_grp = h5_time_grp[plvl_key]
            ilvl = int(plvl_key[2:])
            lvl_cell_width = root_cell_width / 2 ** ilvl
            patches = {}

            for pkey in h5_patch_lvl_grp.keys():

                h5_patch_grp = h5_patch_lvl_grp[pkey]

                if patch_has_datasets(h5_patch_grp):
                    patch_datas = {}
                    layout = make_layout(h5_patch_grp, lvl_cell_width)
                    add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)

                    if ilvl not in patches:
                        patches[ilvl] = []

                    patches[ilvl].append(Patch(patch_datas))

                    patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

        diag_hier = PatchHierarchy(patch_levels, domain_box, 2, t, data_file)

        return diag_hier

    if load_one_time(time, hier):
        if not silent:
            print("loading data at time {} into existing hierarchy".format(time))
        h5_time_grp = data_file[time]
        t = float(time.strip("t"))

        if t in hier.time_hier:
            if not silent:
                print("time already exist, adding data...")

            # time already exists in the hierarchy
            # all we need to do is adding the data
            # as patchDatas in the appropriate patches
            # and levels, if data compatible with hierarchy

            patch_levels = hier.levels(t)

            for plvl_key in h5_time_grp.keys():
                ilvl = int(plvl_key[2:])
                lvl_cell_width = root_cell_width / 2 ** ilvl

                for ipatch, pkey in enumerate(h5_time_grp[plvl_key].keys()):
                    h5_patch_grp = h5_time_grp[plvl_key][pkey]

                    if patch_has_datasets(h5_patch_grp):
                        hier_patch = patch_levels[ilvl].patches[ipatch]
                        origin = float(h5_time_grp[plvl_key][pkey].attrs['origin'])
                        upper = int(h5_time_grp[plvl_key][pkey].attrs['upper'])
                        lower = int(h5_time_grp[plvl_key][pkey].attrs['lower'])
                        file_patch_box = Box(lower, upper)

                        assert file_patch_box == hier_patch.box
                        assert abs(origin - hier_patch.origin[0]) < 1e-6
                        assert abs(lvl_cell_width - hier_patch.dx) < 1e-6

                        layout = make_layout(h5_patch_grp, lvl_cell_width)
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

            lvl_cell_width = root_cell_width / 2 ** ilvl
            lvl_patches = []

            for ipatch, pkey in enumerate(h5_time_grp[plvl_key].keys()):
                h5_patch_grp = h5_time_grp[plvl_key][pkey]

                if patch_has_datasets(h5_patch_grp):
                    layout = make_layout(h5_patch_grp, lvl_cell_width)
                    patch_datas = {}
                    add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)
                    lvl_patches.append(Patch(patch_datas))

            patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

        hier.time_hier[t] = patch_levels
        return hier

    if load_all_times(time, hier):
        if not silent:
            print("loading all times in existing hier")
        for time in data_file.keys():
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

    root_cell_width = float(simulator.cell_width())
    domain_box = Box(0, int(simulator.domain_box()))

    for ilvl in range(nbr_levels):

        lvl_cell_width = root_cell_width / 2 ** ilvl

        patches = {ilvl : [] for ilvl in range(nbr_levels)}
        getters = quantidic(ilvl, dw)


        if isFieldQty(qty):
            wpatches = getters[qty]()
            for patch in wpatches:
                patch_datas = {}
                lower = int(patch.lower[0])
                upper = int(patch.upper[0])
                origin = float(patch.origin)
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

                lower = int(patch.lower[0])
                upper = int(patch.upper[0])
                origin = float(patch.origin)
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

                        lower = int(dwpatch.lower[0])
                        upper = int(dwpatch.upper[0])

                        box = Box(lower, upper)

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

    print(h5_filename)
    raise ValueError("can't make hierarchy")



def merge_particles(hierarchy):

    for time, patch_levels in hierarchy.time_hier.items():
        for ilvl, plvl in hierarchy.patch_levels.items():
            for ip, patch in enumerate(plvl.patches):
                pdatas = patch.patch_datas
                domain_pdata = [(pdname,pd) for pdname, pd in pdatas.items() if "domain" in pdname][0]

                pghost_pdatas = [(pdname,pd) for pdname, pd in pdatas.items() if "patchGhost" in pdname]
                lghost_pdatas = [(pdname,pd) for pdname, pd in pdatas.items() if "levelGhost" in pdname]

                pghost_pdata = pghost_pdatas[0] if pghost_pdatas else None
                lghost_pdata = lghost_pdatas[0]if lghost_pdatas else None

                if pghost_pdata is not None:
                    domain_pdata[1].dataset.add(pghost_pdata[1].dataset)
                    del pdatas[pghost_pdata[0]]

                if lghost_pdata is not None:
                    domain_pdata[1].dataset.add(lghost_pdata[1].dataset)
                    del pdatas[lghost_pdata[0]]

                popname = domain_pdata[0].split('_')[0]
                pdatas[popname+"_particles"] = pdatas[domain_pdata[0]]
                del pdatas[domain_pdata[0]]
