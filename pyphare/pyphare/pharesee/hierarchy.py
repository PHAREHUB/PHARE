
from ..core import box as boxm
from ..core.box import Box
from ..core.gridlayout import GridLayout
from .particles import Particles
import numpy as np
import os

class PatchData:

    def __init__(self, layout, quantity):
        self.quantity = quantity
        self.box = layout.box
        self.origin = layout.origin
        self.layout = layout


class FieldData(PatchData):
    def __init__(self, layout, field_name, data):
        super().__init__(layout, 'field')


        self.field_name = field_name
        self.dx = layout.dl[0]

        centering = layout.centering["X"][field_name]
        self._ghosts_nbr = layout.nbrGhosts(layout.interp_order, centering)
        self.ghost_box = boxm.grow(layout.box, self._ghosts_nbr)

        if centering == "primal":
            self.size = self.ghost_box.size() + 1
            offset = 0
        else:
            self.size = self.ghost_box.size()
            offset = 0.5*self.dx

        self.x = self.origin[0] - self._ghosts_nbr * self.dx + np.arange(self.size) * self.dx + offset
        self.dataset = data



class ParticleData(PatchData):
    def __init__(self, layout, data):
        super().__init__(layout, 'particles')
        self.dataset = data
        if layout.interp_order == 1:
            self._ghosts_nbr = 1
        elif layout.interp_order == 2 or layout.interp_order == 3:
            self._ghosts_nbr = 2
        else:
            raise RuntimeError("invalid interpolation order")
        self.ghost_box = boxm.grow(layout.box, self._ghosts_nbr)


class Patch:
    """
    A patch represents a hyper-rectangular region of space
    """

    def __init__(self, patch_datas):
        # we assume patch data all have the same layout
        keys = list(patch_datas.keys())
        patchData0 = patch_datas[keys[0]]
        layout = patchData0.layout
        self.box = layout.box
        self.origin = layout.origin
        self.dx = layout.dl[0]
        self.patch_datas = patch_datas


class PatchLevel:
    """is a collection of patches """

    def __init__(self, lvl_nbr, patches):
        self.level_number = lvl_nbr
        self.patches = patches


class PatchHierarchy:
    """is a collection of patch levels """


    def __init__(self, patch_levels, domain_box, refinement_ratio, time=0., data_files=None):
        self.patch_levels = patch_levels
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
        return boxm.refine(self.domain_box, self.refinement_ratio ** level_number)

    def __str__(self):
        s = "Hierarchy: \n"
        for t, patch_levels in self.time_hier.items():
            print("t = {}".format(t))
            for ilvl, lvl in enumerate(patch_levels):
                s = s + "Level {}\n".format(ilvl)
                for ip, patch in enumerate(lvl.patches):
                    for qty_name, pd in patch.patch_datas.items():
                        pdstr = "    P{ip} {pdname} box is {box} and ghost box is {gbox}"
                        s = s + pdstr.format(ip=ip, pdname=qty_name,
                                             box=patch.box, gbox=pd.ghost_box)
                        s = s + "\n"
        return s


    def plot(self):

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 3))
        for ilvl, lvl in enumerate(self.time_hier[0.]):
            lvl_offset = ilvl * 0.1
            for patch in lvl.patches:
                dx = patch.dx
                origin = patch.origin
                x0 = patch.box.lower * dx
                x1 = patch.box.upper * dx
                xcells = np.arange(x0, x1 + dx, dx)
                y = lvl_offset + np.zeros_like(xcells)
                ax.plot(xcells, y, marker=".")

        fig.savefig("hierarchy.png")



def is_root_lvl(patch_level):
    return patch_level.level_number == 0






import h5py

field_qties = {"EM_B_x": "Bx",
               "EM_B_z": "By",
               "EM_B_y": "Bz",
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
    return "alpha_domain" from ions_pop_ions_alpha_domain.h5
    """
    popname = basename.strip(".h5").split("_")[-2]
    part_type = basename.strip(".h5").split("_")[-1]
    dataset_name = popname + "_" + part_type

    return dataset_name


def is_pop_fluid_file(basename):
    return (is_particle_file(basename) is False) and "pop" in basename


def pop_name(basename):
    return basename.strip(".h5").split("_")[-2]


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




def add_to_patchdata(patch_datas, h5_patch_grp, basename, layout):
    """
    adds data in the h5_patch_grp in the given PatchData dict
    """

    if is_particle_file(basename):
        particles = Particles(icells=h5_patch_grp["iCell"],
                              deltas=h5_patch_grp["delta"],
                              v=h5_patch_grp["v"],
                              weights=h5_patch_grp["weight"],
                              charges=h5_patch_grp["charge"])

        pdname = particle_dataset_name(basename)
        if pdname in patch_datas:
            raise ValueError("error - {} already in patchdata".format(pdname))

        patch_datas[pdname] = ParticleData(layout, particles)

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






def hierarchy_from(h5_filename, time=None, hier=None):
    """
    this function reads an HDF5 PHARE file and returns a PatchHierarchy from
    which data is accessible.
    if 'time' is None, all times in the file will be read, if a time is given
    then only that time will be read
    if 'hier' is None, then a new hierarchy will be created, if not then the
    given hierarchy 'hier' will be filled.

    The function fails if the data is already in hierarchy
    """

    data_file = h5py.File(h5_filename, "r")
    basename = os.path.basename(h5_filename)
    root_cell_width = float(data_file.attrs["cell_width"])
    domain_box = Box(0, int(data_file.attrs["domain_box"]))


    if create_from_all_times(time, hier):
        # first create from first time
        # then add all other times
        print("creating hierarchy from all times in file")
        times = list(data_file.keys())
        hier = hierarchy_from(h5_filename, time=times[0])
        if len(times)>1:
            for time in times[1:]:
                hierarchy_from(h5_filename, time=time, hier=hier)
        return hier



    if create_from_one_time(time, hier):
        print("creating hierarchy from one time {}".format(time))
        t = float(time.strip("t"))

        h5_time_grp = data_file[time]
        patch_levels = []

        for plvl_key in h5_time_grp.keys():

            h5_patch_lvl_grp = h5_time_grp[plvl_key]
            ilvl = int(plvl_key[2:])
            lvl_cell_width = root_cell_width / 2 ** ilvl
            patches = {}

            for pkey in h5_patch_lvl_grp.keys():

                h5_patch_grp = data_file[time][plvl_key][pkey]
                patch_datas = {}
                layout = make_layout(h5_patch_grp, lvl_cell_width)

                add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)

                if ilvl not in patches:
                    patches[ilvl] = []

                patches[ilvl].append(Patch(patch_datas))

            patch_levels.append(PatchLevel(ilvl, patches[ilvl]))

        diag_hier = PatchHierarchy(patch_levels, domain_box, 2, t, data_file)

        return diag_hier



    if load_one_time(time, hier):
        print("loading data at time {} into existing hierarchy".format(time))
        h5_time_grp = data_file[time]
        t = float(time.strip("t"))

        if t in hier.time_hier:
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
                    #hier_patch.patch_datas.update({dataset_name: pdata})

                    #hier.data_files.update({dataset_name: f})

        else:
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

                    layout = make_layout(h5_patch_grp, lvl_cell_width)

                    patch_datas = {}

                    add_to_patchdata(patch_datas, h5_patch_grp, basename, layout)

                    lvl_patches.append(Patch(patch_datas))

                patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

            hier.time_hier[t] = patch_levels
            return hier


    if load_all_times(time, hier):
        print("loading all times in existing hier")
        for time in data_file.keys():
            hier = hierarchy_from(h5_filename, time=time, hier=hier)

        return hier







