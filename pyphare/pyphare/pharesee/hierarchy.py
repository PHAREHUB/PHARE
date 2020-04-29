
from ..core import box as boxm
from ..core.box import Box
from ..core.gridlayout import GridLayout
import numpy as np

class PatchData:

    def __init__(self, box, layout, quantity):
        self.quantity = quantity
        self._box = box
        self.origin = layout.origin
        self.layout = layout


class FieldData(PatchData):
    def __init__(self, box, layout, field_name, data):
        super().__init__(box, layout, 'field')


        self.field_name = field_name
        self.dx = layout.dl[0]

        centering = layout.centering["X"][field_name]
        self._ghosts_nbr = layout.nbrGhosts(layout.interp_order, centering)
        self.ghost_box = boxm.grow(box, self._ghosts_nbr)

        if centering == "primal":
            self.size = self.ghost_box.size() + 1
            offset = 0
        else:
            self.size = self.ghost_box.size()
            offset = 0.5*self.dx

        self.x = self.origin[0] - self._ghosts_nbr * self.dx + np.arange(self.size) * self.dx + offset
        self.dataset = data



class ParticleData(PatchData):
    def __init__(self, box, layout, data):
        super().__init__(box, layout, 'particles')
        self.domain_particles = data
        if layout.interp_order == 1:
            self._ghosts_nbr = 1
        elif layout.interp_order == 2 or layout.interp_order == 3:
            self._ghosts_nbr = 2
        else:
            raise RuntimeError("invalid interpolation order")
        self.ghost_box = boxm.grow(box, self._ghosts_nbr)


class Patch:
    """
    A patch represents a hyper-rectangular region of space
    """

    def __init__(self, box, layout, patch_datas):
        self.box = box
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
        for ilvl, lvl in enumerate(self.time_hier[0.].values()):
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

qties = {"EM_B_x": "Bx",
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


def hierarchy_from(h5_filename):
    """
    returns a PatchHierarchy from a PHARE hdf5 diagnostic file
    """
    data_file = h5py.File(h5_filename, "r")


    root_cell_width = float(data_file.attrs["cell_width"])
    domain_box = Box(0, int(data_file.attrs["domain_box"]))

    diag_hier = None

    for time in data_file.keys():
        t = float(time.strip("t"))

        patch_levels = []
        for plvl_key in data_file[time].keys():
            ilvl = int(plvl_key[2:])

            lvl_cell_width = root_cell_width / 2 ** ilvl

            patches = {}
            for pkey in data_file[time][plvl_key].keys():

                patch_datas = {}
                nbrCells = data_file[time][plvl_key][pkey].attrs['nbrCells']
                origin = float(data_file[time][plvl_key][pkey].attrs['origin'])
                upper = int(data_file[time][plvl_key][pkey].attrs['upper'])
                lower = int(data_file[time][plvl_key][pkey].attrs['lower'])

                patch_box = Box(lower, upper)
                layout = GridLayout(patch_box, origin, lvl_cell_width)
                print("reading patch of {} cells origin at {}, patchbox = {}".format(nbrCells,
                                                                                     origin,
                                                                                     patch_box))

                patch_grp = data_file[time][plvl_key][pkey]
                datasets_names = list(patch_grp.keys())

                for dataset_name in datasets_names:
                    dataset = patch_grp[dataset_name]

                    pdata = FieldData(layout, qties[dataset_name], dataset)
                    patch_datas[dataset_name] = pdata

                if ilvl not in patches:
                    patches[ilvl] = []

                patches[ilvl].append(Patch(patch_datas))

            patch_levels.append(PatchLevel(ilvl, patches[ilvl]))

        diag_hier = PatchHierarchy(patch_levels, domain_box, 2, t, data_file)

    return diag_hier




def add_data(filename, hier):
    """
    add the datasets contained in the PHARE hdf5 data into
    an existing PatchHierarchy.

    data in given file should be compatible with the given
    patch hierarchy, i.e. same number of levels, same number of
    patches with same boxes, etc.

    data should not already exist in the hierarchy
    """
    f = h5py.File(filename, "r")


    root_cell_width = float(f.attrs["cell_width"])
    domain_box = Box(0, int(f.attrs["domain_box"]))

    if domain_box !=  hier.domain_box:
        err = "filename domain_box {} incompatible with hierarchy domain_box {}"
        raise ValueError(err.format(domain_box, hier.domain_box))


    for time in f.keys():
        t = float(time.strip("t"))

        if t in hier.time_hier:

            # time already exists in the hierarchy
            # all we need to do is adding the data
            # as patchDatas in the appropriate patches
            # and levels, if data compatible with hierarchy

            lvl_patches = hier.levels(t).patches

            for plvl_key in f[time].keys():
                ilvl = int(plvl_key[2:])

                lvl_cell_width = root_cell_width / 2 ** ilvl


                for ipatch, pkey in enumerate(f[time][plvl_key].keys()):

                    file_patch = f[time][plvl_key][pkey]
                    hier_patch = lvl_patches[ipatch]

                    nbrCells = f[time][plvl_key][pkey].attrs['nbrCells']
                    origin = float(f[time][plvl_key][pkey].attrs['origin'])
                    upper = int(f[time][plvl_key][pkey].attrs['upper'])
                    lower = int(f[time][plvl_key][pkey].attrs['lower'])

                    file_patch_box = Box(lower, upper)

                    assert file_patch_box == hier_patch.box
                    assert abs(origin - hier_patch.origin[0]) < 1e-6
                    assert abs(lvl_cell_width - hier_patch.dx) < 1e-6

                    layout = GridLayout(file_patch_box, origin, lvl_cell_width)

                    print("reading patch of {} cells origin at {}, patchbox = {}".format(nbrCells,
                                                                                         origin,
                                                                                         file_patch_box))

                    for key in file_patch.keys():

                        if key in hier_patch.patch_datas:
                            raise ValueError("dataset {} already in the hierarchy".format(key))

                        dataset = file_patch[key]

                        pdata = FieldData(layout, qties[key], dataset)
                        hier_patch.patch_datas.update({key: pdata})
                        hier.data_files.update({key:f})

        else:

            # time does not exist in the hierarchy
            # we have to create a brand new set of patchLevels
            # containing patches, and load data in their patchdatas

            patch_levels = {}

            for plvl_key in f[time].keys():
                ilvl = int(plvl_key[2:])

                lvl_cell_width = root_cell_width / 2 ** ilvl

                lvl_patches = []

                for ipatch, pkey in enumerate(f[time][plvl_key].keys()):

                    file_patch = f[time][plvl_key][pkey]


                    nbrCells = f[time][plvl_key][pkey].attrs['nbrCells']
                    origin = float(f[time][plvl_key][pkey].attrs['origin'])
                    upper = int(f[time][plvl_key][pkey].attrs['upper'])
                    lower = int(f[time][plvl_key][pkey].attrs['lower'])

                    file_patch_box = Box(lower, upper)

                    layout = GridLayout(file_patch_box, origin, lvl_cell_width)

                    print("reading patch of {} cells origin at {}, patchbox = {}".format(nbrCells,
                                                                                         origin,
                                                                                         file_patch_box))

                    patch_datas = {}
                    for key in file_patch.keys():

                        dataset = file_patch[key]

                        pdata = FieldData(layout, qties[key], dataset)
                        patch_datas.update({key: pdata})
                        #hier.data_files.update({key: f})  # TODO, file handles per time? deal with restarts

                    lvl_patches.append(Patch(patch_datas))

                patch_levels[ilvl] = PatchLevel(ilvl, lvl_patches)

            hier.time_hier[t] = patch_levels
