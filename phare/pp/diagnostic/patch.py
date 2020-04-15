REFINEMENT_RATIO = 2
xyz_to_dim = {"x": 0, "y": 1, "z": 2}


def max_amr_index_for(patch_level, direction):
    assert isinstance(patch_level, PatchLevel) and direction in xyz_to_dim

    return (
        REFINEMENT_RATIO ** patch_level.lvlNbr
        * patch_level.diag.cells[xyz_to_dim[direction]]
    )


class PatchLevel:

    """
    Class representing AMR hierarchy PatchLevel for the purposes of forwarding data contained with HDF5 "Patches"

    Parameters:
    -------------------

    diag                 : Diagnostic object from which this PatchLevel is contained
    lvlNbn               : int - level index in the Hierarchy, starts at 0, < Hierarchy.levels.size()

    """

    def __init__(self, diag, lvlNbr):
        self.diag = diag
        self.lvlNbr = lvlNbr
        power = REFINEMENT_RATIO ** lvlNbr
        self.dl = [x / power for x in diag.dl]
        self.patches = {}

    def is_root(self):
        return self.lvlNbr == 0

    def position_to_amr_ndex(self, pos, direction):
        origin = self.diag.origin[xyz_to_dim[direction]]
        return round((pos - origin) / self.cell_width(direction))

    def patch_length(self, cells, direction):
        return self.cell_width(direction) * cells

    def cell_width(self, direction):
        return self.dl[xyz_to_dim[direction]]

    def data_names(self):
        return self._first_patch().patch_data.dataset_names()

    def nGhosts(self, data_name):
        return self._first_patch().patch_data.nGhosts(data_name)

    def _first_patch(self):
        patches = self.patches_list()
        assert len(patches)
        return patches[0]

    def patches_list(self):
        return list(self.patches.values())


class Patch:
    """
    Class representing AMR hierarchy Patch

    Parameters:
    -------------------

    patch_level    : PatchLevel object within which this Patch exists
    h5PatchGroup   : h5py object for the group for this Patch in the HDF5 file
    patch_data     : subclass(PatchData) with access to quantity data for this patch

    """

    def __init__(self, patch_level, h5PatchGroup, patch_data):
        self.patch_level = patch_level
        self.h5PatchGroup = h5PatchGroup
        self.id = h5PatchGroup.name.split("/")[-1][1:]  # samrai patch id e.g. 0x0
        self.patch_data = patch_data
        self.origin = [float(v) for v in h5PatchGroup.attrs["origin"].split(",")]
        self.cells = [int(v) for v in h5PatchGroup.attrs["nbrCells"].split(",")]

    # copy/transform patch origins for periodic overlap calculations
    def copy(self, transform=[0, 0, 0]):
        p = Patch(self.patch_level, self.h5PatchGroup, self.patch_data)
        p.origin = [f + transform[i] for i, f in enumerate(self.origin)]
        return p

    def min_coord(self, direction):
        return self.origin[xyz_to_dim[direction]]

    def max_coord(self, direction):
        return round(
            self.min_coord(direction)
            + self.patch_level.patch_length(
                self.cells[xyz_to_dim[direction]], direction
            ),
            6,
        )

    def min_max_coords(self, direction):
        return (self.min_coord(direction), self.max_coord(direction))

    def data(self, ds_name=""):
        return self.patch_data.data(ds_name)


def aggregate_level0_patch_domain(patch_level, shared_patch_border=False, ds_names=[]):
    """ ds_names is an optional list to aggregate only matching strings
          used for EM, to avoid mixing non/primal ds_names
    """
    assert isinstance(shared_patch_border, bool)  # you never know
    assert patch_level.lvlNbr == 0 and patch_level.diag.dim == 1
    import numpy as np

    patches = sorted(patch_level.patches.values(), key=lambda x: x.min_coord("x"))

    physical_datasets = {}
    for patch in patches:
        for ds_name in patch.patch_data.dataset_names():
            if len(ds_names) and not ds_name in ds_names:
                continue
            nGhosts = patch.patch_data.nGhosts(ds_name)
            end = -nGhosts
            if shared_patch_border:
                # first value of second patch dataset is last of first
                end = end if patch == patches[-1] else end - 1
            if ds_name not in physical_datasets:
                physical_datasets[ds_name] = []

            hdf5_data = patch.data(ds_name)
            physical_datasets[ds_name].append(np.asarray(hdf5_data[nGhosts:end]))

    assert len(
        list(physical_datasets.keys())
    )  # you should expect one if you use this function

    return {
        dataset_key: np.hstack(array_list)
        for dataset_key, array_list in physical_datasets.items()
    }
