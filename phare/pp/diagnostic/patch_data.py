import h5py
from abc import ABC  # abstract class


class PatchData(ABC):
    """
    Abstract Class for fowarding different types of quantity data

    Parameters:
    -------------------

    diag         : Diagnostic object from which this PatchData is contained
    quantity_key : string representing the physical quantity eg. E/B/density
    h5item       : Dataset for single data, Group for multiple, immediate children being data

    """

    def __init__(self):
        raise ValueError("Abstract class not to be instantiated manually")

    def data(self, str = ""):
        raise ValueError("Abstract class not to be instantiated manually")

    def __init__(self, diag, quantity_key, h5item):
        assert isinstance(h5item, h5py.Dataset) or isinstance(h5item, h5py.Group)

        self.diag = diag
        self.quantity_key = quantity_key
        self.h5item = h5item  # group or dataset
        self.name = (
            self.h5item.parent.name
            if isinstance(self.h5item, h5py.Dataset)
            else self.h5item.name
        ).split("/")[-1] # take last entry in path of hdf5 item

    def items(self):
        return self.data().items()

    def nGhosts(self, dataset_key=""):
        if isinstance(self.h5item, h5py.Dataset):  # Fields
            return int(self.data()[dataset_key].attrs["ghosts"])
        if dataset_key in list(self.h5item.keys()):  # VecFields
            return int(self.data()[dataset_key].attrs["ghosts"])
        raise ValueError(
            "nGhosts can not be discerned, no ghosts attributes in hdf5 file"
        )

    def dataset_names(self):
        if isinstance(self.h5item, h5py.Dataset):
            return [self.h5item.name.split("/")[-1]]
        return list(self.h5item.keys())

    def _get_VF(self, group):
        return {v: group[v] for v in list(group.keys())}


class _EMPatchData(PatchData):
    def __init__(self, diag, quantity_key, h5item):
        quantity_key = quantity_key.split("_")[1]  # drop EM_ from EM_B
        PatchData.__init__(self, diag, quantity_key, h5item)

    def data(self, ds_name=""):
        datasets = self._get_VF(self.h5item)
        if len(ds_name):
            return datasets[ds_name]
        return datasets


class _FluidPatchData(PatchData):
    def __init__(self, diag, quantity_key, h5item):
        PatchData.__init__(self, diag, quantity_key, h5item)
        if "pop/" in h5item.name:
            self.pop_name = self.h5item.name.split("/")[-2]

    def data(self, ds_name=""):
        datasets = None
        if self.quantity_key == "density":
            datasets = {"density": self.h5item}
        else:
            datasets = self._get_VF(self.h5item)
        if len(ds_name):
            return datasets[ds_name]
        return datasets



class _ParticlePatchData(PatchData):
    dataset_keys = ["weight", "charge", "iCell", "delta", "v"]

    def __init__(self, diag, quantity_key, h5item):
        PatchData.__init__(self, diag, quantity_key, h5item)
        self.pop_name = self.h5item.name.split("/")[-2]

    def data(self, ds_name=""):
        datasets = {v: self.h5item[v] for v in _ParticlePatchData.dataset_keys}
        if len(ds_name):
            return datasets[ds_name]
        return datasets

    def nGhosts(self, dataset_key=""):
        if "ghosts" in self.h5item.attrs:
            return int(self.h5item.attrs["ghosts"])
        raise ValueError(
            "nGhosts can not be discerned, no ghosts attributes in hdf5 file"
        )


def particlesForDiags(diags):
    """Takes per HDF5 file Diagnostic objects and coallates the particle population
        types into one replica Diagnostic object."""

    particles, particle_pops = [], {}
    for diag in diags:
        if diag.patch_data_type is _ParticlePatchData:
            assert 1 in diag.levels  # levelGhost has no level0
            patches = list(diag.levels[1].patches.values())
            pop_name = patches[0].patch_data.pop_name
            name = patches[0].patch_data.name
            if pop_name not in particle_pops:
                particle_pops[pop_name] = {}
            particle_pops[pop_name][name] = diag
    for pop_name, diags in particle_pops.items():
        particles.append(Particles(pop_name, **diags))
    return particles


class Particles:
    """Wrapper container for domain, patchGhost and levelGhost for the same population
        Is used in place of phare.pp.Diagnostic so may need to forward attributes"""

    def __init__(self, pop_name, **diags):
        assert all([key in diags for key in ["domain", "patchGhost", "levelGhost"]])

        self.pop_name = pop_name
        self.domainDiag = diags["domain"]
        self.pGhostDiag = diags["patchGhost"]
        self.lGhostDiag = diags["levelGhost"]
        self.dim = self.domainDiag.dim
        self.domain_upper = self.domainDiag.domain_upper
        self.origin = self.domainDiag.origin

    def min_coord(self, direction):
        return self.domainDiag.min_coord(direction)

    def max_coord(self, direction):
        return self.domainDiag.max_coord(direction)
