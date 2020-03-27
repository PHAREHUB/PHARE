import h5py


class PatchData:
    def __init__(self, diag, key, h5):
        self.diag = diag
        self.key = key
        self.h5 = h5
        self.name = (
            self.h5.parent.name if isinstance(self.h5, h5py.Dataset) else self.h5.name
        )
        self.name = self.name.split("/")[-1]
        diag.type = self

    def items(self):
        return self.get().items()

    def nGhosts(self, dataset_key=""):
        if isinstance(self.h5, h5py.Dataset):  # Fields
            return self.get()[dataset_key].attrs["ghosts"]
        if dataset_key in list(self.h5.keys()):  # VecFields
            return self.get()[dataset_key].attrs["ghosts"]
        raise ValueError(
            "nGhosts can not be discerned, no ghosts attributes in hdf5 file"
        )

    def keys(self):
        if isinstance(self.h5, h5py.Dataset):
            return [self.h5.name.split("/")[-1]]
        return list(self.h5.keys())

    def _get_VF(self, group):
        return {v: group[v] for v in list(group.keys())}


class _EM(PatchData):
    def __init__(self, diag, key, h5):
        key = key.split("_")[1]  # drop EM_ from EM_B
        PatchData.__init__(self, diag, key, h5)

    def get(self):
        return self._get_VF(self.h5)


class _Fluid(PatchData):
    def __init__(self, diag, key, h5):
        PatchData.__init__(self, diag, key, h5)

    def get(self):
        if self.key == "density":
            return {"density": self.h5}
        return self._get_VF(self.h5)


class _Particle(PatchData):
    datasets = ["weight", "charge", "iCell", "delta", "v"]

    def __init__(self, diag, val, h5):
        PatchData.__init__(self, diag, val, h5)
        self.pop = self.h5.name.split("/")[-2]

    def get(self):
        return {v: self.h5[v] for v in _Particle.datasets}

    def nGhosts(self, dataset_key=""):
        if "ghosts" in self.h5.attrs:
            return self.h5.attrs["ghosts"]
        raise ValueError(
            "nGhosts can not be discerned, no ghosts attributes in hdf5 file"
        )


class Particles:
    """Wrapper container for domain, patchGhost and levelGhost for the same population
        Is used in place of phare.pp.Diagnostic so may need to forward attributes"""

    def __init__(self, pop, domain, pGhost, lGhost):
        self.pop = pop
        self.domain = domain
        self.pGhost = pGhost
        self.lGhost = lGhost
        self.sim = domain.sim
        self.dim = len(domain.dl)
        self.type = domain.type

    def min_coord(self, direction):
        return self.domain.min_coord(direction)

    def max_coord(self, direction):
        return self.domain.max_coord(direction)
