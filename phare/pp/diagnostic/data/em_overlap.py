from phare.pp.diagnostics import _EM
from .overlap import Overlap, getOverlaps


class EMOverlap(Overlap):
    patch_data_type = _EM

    def __init__(self, patch0, patch1, dataset_key, nGhosts, sizes):
        Overlap.__init__(self, patch0, patch1, dataset_key, nGhosts, sizes)

    def get_shared_data(self):
        size = self.sizes[0]
        m_size = int(size * -1)
        return (
            self.p0.patch_data.get()[self.key][m_size:],
            self.p1.patch_data.get()[self.key][:size],
        )


def getEMOverlapsFrom(diags):
    return getOverlaps(EMOverlap, diags)
