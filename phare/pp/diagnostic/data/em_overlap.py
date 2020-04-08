from phare.pp.diagnostics import _EMPatchData
from .overlap import Overlap, getOverlaps

from phare.data.wrangler import DataWrangler, is_primal_dataset


class EMOverlap(Overlap):
    patch_data_type = _EMPatchData

    def __init__(self, patch0, patch1, data_name, nGhosts, sizes):
        Overlap.__init__(self, patch0, patch1, data_name, nGhosts, sizes)

    def shared_data(self):
        data_name = self.data_name
        size = self.sizes[0]
        is_primal = int(
            is_primal_dataset(self.patch0.patch_data.name + "_" + data_name)
        )
        return (
            self.patch0.patch_data.data(data_name)[-size * 2 - is_primal :],
            self.patch1.patch_data.data(data_name)[: size * 2 + is_primal],
        )


def getEMOverlapsFrom(diags):
    return getOverlaps(EMOverlap, diags)
