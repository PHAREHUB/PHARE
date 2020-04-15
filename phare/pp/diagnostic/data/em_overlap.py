from phare.pp.diagnostics import _EMPatchData
from .overlap import Overlap, getOverlaps

from phare.core.gridlayout import yee_element_is_primal


class EMOverlap(Overlap):
    patch_data_type = _EMPatchData

    def __init__(self, patch0, patch1, data_name, nGhosts, sizes):
        Overlap.__init__(self, patch0, patch1, data_name, nGhosts, sizes)

    def shared_data(self):
        data_name = self.data_name
        is_primal = int(
            yee_element_is_primal(self.patch0.patch_data.quantity_key + data_name)
        )
        gap_addition = 2 if self.sizes[0] != self.nGhosts else 0
        size = (self.sizes[0] * 2) + is_primal + gap_addition
        return (
            self.patch0.data(data_name)[-size:],
            self.patch1.data(data_name)[:size],
        )


def getEMOverlapsFrom(diags):
    return getOverlaps(EMOverlap, diags)
