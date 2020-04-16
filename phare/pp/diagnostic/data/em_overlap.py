from phare.pp.diagnostics import _EMPatchData
from .overlap import Overlap, getOverlaps

from phare.core.gridlayout import yee_element_is_primal


class EMOverlap(Overlap):
    """
    Class containing information for PatchGhost/Domain EM overlap on the same level, made from overlap_calculate_1d

    Parameters:
    -------------------

    lowerPatch : patch with lower origin
    upperPatch : patch with upper origin
    data_name  : Unused for particles
    nGhosts    : number of ghosts associated with the Particles datasets
    sizes      : size of the overlap box
    """

    patch_data_type = _EMPatchData

    def __init__(self, lowerPatch, upperPatch, data_name, nGhosts, sizes):
        Overlap.__init__(self, lowerPatch, upperPatch, data_name, nGhosts, sizes)
        self.lowerPatch = lowerPatch
        self.upperPatch = upperPatch

    def shared_data(self):
        data_name = self.data_name
        is_primal = int(
            yee_element_is_primal(self.lowerPatch.patch_data.quantity_key + data_name)
        )
        gap = self.nGhosts - self.sizes[0]
        size = (self.nGhosts * 2) + is_primal - gap
        return (
            self.lowerPatch.data(data_name)[-size:],
            self.upperPatch.data(data_name)[:size],
        )


def getEMOverlapsFrom(diags):
    return getOverlaps(EMOverlap, diags)
