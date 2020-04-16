from phare.pp.diagnostics import _FluidPatchData
from .overlap import Overlap, getOverlaps


class FluidOverlap(Overlap):
    """
    Class containing information for PatchGhost/Domain fluid overlap on the same level, made from overlap_calculate_1d

    Parameters:
    -------------------

    lowerPatch : patch with lower origin
    upperPatch : patch with upper origin
    data_name  : dataset key for patch.data(key) pertaining to this overlap
    nGhosts    : number of ghosts considered in the overlap
    sizes      : size of the overlap box
    """

    patch_data_type = _FluidPatchData

    def __init__(self, lowerPatch, upperPatch, data_name, nGhosts, sizes):
        Overlap.__init__(self, lowerPatch, upperPatch, data_name, nGhosts, sizes)
        self.lowerPatch = lowerPatch
        self.upperPatch = upperPatch

    def shared_data(self):
        """ in the case of fluid patch overlaps, only the border value will be equal.
             and for this to occur the patches need to be touching"""
        patches_touch = self.sizes[0] == self.nGhosts
        if patches_touch:
            data_name = self.data_name
            nGhosts = self.nGhosts

            return (
                [self.lowerPatch.data(data_name)[-nGhosts - 1]],
                [self.upperPatch.data(data_name)[nGhosts]],
            )
        return ([], [])


def getFluidOverlapsFrom(diags):
    return getOverlaps(FluidOverlap, diags)
