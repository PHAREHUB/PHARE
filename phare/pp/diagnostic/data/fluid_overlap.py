from phare.pp.diagnostics import _FluidPatchData
from .overlap import Overlap, getOverlaps


class FluidOverlap(Overlap):
    patch_data_type = _FluidPatchData

    def __init__(self, patch0, patch1, data_name, nGhosts, sizes):
        Overlap.__init__(self, patch0, patch1, data_name, nGhosts, sizes)

    def shared_data(self):
        """ in the case of fluid patch overlaps, only the border value will be equal.
             and for this to occur the patches need to be touching"""
        patches_touch = self.sizes[0] == self.nGhosts
        if patches_touch:
            data_name = self.data_name
            nGhosts = self.nGhosts

            return (
                [self.patch0.data(data_name)[-nGhosts - 1]],
                [self.patch1.data(data_name)[nGhosts]],
            )
        return ([], [])


def getFluidOverlapsFrom(diags):
    return getOverlaps(FluidOverlap, diags)
