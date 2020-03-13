from phare.pp.diagnostics import _Fluid
from .overlap import Overlap, getOverlaps


class FluidOverlap(Overlap):
    dType = _Fluid

    def __init__(self, patch0, patch1, dataset_key, nGhosts, sizes):
        Overlap.__init__(self, patch0, patch1, dataset_key, nGhosts, sizes)

    def get_shared_data(self):
        """ in the case of fluid patch overlaps, only the border value will be equal.
             and for this to occur the patches need to be touching"""
        if self.sizes[0] == self.nGhosts:
            key, nGhosts = self.key, self.nGhosts
            m_nGhosts = int(nGhosts * -1)
            return (
                [self.p0.dtype.get()[key][m_nGhosts - 1]],
                [self.p1.dtype.get()[key][nGhosts]],
            )
        return ([], [])


def getFluidOverlapsFrom(diags):
    return getOverlaps(FluidOverlap, diags)
