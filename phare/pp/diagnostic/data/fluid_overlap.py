from phare.pp.diagnostics import Diagnostic, Patch, _Fluid
from .periodic_overlap import PeriodicOverlap
from .overlap import Overlap


class FluidOverlap(Overlap):
    dType = _Fluid

    def __init__(self, patch0, patch1, dataset_key, nGhosts, offsets, sizes):
        Overlap.__init__(self, patch0, patch1, dataset_key, nGhosts, offsets, sizes)

    @classmethod
    def get(clazz, diags):
        return Overlap._get(diags, clazz)
