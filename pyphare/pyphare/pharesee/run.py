import os
from .hierarchy import hierarchy_from, finest_field
import numpy as np
from pyphare.pharesee.hierarchy import compute_hier_from


def _current1d(by, bz, xby, xbz):
    # jx = 0
    # jy = -dxBz
    # jz = dxBy
    # the following hard-codes yee layout
    # which is not true in general
    # we should at some point provide proper
    # derivation routines in the gridlayout
    dx = xbz[1]-xbz[0]
    jy = np.zeros(by.size+1)
    jy[1:-1] = -(bz[1:]-bz[:-1])/dx
    dx = xby[1]-xby[0]
    jz = np.zeros(bz.size+1)
    jz[1:-1] = (by[1:]-by[:-1])/dx
    jy[0]=jy[1]
    jy[-1]=jy[-2]
    jz[0]=jz[1]
    jz[-1]=jz[-2]
    return jy, jz


def _compute_current(patch):
    By = patch.patch_datas["By"].dataset[:]
    xby  = patch.patch_datas["By"].x
    Bz = patch.patch_datas["Bz"].dataset[:]
    xbz  = patch.patch_datas["Bz"].x
    Jy, Jz =  _current1d(By, Bz, xby, xbz)
    return ({"name":"Jy", "data":Jy,"centering":"primal"},
            {"name":"Jz", "data":Jz,"centering":"primal"})


class Run:
    def __init__(self, path):
        self.path = path

    def _get_hierarchy(self, time, filename, hier=None):
        t = "{:.10f}".format(time)
        return hierarchy_from(h5_filename=os.path.join(self.path, filename), time=t, hier=hier)

    def _get(self, hierarchy, time, merged):
        if merged:
            merged_qties = {}
            for qty in hierarchy.quantities():
                merged_qties[qty] = finest_field(hierarchy, qty, time=time)
            return merged_qties
        else:
            return hierarchy

    def GetB(self, time, merged=False):
        hier = self._get_hierarchy(time, "EM_B.h5")
        return self._get(hier, time, merged)

    def GetE(self, time, merged=False):
        hier = self._get_hierarchy(time, "EM_E.h5")
        return self._get(hier, time, merged)

    def GetNi(self, time, merged=False):
        hier = self._get_hierarchy(time, "ions_density.h5")
        return self._get(hier, time, merged)

    def GetN(self, time, pop_name, merged=False):
        hier =  self._get_hierarchy(time, "ions_{}_density.h5".format(pop_name))
        return self._get(hier, time, merged)

    def GetVi(self, time, merged=False):
        hier =  self._get_hierarchy(time, "ions_bulkVelocity.h5")
        return self._get(hier, time, merged)

    def GetFlux(self, time, pop_name, merged=False):
        hier = self._get_hierarchy(time, "ions_pop_{}_flux.h5".format(pop_name))
        return self._get(hier, time, merged)

    def GetJ(self, time, merged=False):
        B = self.GetB(time)
        J = compute_hier_from(B, _compute_current)
        return self._get(J, time, merged)

    def GetParticles(self, time, pop_name, hier=None):
        def filename(name):
            return f"ions_pop_{name}_domain.h5"
        if isinstance(pop_name, (list, tuple)):
            for pop in pop_name:
                hier = self._get_hierarchy(time, filename(pop), hier=hier)
            return hier
        return self._get_hierarchy(time, filename(pop_name), hier=hier)
