import os
from .hierarchy import hierarchy_from
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
    By = patch.patch_datas["EM_B_y"].dataset[:]
    xby  = patch.patch_datas["EM_B_y"].x
    Bz = patch.patch_datas["EM_B_z"].dataset[:]
    xbz  = patch.patch_datas["EM_B_z"].x
    Jy, Jz =  _current1d(By, Bz, xby, xbz)
    return ({"name":"J_y", "data":Jy,"centering":"primal"},
            {"name":"J_z", "data":Jz,"centering":"primal"})


class Run:
    def __init__(self, path):
        self.path = path

    def _get_hierarchy(self, time, filename, hier=None):
        t = "{:.10f}".format(time)
        return hierarchy_from(h5_filename=os.path.join(self.path, filename), time=t, hier=hier)

    def GetB(self, time):
        return self._get_hierarchy(time, "EM_B.h5")

    def GetE(self, time):
        return self._get_hierarchy(time, "EM_E.h5")

    def GetNi(self, time):
        return self._get_hierarchy(time, "ions_density.h5")

    def GetN(self, time, pop_name):
        return self._get_hierarchy(time, "ions_{}_density.h5".format(pop_name))

    def GetVi(self, time):
        return self._get_hierarchy(time, "ions_bulkVelocity.h5")

    def GetFlux(self, time, pop_name):
        return self._get_hierarchy(time, "ions_pop_{}_flux.h5".format(pop_name))

    def GetJ(self, time):
        B = self.GetB(time)
        return compute_hier_from(B, _compute_current)

    def GetParticles(self, time, pop_name, hier=None):
        def filename(name):
            return f"ions_pop_{name}_domain.h5"
        if isinstance(pop_name, (list, tuple)):
            for pop in pop_name:
                hier = self._get_hierarchy(time, filename(pop), hier=hier)
            return hier
        return self._get_hierarchy(time, filename(pop_name), hier=hier)
