import os

import numpy as np
from pyphare.pharesee.hierarchy import compute_hier_from

from .hierarchy import flat_finest_field, hierarchy_from


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

def _current2d(bx, by, bz, dx, dy):
    # jx = dyBz
    # jy = -dxBz
    # jz = dxBy - dyBx
    # the following hard-codes yee layout
    # which is not true in general
    # we should at some point provide proper
    # derivation routines in the gridlayout
    jx = np.zeros(by.shape)
    jy = np.zeros(bx.shape)
    jz = np.zeros((bx.shape[0], by.shape[1]))

    jx[:,1:-1] = (bz[:,1:] - bz[:,:-1])/dy
    jy[1:-1,:] = -(bz[1:,:] - bz[:-1,:])/dx
    jz[1:-1,1:-1] = (by[1:, 1:-1] - by[:-1 , 1:-1])/dx - (bx[1:-1 , 1:] - bx[1:-1, :-1])/dy

    jy[0,:]  = jy[1,:]
    jy[:,0]  = jy[:,1]
    jy[-1,:] = jy[-2,:]
    jy[:,-1] = jy[:,-2]

    jz[0,:]  = jz[1,:]
    jz[:,0]  = jz[:,1]
    jz[-1,:] = jz[-2,:]
    jz[:,-1] = jz[:,-2]

    jx[0,:]  = jx[1,:]
    jx[:,0]  = jx[:,1]
    jx[-1,:] = jx[-2,:]
    jx[:,-1] = jx[:,-2]

    return jx, jy, jz



def _compute_current(patch):
    if patch.box.ndim ==1 : 
        By = patch.patch_datas["By"].dataset[:]
        xby  = patch.patch_datas["By"].x
        Bz = patch.patch_datas["Bz"].dataset[:]
        xbz  = patch.patch_datas["Bz"].x
        Jy, Jz =  _current1d(By, Bz, xby, xbz)
        return ({"name":"Jy", "data":Jy,"centering":"primal"},
                {"name":"Jz", "data":Jz,"centering":"primal"})

    elif patch.box.ndim ==2 :
        Bx = patch.patch_datas["Bx"].dataset[:]
        By = patch.patch_datas["By"].dataset[:]
        Bz = patch.patch_datas["Bz"].dataset[:]

        dx,dy  = patch.dl

        Jx, Jy, Jz =  _current2d(Bx, By, Bz, dx, dy)
        #return ({"name":"Jx", "data":Jx,"centering":"primal"},
        #        {"name":"Jy", "data":Jy,"centering":"primal"},
        #        {"name":"Jz", "data":Jz,"centering":"primal"})


        components = ("Jx", "Jy", "Jz")
        centering = {component:[patch.layout.centering[direction][component] for direction in ("X", "Y")] for component in components}

        return ({"name":"Jx", "data":Jx, "centering":centering["Jx"]},
                {"name":"Jy", "data":Jy, "centering":centering["Jy"]},
                {"name":"Jz", "data":Jz,"centering":centering["Jz"]})


def _divB2D(Bx,By, xBx, yBy):
    dxbx = (Bx[1:,:]-Bx[:-1,:])/(xBx[1]-xBx[0])
    dyby = (By[:,1:]-By[:,:-1])/(yBy[1]-yBy[0])
    return dxbx + dyby

def _compute_divB(patch):

    if patch.box.ndim ==1:
        raise ValueError("divB is 0 by construction in 1D")

    elif patch.box.ndim == 2:
        By  = patch.patch_datas["By"].dataset[:]
        Bx  = patch.patch_datas["Bx"].dataset[:]
        xBx = patch.patch_datas["Bx"].x
        yBy = patch.patch_datas["By"].y
        divB = _divB2D(Bx,By, xBx, yBy)

        return ({"name":"divB", "data":divB, "centering":["dual","dual"]},)

    else:
        raise RuntimeError("dimension not implemented")



def _get_rank(patch):
    """
    make a field dataset cell centered coding the MPI rank
    rank is obtained from patch global id == "rank#local_patch_id"
    """
    from pyphare.core.box import grow
    layout = patch.layout
    centering = "dual"
    nbrGhosts = layout.nbrGhosts(layout.interp_order, centering)
    shape = grow(patch.box, [nbrGhosts]*2).shape

    if patch.box.ndim == 1:
        pass

    if patch.box.ndim == 2:
        data = np.zeros(shape) + int(patch.id.strip("p").split("#")[0])
        return ({"name":"rank", "data":data, "centering":[centering]*2},)
    else:
        raise RuntimeError("Not Implemented yet")



def make_interpolator(data, coords, interp, domain, dl, qty, nbrGhosts):
    """
    :param data: the values of the data that will be used for making
    the interpolator, defined on coords
    :param coords: coordinates where the data are known. they
    can be define on an irregular grid (eg the finest)

    finest_coords will be the structured coordinates defined on the
    finest grid.
    """
    from pyphare.core.gridlayout import yeeCoordsFor
    dim = coords.ndim

    if dim == 1:
        from scipy.interpolate import interp1d

        interpolator = interp1d(coords, data,\
                                kind=interp,\
                                fill_value="extrapolate",\
                                assume_sorted=False)

        nx = 1+int(domain[0]/dl[0])

        x = yeeCoordsFor([0]*dim, nbrGhosts, dl, [nx],  qty, "x")
        finest_coords = (x,)

    elif dim == 2:
        from scipy.interpolate import (LinearNDInterpolator,
                                       NearestNDInterpolator)

        if interp == 'nearest':
            interpolator = NearestNDInterpolator(coords, data)
        elif interp == 'bilinear':
            interpolator = LinearNDInterpolator(coords, data)
        else:
            raise ValueError("interp can only be 'nearest' or 'bilinear'")

        nCells = [1+int(d/dl) for d,dl in zip(domain, dl)]
        x = yeeCoordsFor([0]*dim, nbrGhosts, dl, nCells,  qty, "x")
        y = yeeCoordsFor([0]*dim, nbrGhosts, dl, nCells,  qty, "y")
        #x = np.arange(0, domain[0]+dl[0], dl[0])
        #y = np.arange(0, domain[1]+dl[1], dl[1])
        finest_coords = (x, y)

    else:
        raise ValueError("make_interpolator is not yet 3d")

    return interpolator, finest_coords


class Run:
    def __init__(self, path):
        self.path = path

    def _get_hierarchy(self, time, filename, hier=None):
        t = "{:.10f}".format(time)
        return hierarchy_from(h5_filename=os.path.join(self.path, filename),\
                              time=t, hier=hier)

    def _get(self, hierarchy, time, merged, interp):
        """
        if merged=True, will return an interpolator and a tuple of 1d arrays
        with the coordinates of the finest grid where the interpolator
        can be calculated (that is the return of flat_finest_field)
        """
        if merged:
            domain = self.GetDomainSize()
            dl = self.GetDl()

            # assumes all qties in the hierarchy have the same ghost width
            # so take the first patch data of the first patch of the first level....
            nbrGhosts = list(hierarchy.level(0).patches[0].patch_datas.values())[0].ghosts_nbr
            merged_qties = {}
            for qty in hierarchy.quantities():
                data, coords = flat_finest_field(hierarchy, qty, time=time)
                merged_qties[qty] = make_interpolator(data, coords,\
                                                      interp, domain, dl, qty, nbrGhosts)
            return merged_qties
        else:
            return hierarchy

    def GetTags(self, time, merged=False):
        hier = self._get_hierarchy(time, "tags.h5")
        return self._get(hier, time, merged, "nearest")

    def GetB(self, time, merged=False, interp='nearest'):
        hier = self._get_hierarchy(time, "EM_B.h5")
        return self._get(hier, time, merged, interp)

    def GetE(self, time, merged=False, interp='nearest'):
        hier = self._get_hierarchy(time, "EM_E.h5")
        return self._get(hier, time, merged, interp)

    def GetNi(self, time, merged=False, interp='nearest'):
        hier = self._get_hierarchy(time, "ions_density.h5")
        return self._get(hier, time, merged, interp)

    def GetN(self, time, pop_name, merged=False, interp='nearest'):
        hier =  self._get_hierarchy(time, "ions_pop_{}_density.h5".format(pop_name))
        return self._get(hier, time, merged, interp)

    def GetVi(self, time, merged=False, interp='nearest'):
        hier =  self._get_hierarchy(time, "ions_bulkVelocity.h5")
        return self._get(hier, time, merged, interp)

    def GetFlux(self, time, pop_name, merged=False, interp='nearest'):
        hier = self._get_hierarchy(time, "ions_pop_{}_flux.h5".format(pop_name))
        return self._get(hier, time, merged, interp)

    def GetJ(self, time, merged=False, interp='nearest'):
        B = self.GetB(time)
        J = compute_hier_from(B, _compute_current)
        return self._get(J, time, merged, interp)

    def GetDivB(self, time, merged=False, interp='nearest'):
        B = self.GetB(time)
        db = compute_hier_from(B, _compute_divB)
        return self._get(db, time, merged, interp)

    def GetRanks(self, time, merged=False, interp='nearest'):
        B = self.GetB(time)
        ranks = compute_hier_from(B, _get_rank)
        return self._get(ranks, time, merged, interp)

    def GetParticles(self, time, pop_name, hier=None):
        def filename(name):
            return f"ions_pop_{name}_domain.h5"
        if isinstance(pop_name, (list, tuple)):
            for pop in pop_name:
                hier = self._get_hierarchy(time, filename(pop), hier=hier)
            return hier
        return self._get_hierarchy(time, filename(pop_name), hier=hier)

    def GetMass(self, pop_name):
        list_of_qty = ['density', 'flux', 'domain', 'levelGhost', 'patchGhost']
        list_of_mass = []

        import h5py
        for qty in list_of_qty:
            file = os.path.join(self.path, "ions_pop_{}_{}.h5".format(pop_name, qty))
            if os.path.isfile(file):
                h5_file = h5py.File(file, "r")
                list_of_mass.append(h5_file.attrs["pop_mass"])

        assert all(m==list_of_mass[0] for m in list_of_mass)

        return list_of_mass[0]

    def GetDomainSize(self):
        h5_filename = "EM_B.h5" # _____ TODO : could be another file

        import h5py
        data_file = h5py.File(os.path.join(self.path, h5_filename), "r")

        root_cell_width = np.asarray(data_file.attrs["cell_width"])

        return (data_file.attrs["domain_box"]+1)*root_cell_width

    def GetDl(self, level='finest', time=None):
        """
        gives the ndarray containing the grid sizes at the given time
        for the hierarchy defined in the given run, and for the given level
        (default is 'finest', but can also be a int)

        :param level: the level at which get the associated grid size
        :param time: the time because level depends on it
        """

        h5_time_grp_key = "t"
        h5_filename = "EM_B.h5" # _____ TODO : could be another file

        import h5py
        data_file = h5py.File(os.path.join(self.path, h5_filename), "r")

        if time is None:
            time = float(list(data_file[h5_time_grp_key].keys())[0])

        hier = self._get_hierarchy(time, h5_filename)

        if level == 'finest':
            level = hier.finest_level(time)
        fac = np.power(hier.refinement_ratio, level)

        root_cell_width = np.asarray(data_file.attrs["cell_width"])

        return root_cell_width/fac
