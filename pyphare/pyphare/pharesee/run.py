import os
import numpy as np

from .hierarchy import (
    compute_hier_from,
    flat_finest_field,
    hierarchy_from,
)
from pyphare.logger import getLogger

logger = getLogger(__name__)


def _current1d(by, bz, xby, xbz):
    # jx = 0
    # jy = -dxBz
    # jz = dxBy
    # the following hard-codes yee layout
    # which is not true in general
    # we should at some point provide proper
    # derivation routines in the gridlayout
    dx = xbz[1] - xbz[0]
    jy = np.zeros(by.size + 1)
    jy[1:-1] = -(bz[1:] - bz[:-1]) / dx
    dx = xby[1] - xby[0]
    jz = np.zeros(bz.size + 1)
    jz[1:-1] = (by[1:] - by[:-1]) / dx
    jy[0] = jy[1]
    jy[-1] = jy[-2]
    jz[0] = jz[1]
    jz[-1] = jz[-2]
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

    jx[:, 1:-1] = (bz[:, 1:] - bz[:, :-1]) / dy
    jy[1:-1, :] = -(bz[1:, :] - bz[:-1, :]) / dx
    jz[1:-1, 1:-1] = (by[1:, 1:-1] - by[:-1, 1:-1]) / dx - (
        bx[1:-1, 1:] - bx[1:-1, :-1]
    ) / dy

    jy[0, :] = jy[1, :]
    jy[:, 0] = jy[:, 1]
    jy[-1, :] = jy[-2, :]
    jy[:, -1] = jy[:, -2]

    jz[0, :] = jz[1, :]
    jz[:, 0] = jz[:, 1]
    jz[-1, :] = jz[-2, :]
    jz[:, -1] = jz[:, -2]

    jx[0, :] = jx[1, :]
    jx[:, 0] = jx[:, 1]
    jx[-1, :] = jx[-2, :]
    jx[:, -1] = jx[:, -2]

    return jx, jy, jz


def _compute_current(patchdatas, **kwargs):
    reference_pd = patchdatas["Bx"]  # take Bx as a reference, but could be any other

    ndim = reference_pd.box.ndim
    if ndim == 1:
        By = patchdatas["By"].dataset[:]
        xby = patchdatas["By"].x
        Bz = patchdatas["Bz"].dataset[:]
        xbz = patchdatas["Bz"].x
        Jy, Jz = _current1d(By, Bz, xby, xbz)
        return (
            {"name": "Jy", "data": Jy, "centering": "primal"},
            {"name": "Jz", "data": Jz, "centering": "primal"},
        )

    elif ndim == 2:
        Bx = patchdatas["Bx"].dataset[:]
        By = patchdatas["By"].dataset[:]
        Bz = patchdatas["Bz"].dataset[:]

        dx, dy = reference_pd.dl

        Jx, Jy, Jz = _current2d(Bx, By, Bz, dx, dy)

        components = ("Jx", "Jy", "Jz")
        centering = {
            component: [
                reference_pd.layout.centering[direction][component]
                for direction in ("X", "Y")
            ]
            for component in components
        }

        return (
            {"name": "Jx", "data": Jx, "centering": centering["Jx"]},
            {"name": "Jy", "data": Jy, "centering": centering["Jy"]},
            {"name": "Jz", "data": Jz, "centering": centering["Jz"]},
        )


def _divB2D(Bx, By, xBx, yBy):
    dxbx = (Bx[1:, :] - Bx[:-1, :]) / (xBx[1] - xBx[0])
    dyby = (By[:, 1:] - By[:, :-1]) / (yBy[1] - yBy[0])
    return dxbx + dyby


def _compute_divB(patchdatas, **kwargs):
    reference_pd = patchdatas["Bx"]  # take Bx as a reference, but could be any other
    ndim = reference_pd.box.ndim

    if ndim == 1:
        raise ValueError("divB is 0 by construction in 1D")

    elif ndim == 2:
        By = patchdatas["By"].dataset[:]
        Bx = patchdatas["Bx"].dataset[:]
        xBx = patchdatas["Bx"].x
        yBy = patchdatas["By"].y
        divB = _divB2D(Bx, By, xBx, yBy)

        return ({"name": "divB", "data": divB, "centering": ["dual", "dual"]},)

    else:
        raise RuntimeError("dimension not implemented")


def _get_rank(patchdatas, **kwargs):
    """
    make a field dataset cell centered coding the MPI rank
    rank is obtained from patch global id == "rank#local_patch_id"
    """
    from pyphare.core.box import grow

    reference_pd = patchdatas["Bx"]  # take Bx as a reference, but could be any other
    ndim = reference_pd.box.ndim
    pid = kwargs["id"]

    layout = reference_pd.layout
    centering = "dual"
    nbrGhosts = layout.nbrGhosts(layout.interp_order, centering)
    shape = grow(reference_pd.box, [nbrGhosts] * 2).shape

    if ndim == 1:
        pass

    elif ndim == 2:
        data = np.zeros(shape) + int(pid.strip("p").split("#")[0])
        return ({"name": "rank", "data": data, "centering": [centering] * 2},)
    else:
        raise RuntimeError("Not Implemented yet")


def _compute_pressure(patch_datas, **kwargs):
    Mxx = patch_datas["Mxx"].dataset[:]
    Mxy = patch_datas["Mxy"].dataset[:]
    Mxz = patch_datas["Mxz"].dataset[:]
    Myy = patch_datas["Myy"].dataset[:]
    Myz = patch_datas["Myz"].dataset[:]
    Mzz = patch_datas["Mzz"].dataset[:]
    massDensity = patch_datas["rho"].dataset[:]
    Vix = patch_datas["Vx"].dataset[:]
    Viy = patch_datas["Vy"].dataset[:]
    Viz = patch_datas["Vz"].dataset[:]

    Pxx = Mxx - Vix * Vix * massDensity
    Pxy = Mxy - Vix * Viy * massDensity
    Pxz = Mxz - Vix * Viz * massDensity
    Pyy = Myy - Viy * Viy * massDensity
    Pyz = Myz - Viy * Viz * massDensity
    Pzz = Mzz - Viz * Viz * massDensity

    return (
        {"name": "Pxx", "data": Pxx, "centering": ["primal", "primal"]},
        {"name": "Pxy", "data": Pxy, "centering": ["primal", "primal"]},
        {"name": "Pxz", "data": Pxz, "centering": ["primal", "primal"]},
        {"name": "Pyy", "data": Pyy, "centering": ["primal", "primal"]},
        {"name": "Pyz", "data": Pyz, "centering": ["primal", "primal"]},
        {"name": "Pzz", "data": Pzz, "centering": ["primal", "primal"]},
    )


def _compute_pop_pressure(patch_datas, **kwargs):
    """
    computes the pressure tensor for a given population
    this method is different from _compute_pressure in that:
        P = M - V*V*n*mass
    where V is the bulk velocity, n is the density and mass is the particle mass
    but for populations we don't have V we have the Flux.
    so the formula becomes:
        P = M - F*F/N * mass
    """
    popname = kwargs["popname"]
    Mxx = patch_datas[popname + "_Mxx"].dataset[:]
    Mxy = patch_datas[popname + "_Mxy"].dataset[:]
    Mxz = patch_datas[popname + "_Mxz"].dataset[:]
    Myy = patch_datas[popname + "_Myy"].dataset[:]
    Myz = patch_datas[popname + "_Myz"].dataset[:]
    Mzz = patch_datas[popname + "_Mzz"].dataset[:]
    Fx = patch_datas[popname + "_Fx"].dataset[:]
    Fy = patch_datas[popname + "_Fy"].dataset[:]
    Fz = patch_datas[popname + "_Fz"].dataset[:]
    N = patch_datas[popname + "_rho"].dataset[:]

    mass = kwargs["mass"]

    Pxx = Mxx - Fx * Fx * mass / N
    Pxy = Mxy - Fx * Fy * mass / N
    Pxz = Mxz - Fx * Fz * mass / N
    Pyy = Myy - Fy * Fy * mass / N
    Pyz = Myz - Fy * Fz * mass / N
    Pzz = Mzz - Fz * Fz * mass / N

    return (
        {"name": popname + "_Pxx", "data": Pxx, "centering": ["primal", "primal"]},
        {"name": popname + "_Pxy", "data": Pxy, "centering": ["primal", "primal"]},
        {"name": popname + "_Pxz", "data": Pxz, "centering": ["primal", "primal"]},
        {"name": popname + "_Pyy", "data": Pyy, "centering": ["primal", "primal"]},
        {"name": popname + "_Pyz", "data": Pyz, "centering": ["primal", "primal"]},
        {"name": popname + "_Pzz", "data": Pzz, "centering": ["primal", "primal"]},
    )


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

        interpolator = interp1d(
            coords, data, kind=interp, fill_value="extrapolate", assume_sorted=False
        )

        nx = 1 + int(domain[0] / dl[0])

        x = yeeCoordsFor([0] * dim, nbrGhosts, dl, [nx], qty, "x")
        finest_coords = (x,)

    elif dim == 2:
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

        if interp == "nearest":
            interpolator = NearestNDInterpolator(coords, data)
        elif interp == "bilinear":
            interpolator = LinearNDInterpolator(coords, data)
        else:
            raise ValueError("interp can only be 'nearest' or 'bilinear'")

        nCells = [1 + int(d / dl) for d, dl in zip(domain, dl)]
        x = yeeCoordsFor([0] * dim, nbrGhosts, dl, nCells, qty, "x")
        y = yeeCoordsFor([0] * dim, nbrGhosts, dl, nCells, qty, "y")
        # x = np.arange(0, domain[0]+dl[0], dl[0])
        # y = np.arange(0, domain[1]+dl[1], dl[1])
        finest_coords = (x, y)

    else:
        raise ValueError("make_interpolator is not yet 3d")

    return interpolator, finest_coords


class Run:
    def __init__(self, path, single_hier_for_all_quantities=False):
        self.path = path
        self.single_hier_for_all_quantities = single_hier_for_all_quantities
        self.hier = None  # only used if single_hier_for_all_quantities == True

    def _get_hierarchy(self, time, filename, hier=None):
        t = "{:.10f}".format(time)

        def _get_hier(h):
            return hierarchy_from(
                h5_filename=os.path.join(self.path, filename), time=t, hier=h
            )

        if self.single_hier_for_all_quantities:
            self.hier = _get_hier(self.hier)
            return self.hier
        return _get_hier(hier)

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
            nbrGhosts = list(hierarchy.level(0).patches[0].patch_datas.values())[
                0
            ].ghosts_nbr
            merged_qties = {}
            for qty in hierarchy.quantities():
                data, coords = flat_finest_field(hierarchy, qty, time=time)
                merged_qties[qty] = make_interpolator(
                    data, coords, interp, domain, dl, qty, nbrGhosts
                )
            return merged_qties
        else:
            return hierarchy

    def GetTags(self, time, merged=False):
        hier = self._get_hierarchy(time, "tags.h5")
        return self._get(hier, time, merged, "nearest")

    def GetB(self, time, merged=False, interp="nearest"):
        hier = self._get_hierarchy(time, "EM_B.h5")
        return self._get(hier, time, merged, interp)

    def GetE(self, time, merged=False, interp="nearest"):
        hier = self._get_hierarchy(time, "EM_E.h5")
        return self._get(hier, time, merged, interp)

    def GetMassDensity(self, time, merged=False, interp="nearest"):
        hier = self._get_hierarchy(time, "ions_mass_density.h5")
        return self._get(hier, time, merged, interp)

    def GetNi(self, time, merged=False, interp="nearest"):
        hier = self._get_hierarchy(time, "ions_density.h5")
        return self._get(hier, time, merged, interp)

    def GetN(self, time, pop_name, merged=False, interp="nearest"):
        hier = self._get_hierarchy(time, f"ions_pop_{pop_name}_density.h5")
        return self._get(hier, time, merged, interp)

    def GetVi(self, time, merged=False, interp="nearest"):
        hier = self._get_hierarchy(time, "ions_bulkVelocity.h5")
        return self._get(hier, time, merged, interp)

    def GetFlux(self, time, pop_name, merged=False, interp="nearest"):
        hier = self._get_hierarchy(time, f"ions_pop_{pop_name}_flux.h5")
        return self._get(hier, time, merged, interp)

    def GetPressure(self, time, pop_name, merged=False, interp="nearest"):
        M = self._get_hierarchy(time, f"ions_pop_{pop_name}_momentum_tensor.h5")
        V = self.GetFlux(time, pop_name)
        N = self.GetN(time, pop_name)
        P = compute_hier_from(
            _compute_pop_pressure,
            (M, V, N),
            popname=pop_name,
            mass=self.GetMass(pop_name),
        )
        return self._get(P, time, merged, interp)

    def GetPi(self, time, merged=False, interp="nearest"):
        M = self._get_hierarchy(time, f"ions_momentum_tensor.h5")
        massDensity = self.GetMassDensity(time)
        Vi = self._get_hierarchy(time, f"ions_bulkVelocity.h5")
        Pi = compute_hier_from(_compute_pressure, (M, massDensity, Vi))
        return self._get(Pi, time, merged, interp)

    def GetJ(self, time, merged=False, interp="nearest"):
        B = self.GetB(time)
        J = compute_hier_from(_compute_current, B)
        return self._get(J, time, merged, interp)

    def GetDivB(self, time, merged=False, interp="nearest"):
        B = self.GetB(time)
        db = compute_hier_from(_compute_divB, B)
        return self._get(db, time, merged, interp)

    def GetRanks(self, time, merged=False, interp="nearest"):
        """
        returns a hierarchy of MPI ranks
        takes the information from magnetic field diagnostics arbitrarily
        this fails if the magnetic field is not written and could at some point
        be replace by a search of any available diag at the requested time.
        """
        B = self.GetB(time)
        ranks = compute_hier_from(_get_rank, B)
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
        list_of_qty = ["density", "flux", "domain", "levelGhost", "patchGhost"]
        list_of_mass = []

        import h5py

        for qty in list_of_qty:
            file = os.path.join(self.path, "ions_pop_{}_{}.h5".format(pop_name, qty))
            if os.path.isfile(file):
                h5_file = h5py.File(file, "r")
                list_of_mass.append(h5_file.attrs["pop_mass"])

        assert all(m == list_of_mass[0] for m in list_of_mass)

        return list_of_mass[0]

    def GetDomainSize(self):
        h5_filename = "EM_B.h5"  # _____ TODO : could be another file

        import h5py

        data_file = h5py.File(os.path.join(self.path, h5_filename), "r")

        root_cell_width = np.asarray(data_file.attrs["cell_width"])

        return (data_file.attrs["domain_box"] + 1) * root_cell_width

    def GetDl(self, level="finest", time=None):
        """
        gives the ndarray containing the grid sizes at the given time
        for the hierarchy defined in the given run, and for the given level
        (default is 'finest', but can also be a int)

        :param level: the level at which get the associated grid size
        :param time: the time because level depends on it
        """

        h5_time_grp_key = "t"
        h5_filename = "EM_B.h5"  # _____ TODO : could be another file

        import h5py

        data_file = h5py.File(os.path.join(self.path, h5_filename), "r")

        if time is None:
            time = float(list(data_file[h5_time_grp_key].keys())[0])

        hier = self._get_hierarchy(time, h5_filename)

        if level == "finest":
            level = hier.finest_level(time)
        fac = np.power(hier.refinement_ratio, level)

        root_cell_width = np.asarray(data_file.attrs["cell_width"])

        return root_cell_width / fac

    def GetAllAvailableQties(self, time=0, pops=[]):
        assert self.single_hier_for_all_quantities  # can't work otherwise

        def _try(fn, *args, **kwargs):
            try:
                fn(*args, **kwargs)
            except FileNotFoundError:
                # normal to not have a diagnostic if not requested
                logger.debug(f"No file for function {fn.__name__}")

        _try(self.GetParticles, time, pops)
        _try(self.GetB, time)
        _try(self.GetE, time)
        _try(self.GetNi, time)
        _try(self.GetVi, time)

        for pop in pops:
            _try(self.GetFlux, time, pop)
            _try(self.GetN, time, pop)

        return self.hier
