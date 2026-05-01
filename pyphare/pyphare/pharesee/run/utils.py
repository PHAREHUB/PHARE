#
#
#

import numpy as np


from pyphare.pharesee.hierarchy.interpolation import make_interpolator  # noqa: F401


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


def _compute_current(patch, **kwargs):
    reference_pd = patch["Bx"]  # take Bx as a reference, but could be any other

    ndim = reference_pd.box.ndim
    if ndim == 1:
        By = patch["By"]
        xby = By.x
        Bz = patch["Bz"]
        xbz = Bz.x
        Jy, Jz = _current1d(By[:], Bz[:], xby, xbz)
        return (
            {"name": "Jy", "data": By.copy_as(Jy, centering="primal")},
            {"name": "Jz", "data": Bz.copy_as(Jz, centering="primal")},
        )

    elif ndim == 2:
        Bx = patch["Bx"]
        By = patch["By"]
        Bz = patch["Bz"]

        dx, dy = reference_pd.dl
        Jx, Jy, Jz = _current2d(Bx[:], By[:], Bz[:], dx, dy)

        centering = {
            component: [
                reference_pd.layout.centering[direction][component]
                for direction in ("X", "Y")
            ]
            for component in ("Jx", "Jy", "Jz")
        }

        return (
            {"name": "Jx", "data": Bx.copy_as(Jx, centering=centering["Jx"])},
            {"name": "Jy", "data": By.copy_as(Jy, centering=centering["Jy"])},
            {"name": "Jz", "data": Bz.copy_as(Jz, centering=centering["Jz"])},
        )

    raise RuntimeError("dimension not implemented")


def _divB2D(Bx, By, xBx, yBy):
    dxbx = (Bx[1:, :] - Bx[:-1, :]) / (xBx[1] - xBx[0])
    dyby = (By[:, 1:] - By[:, :-1]) / (yBy[1] - yBy[0])
    return dxbx + dyby


def _compute_divB(patch, **kwargs):
    reference_pd = patch["Bx"]  # take Bx as a reference, but could be any other
    ndim = reference_pd.box.ndim

    if ndim == 1:
        raise ValueError("divB is 0 by construction in 1D")

    centering = ["dual"] * ndim

    if ndim == 2:
        By = patch["By"].dataset[:]
        Bx = patch["Bx"].dataset[:]
        xBx = patch["Bx"].x
        yBy = patch["By"].y
        divB = reference_pd.copy_as(_divB2D(Bx, By, xBx, yBy), centering=centering)
        return ({"name": "divB", "data": divB},)

    raise RuntimeError("dimension not implemented")


from pyphare.pharesee.hierarchy.primal import (  # noqa: F401, E402
    _inner_slices,
    _ppp_to_ppp_domain_slicing,
    _pdd_to_ppp_domain_slicing,
    _dpd_to_ppp_domain_slicing,
    _ddp_to_ppp_domain_slicing,
    _dpp_to_ppp_domain_slicing,
    _pdp_to_ppp_domain_slicing,
    _ppd_to_ppp_domain_slicing,
    slices_to_primal_,
    merge_centerings,
    slices_to_primal,
    _compute_to_primal_patch_data,
    _compute_to_primal,
)


def _get_rank(patch, **kwargs):
    """
    make a field dataset cell centered coding the MPI rank
    rank is obtained from patch global id == "rank#local_patch_id"
    """
    from pyphare.core.box import grow

    reference_pd = patch["Bx"]  # Bx as a ref, but could be any other
    ndim = reference_pd.box.ndim

    layout = reference_pd.layout
    centering = ["dual"] * ndim
    nbrGhosts = layout.nbrGhosts(layout.interp_order, centering)
    shape = grow(reference_pd.box, [nbrGhosts] * 2).shape

    if ndim == 1:
        pass

    elif ndim == 2:
        data = np.zeros(shape) + int(patch.id.strip("p").split("#")[0])
        pd = reference_pd.copy_as(data, centering=centering)
        return ({"name": "rank", "data": pd},)
    else:
        raise RuntimeError("Not Implemented yet")


def _compute_pressure(patch, **kwargs):
    Mxx = patch["Mxx"]
    Mxy = patch["Mxy"]
    Mxz = patch["Mxz"]
    Myy = patch["Myy"]
    Myz = patch["Myz"]
    Mzz = patch["Mzz"]
    massDensity = patch["rho"][:]
    Vix = patch["Vx"][:]
    Viy = patch["Vy"][:]
    Viz = patch["Vz"][:]

    Pxx = Mxx[:] - Vix * Vix * massDensity
    Pxy = Mxy[:] - Vix * Viy * massDensity
    Pxz = Mxz[:] - Vix * Viz * massDensity
    Pyy = Myy[:] - Viy * Viy * massDensity
    Pyz = Myz[:] - Viy * Viz * massDensity
    Pzz = Mzz[:] - Viz * Viz * massDensity

    return (
        {"name": "Pxx", "data": Mxx.copy_as(Pxx)},
        {"name": "Pxy", "data": Mxy.copy_as(Pxy)},
        {"name": "Pxz", "data": Mxz.copy_as(Pxz)},
        {"name": "Pyy", "data": Myy.copy_as(Pyy)},
        {"name": "Pyz", "data": Myz.copy_as(Pyz)},
        {"name": "Pzz", "data": Mzz.copy_as(Pzz)},
    )


def _compute_pop_pressure_xx(patch, mass, popname):
    Mxx = patch[popname + "_Mxx"]
    Fx = patch[popname + "_Fx"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pxx": Mxx.copy_as(Mxx[:] - Fx * Fx * mass / N)}


def _compute_pop_pressure_xy(patch, mass, popname):
    Mxy = patch[popname + "_Mxy"]
    Fx = patch[popname + "_Fx"][:]
    Fy = patch[popname + "_Fy"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pxy": Mxy.copy_as(Mxy[:] - Fx * Fy * mass / N)}


def _compute_pop_pressure_xz(patch, mass, popname):
    Mxz = patch[popname + "_Mxz"]
    Fx = patch[popname + "_Fx"][:]
    Fz = patch[popname + "_Fz"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pxz": Mxz.copy_as(Mxz[:] - Fx * Fz * mass / N)}


def _compute_pop_pressure_yy(patch, mass, popname):
    Myy = patch[popname + "_Myy"]
    Fy = patch[popname + "_Fy"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pyy": Myy.copy_as(Myy[:] - Fy * Fy * mass / N)}


def _compute_pop_pressure_yz(patch, mass, popname):
    Myz = patch[popname + "_Myz"]
    Fy = patch[popname + "_Fy"][:]
    Fz = patch[popname + "_Fz"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pyz": Myz.copy_as(Myz[:] - Fy * Fz * mass / N)}


def _compute_pop_pressure_zz(patch, mass, popname):
    Mzz = patch[popname + "_Mzz"]
    Fz = patch[popname + "_Fz"][:]
    N = patch[popname + "_rho"][:]
    return {popname + "_Pzz": Mzz.copy_as(Mzz[:] - Fz * Fz * mass / N)}


def _compute_pop_pressure(patch, **kwargs):
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
    Mxx = patch[popname + "_Mxx"]
    Mxy = patch[popname + "_Mxy"]
    Mxz = patch[popname + "_Mxz"]
    Myy = patch[popname + "_Myy"]
    Myz = patch[popname + "_Myz"]
    Mzz = patch[popname + "_Mzz"]
    Fx = patch["x"][:]
    Fy = patch["y"][:]
    Fz = patch["z"][:]
    N = patch["rho"][:]

    mass = kwargs["mass"]

    Pxx = Mxx[:] - Fx * Fx * mass / N
    Pxy = Mxy[:] - Fx * Fy * mass / N
    Pxz = Mxz[:] - Fx * Fz * mass / N
    Pyy = Myy[:] - Fy * Fy * mass / N
    Pyz = Myz[:] - Fy * Fz * mass / N
    Pzz = Mzz[:] - Fz * Fz * mass / N

    return (
        {"name": popname + "_Pxx", "data": Mxx.copy_as(Pxx)},
        {"name": popname + "_Pxy", "data": Mxy.copy_as(Pxy)},
        {"name": popname + "_Pxz", "data": Mxz.copy_as(Pxz)},
        {"name": popname + "_Pyy", "data": Myy.copy_as(Pyy)},
        {"name": popname + "_Pyz", "data": Myz.copy_as(Pyz)},
        {"name": popname + "_Pzz", "data": Mzz.copy_as(Pzz)},
    )
