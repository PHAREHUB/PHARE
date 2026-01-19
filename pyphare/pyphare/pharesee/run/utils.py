#
#
#

import numpy as np


from pyphare.core.gridlayout import yee_centering


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


def _ppp_to_ppp_domain_slicing(**kwargs):
    """
    return the slicing for (primal,primal,primal) to (primal,primal,primal)
    centering that is the centering of moments on a Yee grid
    """

    nb_ghosts = kwargs["nb_ghosts"]
    ndim = kwargs["ndim"]

    inner, _, _ = _inner_slices(nb_ghosts)

    inner_all = tuple([inner] * ndim)
    return inner_all, (inner_all,)


def _pdd_to_ppp_domain_slicing(**kwargs):
    """
    return the slicing for (dual,primal,primal) to (primal,primal,primal)
    centering that is the centering of Bx on a Yee grid
    """

    nb_ghosts = kwargs["nb_ghosts"]
    ndim = kwargs["ndim"]

    inner, inner_shift_left, inner_shift_right = _inner_slices(nb_ghosts)

    if ndim == 1:
        inner_all = tuple([inner] * ndim)
        return inner_all, (inner_all,)

    if ndim == 2:
        inner_all = tuple([inner] * ndim)
        return inner_all, ((inner, inner_shift_left), (inner, inner_shift_right))

    raise RuntimeError("dimension not yet implemented")


def _dpd_to_ppp_domain_slicing(**kwargs):
    """
    return the slicing for (dual,primal,primal) to (primal,primal,primal)
    centering that is the centering of By on a Yee grid
    """

    nb_ghosts = kwargs["nb_ghosts"]
    ndim = kwargs["ndim"]

    inner, inner_shift_left, inner_shift_right = _inner_slices(nb_ghosts)

    if ndim == 1:
        inner_all = tuple([inner] * ndim)
        return inner_all, (inner_shift_left, inner_shift_right)
    elif ndim == 2:
        inner_all = tuple([inner] * ndim)
        return inner_all, ((inner_shift_left, inner), (inner_shift_right, inner))
    else:
        raise RuntimeError("dimension not yet implemented")


def _ddp_to_ppp_domain_slicing(**kwargs):
    """
    return the slicing for (dual,primal,primal) to (primal,primal,primal)
    centering that is the centering of Bz on a Yee grid
    """

    nb_ghosts = kwargs["nb_ghosts"]
    ndim = kwargs["ndim"]

    inner, inner_shift_left, inner_shift_right = _inner_slices(nb_ghosts)

    if ndim == 1:
        inner_all = tuple([inner] * ndim)
        return inner_all, (inner_shift_left, inner_shift_right)
    elif ndim == 2:
        inner_all = tuple([inner] * ndim)
        return inner_all, (
            (inner_shift_left, inner_shift_left),
            (inner_shift_left, inner_shift_right),
            (inner_shift_right, inner_shift_left),
            (inner_shift_right, inner_shift_right),
        )
    else:
        raise RuntimeError("dimension not yet implemented")


def _dpp_to_ppp_domain_slicing(**kwargs):
    """
    return the slicing for (dual,primal,primal) to (primal,primal,primal)
    centering that is the centering of Ex on a Yee grid
    """

    nb_ghosts = kwargs["nb_ghosts"]
    ndim = kwargs["ndim"]

    inner, inner_shift_left, inner_shift_right = _inner_slices(nb_ghosts)

    if ndim == 1:
        inner_all = tuple([inner] * ndim)
        return inner_all, (inner_shift_left, inner_shift_right)
    elif ndim == 2:
        inner_all = tuple([inner] * ndim)
        return inner_all, ((inner_shift_left, inner), (inner_shift_right, inner))
    else:
        raise RuntimeError("dimension not yet implemented")


def _pdp_to_ppp_domain_slicing(**kwargs):
    """
    return the slicing for (dual,primal,primal) to (primal,primal,primal)
    centering that is the centering of Ey on a Yee grid
    """

    nb_ghosts = kwargs["nb_ghosts"]
    ndim = kwargs["ndim"]

    inner, inner_shift_left, inner_shift_right = _inner_slices(nb_ghosts)

    if ndim == 1:
        inner_all = tuple([inner] * ndim)
        return inner_all, (inner_all,)
    elif ndim == 2:
        inner_all = tuple([inner] * ndim)
        return inner_all, ((inner, inner_shift_left), (inner, inner_shift_right))
    else:
        raise RuntimeError("dimension not yet implemented")


def _ppd_to_ppp_domain_slicing(**kwargs):
    """
    return the slicing for (dual,primal,primal) to (primal,primal,primal)
    centering that is the centering of Ez on a Yee grid
    """

    nb_ghosts = kwargs["nb_ghosts"]
    ndim = kwargs["ndim"]

    inner, _, _ = _inner_slices(nb_ghosts)

    if ndim == 1:
        inner_all = tuple([inner] * ndim)
        return inner_all, (inner_all,)
    elif ndim == 2:
        inner_all = tuple([inner] * ndim)
        return inner_all, (inner_all,)
    else:
        raise RuntimeError("dimension not yet implemented")


slices_to_primal_ = {
    "primal_primal_primal": _ppp_to_ppp_domain_slicing,
    "primal_dual_dual": _pdd_to_ppp_domain_slicing,
    "dual_primal_dual": _dpd_to_ppp_domain_slicing,
    "dual_dual_primal": _ddp_to_ppp_domain_slicing,
    "dual_primal_primal": _dpp_to_ppp_domain_slicing,
    "primal_dual_primal": _pdp_to_ppp_domain_slicing,
    "primal_primal_dual": _ppd_to_ppp_domain_slicing,
}


def merge_centerings(pdname):
    from pyphare.core.gridlayout import directions

    return "_".join([yee_centering[d][pdname] for d in directions])


def slices_to_primal(pdname, **kwargs):
    return slices_to_primal_[merge_centerings(pdname)](**kwargs)


def _compute_to_primal(patch, **kwargs):
    """
    datasets have NaN in their ghosts... might need to be properly filled
    with their neighbors already properly projected on primal
    """

    ndim = patch.box.ndim

    pd_attrs = []
    for name, ref_pd in patch.patch_datas.items():
        nb_ghosts = ref_pd.layout.nbrGhosts(ref_pd.layout.interp_order, "primal")
        ref_ds = ref_pd.dataset

        if all(centering == "primal" for centering in ref_pd.centerings):
            pd_attrs.append({"name": name, "data": ref_pd})
            continue

        ds_shape = list(ref_ds.shape)
        for i in range(ndim):
            if ref_pd.centerings[i] == "dual":
                ds_shape[i] += 1

        # should be something else than nan values when the ghosts cells
        # will be filled with correct values coming from the neighbors
        ds_ = np.zeros(ds_shape)

        # inner is the slice containing the points that are updated
        # in the all_primal dataset
        # chunks is a tupls of all the slices coming from the initial dataset
        # that are needed to calculate the average for the all_primal dataset
        inner, chunks = slices_to_primal(name, nb_ghosts=nb_ghosts, ndim=ndim)

        for chunk in chunks:
            ds_[inner] = np.add(ds_[inner], ref_ds[chunk] / len(chunks))

        copy_pd = ref_pd.copy_as(
            ds_[inner], centering=["primal"] * ndim, ghosts_nbr=[0] * ndim
        )

        pd_attrs.append({"name": name, "data": copy_pd})

    return tuple(pd_attrs)


def _inner_slices(nb_ghosts):
    inner = slice(nb_ghosts, -nb_ghosts)
    inner_shift_left = slice(nb_ghosts - 1, -nb_ghosts)
    inner_shift_right = slice(nb_ghosts, -nb_ghosts + 1)

    return inner, inner_shift_left, inner_shift_right


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
    massDensity = patch["value"][:]
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
    N = patch["value"][:]

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
