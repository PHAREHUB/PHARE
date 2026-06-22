#
#
#

import numpy as np
from dataclasses import dataclass


@dataclass
class FlatField:
    data: np.ndarray
    coords: np.ndarray
    qty: str
    ghosts_nbr: list
    centerings: list


def overlap_mask_1d(x, dl, level, qty):
    """
    return the mask for x where x is overlaped by the qty patch datas
    on the given level, assuming that this level is finer than the one of x

    :param x: 1d array containing the [x] position
    :param dl: list containing the grid steps where x is defined
    :param level: a given level associated to a hierarchy
    :param qty: ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Fx', 'Fy', 'Fz', 'Vx', 'Vy', 'Vz', 'rho']
    """

    is_overlaped = np.ones(x.shape[0], dtype=bool) * False

    for patch in level.patches:
        pdata = patch.patch_datas[qty]
        ghosts_nbr = pdata.ghosts_nbr

        fine_x = pdata.x[ghosts_nbr[0] - 1 : -ghosts_nbr[0] + 1]

        fine_dl = pdata.dl
        local_dl = dl

        if fine_dl[0] < local_dl[0]:
            xmin, xmax = fine_x.min(), fine_x.max()

            overlaped_idx = np.where((x > xmin) & (x < xmax))[0]

            is_overlaped[overlaped_idx] = True

        else:
            raise ValueError("level needs to have finer grid resolution than that of x")

    return is_overlaped


def overlap_mask_2d(x, y, dl, level, qty):
    """
    return the mask for x & y where ix & y are overlaped by the qty patch datas
    on the given level, assuming that this level is finer than the one of x & y
    important note : this mask is flatten

    :param x: 1d array containing the [x] position
    :param y: 1d array containing the [y] position
    :param dl: list containing the grid steps where x and y are defined
    :param level: a given level associated to a hierarchy
    :param qty: ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez', 'Fx', 'Fy', 'Fz', 'Vx', 'Vy', 'Vz', 'rho']
    """

    is_overlaped = np.ones([x.shape[0] * y.shape[0]], dtype=bool) * False

    for patch in level.patches:
        pdata = patch.patch_datas[qty]
        ghosts_nbr = pdata.ghosts_nbr

        if any(ghosts_nbr > 0):
            fine_x = pdata.x[ghosts_nbr[0] - 1 : -ghosts_nbr[0] + 1]
            fine_y = pdata.y[ghosts_nbr[1] - 1 : -ghosts_nbr[1] + 1]
        else:
            fine_x = pdata.x
            fine_y = pdata.y

        fine_dl = pdata.dl
        local_dl = dl

        if (fine_dl[0] < local_dl[0]) and (fine_dl[1] < local_dl[1]):
            xmin, xmax = fine_x.min(), fine_x.max()
            ymin, ymax = fine_y.min(), fine_y.max()

            xv, yv = np.meshgrid(x, y, indexing="ij")
            xf = xv.flatten()
            yf = yv.flatten()

            overlaped_idx = np.where(
                (xf > xmin) & (xf < xmax) & (yf > ymin) & (yf < ymax)
            )[0]

            is_overlaped[overlaped_idx] = True

        else:
            raise ValueError(
                "level needs to have finer grid resolution than that of x or y"
            )

    return is_overlaped


def flat_finest_field(hierarchy, qty, time=None):
    """
    Returns a FlatField with flattened data (shape [Npoints]), coordinates
    (shape [Npoints] in 1D, [Npoints, Ndim] in 2D), and field metadata.
    """
    dim = hierarchy.ndim

    if dim == 1:
        return flat_finest_field_1d(hierarchy, qty, time)
    elif dim == 2:
        return flat_finest_field_2d(hierarchy, qty, time)
    elif dim == 3:
        raise RuntimeError("Not implemented")
    else:
        raise ValueError("the dim of a hierarchy should be 1, 2 or 3")


def flat_finest_field_1d(hierarchy, qty, time=None):
    lvl = hierarchy.levels(time)
    finest_level = hierarchy.finest_level(time)
    flat_ghosts_nbr = None
    flat_centerings = None

    for ilvl in range(finest_level + 1)[::-1]:
        patches = lvl[ilvl].patches

        for ip, patch in enumerate(patches):
            pdata = patch[qty]
            data = pdata[:]
            x = pdata.x

            if ilvl == finest_level:
                if ip == 0:
                    final_data = data
                    final_x = x
                    flat_ghosts_nbr = list(pdata.ghosts_nbr)
                    flat_centerings = list(pdata.centerings)
                else:
                    final_data = np.concatenate((final_data, data))
                    final_x = np.concatenate((final_x, x))

            else:
                is_overlaped = overlap_mask_1d(
                    x, pdata.dl, hierarchy.level(ilvl + 1, time), qty
                )
                final_data = np.concatenate((final_data, data[~is_overlaped]))
                final_x = np.concatenate((final_x, x[~is_overlaped]))

    return FlatField(final_data, final_x, qty, flat_ghosts_nbr, flat_centerings)


def flat_finest_field_2d(hierarchy, qty, time=None):
    lvl = hierarchy.levels(time)
    finest_level = hierarchy.finest_level(time)
    flat_ghosts_nbr = None
    flat_centerings = None

    for ilvl in range(finest_level + 1)[::-1]:
        patches = lvl[ilvl].patches

        for ip, patch in enumerate(patches):
            pdata = patch.patch_datas[qty]
            data = pdata.dataset[:]
            x = pdata.x
            y = pdata.y

            xv, yv = np.meshgrid(x, y, indexing="ij")
            data_f = data.flatten()
            xv_f = xv.flatten()
            yv_f = yv.flatten()

            if ilvl == finest_level:
                if ip == 0:
                    final_data = data_f
                    tmp_x = xv_f
                    tmp_y = yv_f
                    flat_ghosts_nbr = list(pdata.ghosts_nbr)
                    flat_centerings = list(pdata.centerings)
                else:
                    final_data = np.concatenate((final_data, data_f))
                    tmp_x = np.concatenate((tmp_x, xv_f))
                    tmp_y = np.concatenate((tmp_y, yv_f))

            else:
                is_overlaped = overlap_mask_2d(
                    x, y, pdata.dl, hierarchy.level(ilvl + 1, time), qty
                )
                final_data = np.concatenate((final_data, data_f[~is_overlaped]))
                tmp_x = np.concatenate((tmp_x, xv_f[~is_overlaped]))
                tmp_y = np.concatenate((tmp_y, yv_f[~is_overlaped]))

    final_xy = np.stack((tmp_x, tmp_y), axis=1)

    return FlatField(final_data, final_xy, qty, flat_ghosts_nbr, flat_centerings)


def make_interpolator(flat_field, interp, domain, dl):
    """
    :param flat_field: FlatField with scattered data, coords, and field metadata
    :param interp: interpolation method
    :param domain: domain size per dimension
    :param dl: cell width per dimension (at the target finest level)

    Returns (interpolator, finest_coords) where finest_coords is the structured
    coordinate arrays on the finest regular grid.
    """
    from pyphare.core.gridlayout import yeeCoordsFor

    data = flat_field.data
    coords = flat_field.coords
    qty = flat_field.qty
    nbrGhosts = flat_field.ghosts_nbr
    centering = flat_field.centerings

    dim = coords.ndim
    withGhosts = any(g > 0 for g in nbrGhosts)

    if dim == 1:
        from scipy.interpolate import interp1d

        interpolator = interp1d(
            coords, data, kind=interp, fill_value="extrapolate", assume_sorted=False
        )

        nx = int(domain[0] / dl[0])
        kw_x = dict(centering=centering[0]) if centering else {}
        x = yeeCoordsFor(
            [0] * dim, nbrGhosts[0], dl, [nx], qty, "x", withGhosts=withGhosts, **kw_x
        )

        finest_coords = (x,)

    elif dim == 2:
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

        if interp == "nearest":
            interpolator = NearestNDInterpolator(coords, data)
        elif interp == "bilinear":
            interpolator = LinearNDInterpolator(coords, data)
        else:
            raise ValueError("interp can only be 'nearest' or 'bilinear'")

        nCells = [int(d / dl) for d, dl in zip(domain, dl)]

        kw_x = dict(centering=centering[0]) if centering else {}
        kw_y = dict(centering=centering[1]) if centering else {}
        x = yeeCoordsFor(
            [0] * dim, nbrGhosts[0], dl, nCells, qty, "x", withGhosts=withGhosts, **kw_x
        )
        y = yeeCoordsFor(
            [0] * dim, nbrGhosts[1], dl, nCells, qty, "y", withGhosts=withGhosts, **kw_y
        )

        finest_coords = (x, y)

    else:
        raise ValueError("make_interpolator is not yet 3d")

    return interpolator, finest_coords
