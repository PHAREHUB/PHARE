#
#
#

import numpy as np


def _compute_gaussian_filter_on_scalarfield(patch_datas, **kwargs):
    from scipy.ndimage import gaussian_filter

    ndim = patch_datas["value"].box.ndim
    nb_ghosts = patch_datas["value"].ghosts_nbr[0]
    sigma = kwargs["sigma"]
    ds = np.asarray(patch_datas["value"][:])

    ds_ = np.full(list(ds.shape), np.nan)

    gf_ = gaussian_filter(ds, sigma=sigma)
    select = tuple([slice(nb_ghosts, -nb_ghosts) for _ in range(ndim)])
    ds_[select] = np.asarray(gf_[select])

    return (
        {"name": "value", "data": ds_, "centering": patch_datas["value"].centerings},
    )


def _compute_gaussian_filter_on_vectorfield(patch_datas, **kwargs):
    from scipy.ndimage import gaussian_filter

    ref_name = next(iter(patch_datas.keys()))

    ndim = patch_datas[ref_name].box.ndim
    nb_ghosts = patch_datas[ref_name].ghosts_nbr[0]
    sigma = kwargs["sigma"]
    ds_x = np.asarray(patch_datas["x"][:])
    ds_y = np.asarray(patch_datas["y"][:])
    ds_z = np.asarray(patch_datas["z"][:])

    dsx_ = np.full(list(ds_x.shape), np.nan)
    dsy_ = np.full(list(ds_y.shape), np.nan)
    dsz_ = np.full(list(ds_z.shape), np.nan)

    gfx_ = gaussian_filter(ds_x, sigma=sigma)
    gfy_ = gaussian_filter(ds_y, sigma=sigma)
    gfz_ = gaussian_filter(ds_z, sigma=sigma)

    select = tuple([slice(nb_ghosts, -nb_ghosts) for _ in range(ndim)])

    dsx_[select] = np.asarray(gfx_[select])
    dsy_[select] = np.asarray(gfy_[select])
    dsz_[select] = np.asarray(gfz_[select])

    return (
        {"name": "x", "data": dsx_, "centering": patch_datas["x"].centerings},
        {"name": "y", "data": dsy_, "centering": patch_datas["y"].centerings},
        {"name": "z", "data": dsz_, "centering": patch_datas["z"].centerings},
    )
