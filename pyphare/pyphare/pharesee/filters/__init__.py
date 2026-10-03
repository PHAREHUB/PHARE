from copy import deepcopy

import numpy as np


def gaussian(hier, qty=None, sigma=2, time=None):
    from pyphare.pharesee.hierarchy import func

    finest = func.GetFinest(hier, time, qty)
    grids = deepcopy(finest)
    for key, grid in finest.items():
        grids[key] = gaussian_filter_uniform_grid(grid, sigma=sigma)
    return grids if len(grids) > 1 else next(iter(grids.values()))


def gaussian_filter_uniform_grid(grid, sigma=2):
    from scipy.ndimage import gaussian_filter

    ndim = grid.box.ndim
    nb_ghosts = grid.ghosts_nbr[0]
    ds = _fill_nan_ghosts(np.asarray(grid[:]))
    ds_ = np.full(list(ds.shape), np.nan)
    gf_ = gaussian_filter(ds, sigma=sigma)
    select = tuple([slice(nb_ghosts or None, -nb_ghosts or None) for _ in range(ndim)])
    ds_[select] = np.asarray(gf_[select])
    copy = deepcopy(grid)
    copy.dataset = ds_
    return copy


def _fill_nan_ghosts(ds):
    """
    ghost layers can be entirely NaN, e.g. outside the convex hull of a bilinear
    interpolation, they are replaced by the first non NaN layer
    """
    nan = np.isnan(ds)
    bounds = []
    for axis in range(ds.ndim):
        others = tuple(i for i in range(ds.ndim) if i != axis)
        valid = np.where(~nan.all(axis=others))[0]
        if len(valid) == 0:
            return ds
        bounds.append((valid[0], valid[-1]))

    inner = ds[tuple(slice(lo, hi + 1) for lo, hi in bounds)]
    pad = [(lo, n - 1 - hi) for (lo, hi), n in zip(bounds, ds.shape)]
    return np.pad(inner, pad, mode="edge")
