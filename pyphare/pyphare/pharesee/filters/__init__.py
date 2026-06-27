#
#
#

import numpy as np
from copy import deepcopy


def gaussian(hier, qty=None, sigma=2, **kwargs):
    from pyphare.pharesee import hierarchy as harch

    time = harch.func.GetTime(hier)
    finest = harch.func.GetFinest(hier, time, qty)
    grids = deepcopy(finest)
    for key, grid in finest.items():
        grids[key] = gaussian_filter_uniform_grid(grid, sigma=sigma)
    return grids if len(grids) > 1 else next(iter(grids.values()))


def gaussian_filter_uniform_grid(grid, sigma=2):
    from scipy.ndimage import gaussian_filter

    ndim = grid.box.ndim
    nb_ghosts = grid.ghosts_nbr[0]
    ds = np.asarray(grid[:])
    ds_ = np.full(list(ds.shape), np.nan)
    gf_ = gaussian_filter(ds, sigma=sigma)
    select = tuple([slice(nb_ghosts or None, -nb_ghosts or None) for _ in range(ndim)])
    ds_[select] = np.asarray(gf_[select])
    copy = deepcopy(grid)
    copy.dataset = ds_
    return copy
