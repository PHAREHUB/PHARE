import numpy as np

from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from

from scipy.ndimage import gaussian_filter



def _compute_gaussian_filter_on_scalarfield(patch_datas, **kwargs):
    from scipy.ndimage import gaussian_filter

    ndim = patch_datas["value"].box.ndim
    nb_ghosts = kwargs["nb_ghosts"]
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
    nb_ghosts = kwargs["nb_ghosts"]
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


def gFilt(hier, **kwargs):
    sigma = kwargs.get("sigma", 2)

    # time0 = list(hier.times())[0]
    # level0 = 0
    # p0 = 0
    # pd0 = hier.levels(time0)[level0].patches[p0].patch_datas
    # key0 = list(pd0.keys())[0]
    # nb_ghosts = np.max(pd0[key0].ghosts_nbr)

    nb_ghosts = np.max(list(hier.level(0).patches[0].patch_datas.values())[0].ghosts_nbr)

    if nb_ghosts < sigma :
        print("nb_ghosts ({0}) < sigma ({1}) : your gaussian filter might be dirty".format(nb_ghosts, sigma))

    if hier.ndim == 1:
        if isinstance(hier, ScalarField) :
            h = compute_hier_from(_compute_gaussian_filter_on_scalarfield, hier, nb_ghosts=nb_ghosts, sigma=sigma)
            return ScalarField(h)
        elif isinstance(hier, VectorField) :
            h = compute_hier_from(_compute_gaussian_filter_on_vectorfield, hier, nb_ghosts=nb_ghosts, sigma=sigma)
            return VectorField(h)
        else:
            return NotImplemented
    else:
        return NotImplemented


def gF(hier, **kwargs):
    sigma = kwargs.get("sigma", 4)
    if sigma == 1:
        raise ValueError("sigma value has to be at least 2")
    h_ = hier.__deepcopy__(memo={})
    ndim = hier.ndim
    n_pad = 4*sigma+1
    # The gaussian filter is calculated on the box extended by
    # n_pad. Hence, the number of points is large enough so that the value
    # at the last point of the real box is as equal as possible to the one
    # at the first point of the next box...

    for time in h_.times():
        interp_ = hier.interpol(time)
        for lvl in h_.levels(time).values():
            for patch in lvl.patches:
                names = list(patch.patch_datas.keys())
                box = patch.box

                for name in names:
                    pdata = patch.patch_datas[name]
                    nb_ghosts = pdata.ghosts_nbr
                    if not n_pad > nb_ghosts:
                        raise ValueError('sigma value is too small')

                    r_ = []
                    for i in range(ndim):
                        s_ = np.arange(box.lower[i]-n_pad, box.upper[i]+2+n_pad)*pdata.dl[i]
                        r_.append(s_)

                    func, _ = interp_[name]

                    if ndim == 1:
                        data = func(r_[0])
                    elif ndim == 2:
                        data = func(r_[0], r_[1])
                    elif ndim == 3:
                        data = func(r_[0], r_[1], r_[2])
                    else:
                        raise ValueError('unvalid dimension')

                    gf_ = gaussian_filter(data, sigma=sigma)
                    select = tuple([slice(n_pad-nb_ghosts[i], -n_pad+nb_ghosts[i]) for i in range(ndim)])
                    pdata.dataset = np.asarray(gf_[select])

    return h_


def peakIds(hier, **kwargs):
    from scipy.signal import find_peaks

    times = list(hier.times())
    if len(times) == 1:
        time = times[0]
    else:
        raise ValueError('multiple time is not yet implemented')

    pks_ = np.array([])
    hgs_ = np.array([])

    names_ = kwargs.pop("names", None)
    if names_ is None:
        names_ = list(hier.levels(time)[0].patches[0].patch_datas.keys())

    ph_ = kwargs.get('peak_heights', None)
    if ph_ is None:
        raise ValueError("the kwarg 'peak_heights' is mandatory for now...")

    for lvl in hier.levels(time).values():
        for patch in lvl.patches:
            for name in names_:
                pdata = patch.patch_datas[name]
                ds = np.asarray(pdata.dataset)
                pks = find_peaks(ds, **kwargs)
                for pk, hg in zip(pks[0], pks[1]['peak_heights']):
                    pks_ = np.append(pks_, np.add(np.multiply(pk, patch.dl), patch.origin))
                    hgs_ = np.append(hgs_, hg)

    return pks_, hgs_
