#
#
#

import numpy as np

from pyphare.core.gridlayout import yee_centering


def _inner_slices(nb_ghosts):
    if nb_ghosts < 1:
        raise ValueError("Do not call with zero or negative ghosts, redundant!")

    inner = slice(nb_ghosts, -nb_ghosts)
    inner_shift_left = slice(nb_ghosts - 1 or None, -nb_ghosts)
    inner_shift_right = slice(nb_ghosts, -nb_ghosts + 1 or None)
    return inner, inner_shift_left, inner_shift_right


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


def _compute_to_primal_patch_data(ref_pd, name=None):
    ndim = ref_pd.box.ndim
    name = name or ref_pd.name
    nb_ghosts = ref_pd.layout.nbrGhosts(ref_pd.layout.interp_order, "primal")

    ref_ds = ref_pd.dataset

    should_skip = all(  # vtkhdf is all primal with no ghosts
        [ref_pd.centerings == ["primal"] * ndim, not any(ref_pd.ghosts_nbr)]
    )
    if should_skip:
        return ref_pd

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

    return copy_pd


def _compute_to_primal(patch, **kwargs):
    """
    datasets have NaN in their ghosts... might need to be properly filled
    with their neighbors already properly projected on primal
    """

    if not hasattr(patch, "patch_datas"):  # assume is patch_data
        return _compute_to_primal_patch_data(patch)

    pd_attrs = []
    for name, ref_pd in patch.patch_datas.items():
        pd_attrs.append(
            {"name": name, "data": _compute_to_primal_patch_data(ref_pd, name)}
        )

    return tuple(pd_attrs)
