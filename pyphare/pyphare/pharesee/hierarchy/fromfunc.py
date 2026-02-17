#
#
#

import numpy as np

from pyphare.pharesee.hierarchy.hierarchy_utils import compute_hier_from


def ions_mass_density_func1d(x, **kwargs):
    masses = kwargs["masses"]  # list of float : the ion pop masses
    densities = kwargs["densities"]  # list of callable : the ion pop density profiles

    assert len(masses) == len(densities)
    funcs = np.zeros((x.size, len(masses)))

    for i, (mass, density) in enumerate(zip(masses, densities)):
        funcs[:, i] = mass * density(x)

    return funcs.sum(axis=1)


def ions_charge_density_func1d(x, **kwargs):
    charges = kwargs["charges"]  # list of float : the ion pop charges
    densities = kwargs["densities"]  # list of callable : the ion pop density profiles

    assert len(charges) == len(densities)

    funcs = np.zeros((x.size, len(charges)))

    for i, (charge, density) in enumerate(zip(charges, densities)):
        funcs[:, i] = charge * density(x)

    return funcs.sum(axis=1)


def hierarchy_from_func1d(func, hier, **kwargs):
    assert hier.ndim == 1

    def compute_(patch, **kwargs):
        ref_name = next(iter(patch.patch_datas.keys()))
        x_ = patch[ref_name].x

        new_patch_data = patch[ref_name].copy_as(
            func(x_, **kwargs), ghosts_nbr=patch[ref_name].ghosts_nbr
        )

        return ({"name": "value", "data": new_patch_data},)

    return compute_hier_from(compute_, hier, **kwargs)


def ions_mass_density_func2d(x, y, **kwargs):
    masses = kwargs["masses"]  # list of float : the ion pop masses
    densities = kwargs["densities"]  # list of callable : the ion pop density profiles

    yv, xv = np.meshgrid(y, x)

    assert len(masses) == len(densities)
    funcs = np.zeros((x.size, y.size, len(masses)))

    for i, (mass, density) in enumerate(zip(masses, densities)):
        funcs[:, :, i] = mass * density(xv, yv)

    return funcs.sum(axis=2)


def ions_charge_density_func2d(x, y, **kwargs):
    charges = kwargs["charges"]  # list of float : the ion pop charges
    densities = kwargs["densities"]  # list of callable : the ion pop density profiles

    yv, xv = np.meshgrid(y, x)

    assert len(charges) == len(densities)
    funcs = np.zeros((x.size, y.size, len(charges)))

    for i, (charge, density) in enumerate(zip(charges, densities)):
        funcs[:, :, i] = charge * density(xv, yv)

    return funcs.sum(axis=2)


def hierarchy_from_func2d(func, hier, **kwargs):
    assert hier.ndim == 2

    def compute_(patch, **kwargs):
        ref_name = next(iter(patch.patch_datas.keys()))
        x_ = patch[ref_name].x
        y_ = patch[ref_name].y

        new_patch_data = patch[ref_name].copy_as(
            func(x_, y_, **kwargs), ghosts_nbr=patch[ref_name].ghosts_nbr
        )

        return ({"name": "value", "data": new_patch_data},)

    return compute_hier_from(compute_, hier, **kwargs)


def hierarchy_from_func(func, hier, **kwargs):
    if hier.ndim == 1:
        return hierarchy_from_func1d(func, hier, **kwargs)
    if hier.ndim == 2:
        return hierarchy_from_func2d(func, hier, **kwargs)
