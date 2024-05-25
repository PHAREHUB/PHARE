from .hierarchy import PatchHierarchy
from .patchdata import FieldData
from .patchlevel import PatchLevel
from .patch import Patch
from ...core.phare_utilities import listify
from ...core.phare_utilities import refinement_ratio
from ...core import box as boxm

import numpy as np

field_qties = {
    "EM_B_x": "Bx",
    "EM_B_y": "By",
    "EM_B_z": "Bz",
    "EM_E_x": "Ex",
    "EM_E_y": "Ey",
    "EM_E_z": "Ez",
    "flux_x": "Fx",
    "flux_y": "Fy",
    "flux_z": "Fz",
    "bulkVelocity_x": "Vx",
    "bulkVelocity_y": "Vy",
    "bulkVelocity_z": "Vz",
    "momentum_tensor_xx": "Mxx",
    "momentum_tensor_xy": "Mxy",
    "momentum_tensor_xz": "Mxz",
    "momentum_tensor_yy": "Myy",
    "momentum_tensor_yz": "Myz",
    "momentum_tensor_zz": "Mzz",
    "density": "rho",
    "mass_density": "rho",
    "tags": "tags",
}


def are_compatible_hierarchies(hierarchies):
    return True


def merge_particles(hierarchy):
    for time, patch_levels in hierarchy.time_hier.items():
        for ilvl, plvl in patch_levels.items():
            for ip, patch in enumerate(plvl.patches):
                pdatas = patch.patch_datas
                domain_pdata = [
                    (pdname, pd) for pdname, pd in pdatas.items() if "domain" in pdname
                ][0]

                pghost_pdatas = [
                    (pdname, pd)
                    for pdname, pd in pdatas.items()
                    if "patchGhost" in pdname
                ]
                lghost_pdatas = [
                    (pdname, pd)
                    for pdname, pd in pdatas.items()
                    if "levelGhost" in pdname
                ]

                pghost_pdata = pghost_pdatas[0] if pghost_pdatas else None
                lghost_pdata = lghost_pdatas[0] if lghost_pdatas else None

                if pghost_pdata is not None:
                    domain_pdata[1].dataset.add(pghost_pdata[1].dataset)
                    del pdatas[pghost_pdata[0]]

                if lghost_pdata is not None:
                    domain_pdata[1].dataset.add(lghost_pdata[1].dataset)
                    del pdatas[lghost_pdata[0]]

                popname = domain_pdata[0].split("_")[0]
                pdatas[popname + "_particles"] = pdatas[domain_pdata[0]]
                del pdatas[domain_pdata[0]]


def getPatch(hier, point):
    """
    convenience function mainly for debug. returns a dict
    {ilevel:patch}  for patches in which the given point is
    """
    patches = {}
    counts = {ilvl: 0 for ilvl in hier.levels().keys()}
    for ilvl, lvl in hier.levels().items():
        for p in lvl.patches:
            px, py = point
            dx, dy = p.layout.dl
            ix = int(px / dx)
            iy = int(py / dy)
            if (ix, iy) in p.box:
                patches[ilvl] = p
                counts[ilvl] += 1
    for k, v in counts.items():
        if v > 1:
            print("error : ", k, v)
            raise RuntimeError("more than one patch found for point")
    return patches


def compute_hier_from(compute, hierarchies, **kwargs):
    """return a new hierarchy using the callback 'compute' on all patchdatas
    of the given hierarchies
    """
    if not are_compatible_hierarchies(hierarchies):
        raise RuntimeError("hierarchies are not compatible")

    hierarchies = listify(hierarchies)
    reference_hier = hierarchies[0]
    domain_box = reference_hier.domain_box
    patch_levels = {}
    for ilvl in range(reference_hier.levelNbr()):
        patch_levels[ilvl] = PatchLevel(
            ilvl, new_patches_from(compute, hierarchies, ilvl, **kwargs)
        )

    assert len(reference_hier.time_hier) == 1  # only single time hierarchies now
    t = list(reference_hier.time_hier.keys())[0]
    return PatchHierarchy(patch_levels, domain_box, refinement_ratio, time=t)


def extract_patchdatas(hierarchies, ilvl, ipatch):
    """
    returns a dict {patchdata_name:patchdata} from a list of hierarchies for patch ipatch at level ilvl
    """
    patches = [h.patch_levels[ilvl].patches[ipatch] for h in hierarchies]
    patch_datas = {
        pdname: pd for p in patches for pdname, pd in list(p.patch_datas.items())
    }
    return patch_datas


def new_patchdatas_from(compute, patchdatas, layout, id, **kwargs):
    new_patch_datas = {}
    datas = compute(patchdatas, id=id, **kwargs)
    for data in datas:
        pd = FieldData(layout, data["name"], data["data"], centering=data["centering"])
        new_patch_datas[data["name"]] = pd
    return new_patch_datas


def new_patches_from(compute, hierarchies, ilvl, **kwargs):
    reference_hier = hierarchies[0]
    new_patches = []
    patch_nbr = len(reference_hier.patch_levels[ilvl].patches)
    for ip in range(patch_nbr):
        current_patch = reference_hier.patch_levels[ilvl].patches[ip]
        layout = current_patch.layout
        patch_datas = extract_patchdatas(hierarchies, ilvl, ip)
        new_patch_datas = new_patchdatas_from(
            compute, patch_datas, layout, id=current_patch.id, **kwargs
        )
        new_patches.append(Patch(new_patch_datas, current_patch.id))
    return new_patches


def is_root_lvl(patch_level):
    return patch_level.level_number == 0


def quantidic(ilvl, wrangler):
    pl = wrangler.getPatchLevel(ilvl)

    return {
        "density": pl.getDensity,
        "bulkVelocity_x": pl.getVix,
        "bulkVelocity_y": pl.getViy,
        "bulkVelocity_z": pl.getViz,
        "EM_B_x": pl.getBx,
        "EM_B_y": pl.getBy,
        "EM_B_z": pl.getBz,
        "EM_E_x": pl.getEx,
        "EM_E_y": pl.getEy,
        "EM_E_z": pl.getEz,
        "flux_x": pl.getFx,
        "flux_y": pl.getFy,
        "flux_z": pl.getFz,
        "particles": pl.getParticles,
    }


def isFieldQty(qty):
    return qty in (
        "density",
        "bulkVelocity_x",
        "bulkVelocity_y",
        "bulkVelocity_z",
        "EM_B_x",
        "EM_B_y",
        "EM_B_z",
        "EM_E_x",
        "EM_E_y",
        "EM_E_z",
        "flux_x",
        "flux_y",
        "flux_z",
        "momentumTensor_xx",
        "momentumTensor_xy",
        "momentumTensor_xz",
        "momentumTensor_yy",
        "momentumTensor_yz",
        "momentumTensor_zz",
    )


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

        fine_x = pdata.x[ghosts_nbr[0] - 1 : -ghosts_nbr[0] + 1]
        fine_y = pdata.y[ghosts_nbr[1] - 1 : -ghosts_nbr[1] + 1]

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
    returns 2 flattened arrays containing the data (with shape [Npoints])
    and the coordinates (with shape [Npoints, Ndim]) for the given
    hierarchy of qty.

    :param hierarchy: the hierarchy where qty is defined
    :param qty: the field (eg "Bx") that we want
    """

    dim = hierarchy.ndim

    if dim == 1:
        return flat_finest_field_1d(hierarchy, qty, time)
    elif dim == 2:
        return flat_finest_field_2d(hierarchy, qty, time)
    elif dim == 3:
        raise RuntimeError("Not implemented")
        # return flat_finest_field_3d(hierarchy, qty, time)
    else:
        raise ValueError("the dim of a hierarchy should be 1, 2 or 3")


def flat_finest_field_1d(hierarchy, qty, time=None):
    lvl = hierarchy.levels(time)

    for ilvl in range(hierarchy.finest_level(time) + 1)[::-1]:
        patches = lvl[ilvl].patches

        for ip, patch in enumerate(patches):
            pdata = patch.patch_datas[qty]

            # all but 1 ghost nodes are removed in order to limit
            # the overlapping, but to keep enough point to avoid
            # any extrapolation for the interpolator
            needed_points = pdata.ghosts_nbr - 1

            # data = pdata.dataset[patch.box] # TODO : once PR 551 will be merged...
            data = pdata.dataset[needed_points[0] : -needed_points[0]]
            x = pdata.x[needed_points[0] : -needed_points[0]]

            if ilvl == hierarchy.finest_level(time):
                if ip == 0:
                    final_data = data
                    final_x = x
                else:
                    final_data = np.concatenate((final_data, data))
                    final_x = np.concatenate((final_x, x))

            else:
                is_overlaped = overlap_mask_1d(
                    x, pdata.dl, hierarchy.level(ilvl + 1, time), qty
                )

                finest_data = data[~is_overlaped]
                finest_x = x[~is_overlaped]

                final_data = np.concatenate((final_data, finest_data))
                final_x = np.concatenate((final_x, finest_x))

    return final_data, final_x


def flat_finest_field_2d(hierarchy, qty, time=None):
    lvl = hierarchy.levels(time)

    for ilvl in range(hierarchy.finest_level(time) + 1)[::-1]:
        patches = lvl[ilvl].patches

        for ip, patch in enumerate(patches):
            pdata = patch.patch_datas[qty]

            # all but 1 ghost nodes are removed in order to limit
            # the overlapping, but to keep enough point to avoid
            # any extrapolation for the interpolator
            needed_points = pdata.ghosts_nbr - 1

            # data = pdata.dataset[patch.box] # TODO : once PR 551 will be merged...
            data = pdata.dataset[
                needed_points[0] : -needed_points[0],
                needed_points[1] : -needed_points[1],
            ]
            x = pdata.x[needed_points[0] : -needed_points[0]]
            y = pdata.y[needed_points[1] : -needed_points[1]]

            xv, yv = np.meshgrid(x, y, indexing="ij")

            data_f = data.flatten()
            xv_f = xv.flatten()
            yv_f = yv.flatten()

            if ilvl == hierarchy.finest_level(time):
                if ip == 0:
                    final_data = data_f
                    tmp_x = xv_f
                    tmp_y = yv_f
                else:
                    final_data = np.concatenate((final_data, data_f))
                    tmp_x = np.concatenate((tmp_x, xv_f))
                    tmp_y = np.concatenate((tmp_y, yv_f))

            else:
                is_overlaped = overlap_mask_2d(
                    x, y, pdata.dl, hierarchy.level(ilvl + 1, time), qty
                )

                finest_data = data_f[~is_overlaped]
                finest_x = xv_f[~is_overlaped]
                finest_y = yv_f[~is_overlaped]

                final_data = np.concatenate((final_data, finest_data))
                tmp_x = np.concatenate((tmp_x, finest_x))
                tmp_y = np.concatenate((tmp_y, finest_y))

    final_xy = np.stack((tmp_x, tmp_y), axis=1)

    return final_data, final_xy
