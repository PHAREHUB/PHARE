from dataclasses import dataclass, field
from copy import deepcopy
import numpy as np

from typing import Any, List, Tuple

from .hierarchy import PatchHierarchy, format_timestamp
from .patchdata import FieldData, ParticleData
from .patchlevel import PatchLevel
from .patch import Patch
from ...core.box import Box
from ...core.gridlayout import GridLayout
from ...core.phare_utilities import listify
from ...core.phare_utilities import refinement_ratio
from pyphare.core import phare_utilities as phut


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


def nbr_ranks(hier):
    """
    returns the number of mpi ranks used in the given hierarchy
    """
    max_rank = 0
    t0 = hier.times()[0]
    for _, lvl in hier.levels(t0).items():
        for patch in lvl.patches:
            rank = patch.attrs["mpi_rank"]
            if rank > max_rank:
                max_rank = rank
    return max_rank


def patch_per_rank(hier):
    """
    returns the number of patch per mpi rank for each time step
    """
    nbranks = nbr_ranks(hier)
    ppr = {}
    for t in hier.times():
        ppr[t] = {ir: 0 for ir in np.arange(nbranks + 1)}
        for _, lvl in hier.levels(t).items():
            for patch in lvl.patches:
                ppr[t][patch.attrs["mpi_rank"]] += 1

    return ppr


def are_compatible_hierarchies(hierarchies):
    ref = hierarchies[0]
    same_box = True
    same_selection = True
    same_files = True
    same_times = True
    for hier in hierarchies[1:]:
        same_box = same_box and (hier.domain_box == ref.domain_box)
        same_selection = same_selection and (hier.selection_box == ref.selection_box)
        same_files = same_files and (hier.data_files.keys() == ref.data_files.keys())
        same_times = same_times and (hier.time_hier.keys() == ref.time_hier.keys())
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
    return hierarchy


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

    hierarchies = listify(hierarchies)
    if not are_compatible_hierarchies(hierarchies):
        raise RuntimeError("hierarchies are not compatible")
    reference_hier = hierarchies[0]
    domain_box = reference_hier.domain_box
    patch_levels_per_time = []
    for t in reference_hier.times():
        patch_levels = {}
        for ilvl in range(reference_hier.levelNbr()):
            patch_levels[ilvl] = PatchLevel(
                ilvl, new_patches_from(compute, hierarchies, ilvl, t, **kwargs)
            )
        patch_levels_per_time.append(patch_levels)
    return PatchHierarchy(
        patch_levels_per_time,
        domain_box,
        refinement_ratio,
        times=reference_hier.times(),
    )


def extract_patchdatas(hierarchies, ilvl, t, ipatch):
    """
    returns a dict {patchdata_name:patchdata} from a list of hierarchies for patch ipatch at level ilvl
    """
    patches = [h.level(ilvl, time=t).patches[ipatch] for h in hierarchies]
    patch_datas = {
        pdname: pd for p in patches for pdname, pd in list(p.patch_datas.items())
    }
    return patch_datas


def new_patchdatas_from(compute, patchdatas, layout, id, **kwargs):
    new_patch_datas = {}
    datas = compute(patchdatas, patch_id=id, **kwargs)
    for data in datas:
        pd = FieldData(layout, data["name"], data["data"], centering=data["centering"])
        new_patch_datas[data["name"]] = pd
    return new_patch_datas


def new_patches_from(compute, hierarchies, ilvl, t, **kwargs):
    reference_hier = hierarchies[0]
    new_patches = []
    ref_patches = reference_hier.level(ilvl, time=t).patches
    for ip, current_patch in enumerate(ref_patches):
        layout = current_patch.layout
        patch_datas = extract_patchdatas(hierarchies, ilvl, t, ip)
        new_patch_datas = new_patchdatas_from(
            compute, patch_datas, layout, id=current_patch.id, **kwargs
        )
        new_patches.append(
            Patch(new_patch_datas, current_patch.id, attrs=current_patch.attrs)
        )
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


def flat_finest_field(hierarchy, qty, time=None, neghosts=1):
    """
    returns 2 flattened arrays containing the data (with shape [Npoints])
    and the coordinates (with shape [Npoints, Ndim]) for the given
    hierarchy of qty.

    :param hierarchy: the hierarchy where qty is defined
    :param qty: the field (eg "Bx") that we want
    """

    dim = hierarchy.ndim

    if dim == 1:
        return flat_finest_field_1d(hierarchy, qty, time, neghosts)
    elif dim == 2:
        return flat_finest_field_2d(hierarchy, qty, time)
    elif dim == 3:
        raise RuntimeError("Not implemented")
        # return flat_finest_field_3d(hierarchy, qty, time)
    else:
        raise ValueError("the dim of a hierarchy should be 1, 2 or 3")


def flat_finest_field_1d(hierarchy, qty, time=None, neghosts=1):
    lvl = hierarchy.levels(time)

    for ilvl in range(hierarchy.finest_level(time) + 1)[::-1]:
        patches = lvl[ilvl].patches

        for ip, patch in enumerate(patches):
            pdata = patch.patch_datas[qty]

            # all but 1 ghost nodes are removed in order to limit
            # the overlapping, but to keep enough point to avoid
            # any extrapolation for the interpolator
            needed_points = pdata.ghosts_nbr - neghosts

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


def compute_rename(patch_datas, **kwargs):
    new_names = kwargs["new_names"]
    pd_attrs = []

    for new_name, pd_name in zip(new_names, patch_datas):
        pd_attrs.append(
            {
                "name": new_name,
                "data": patch_datas[pd_name].dataset,
                "centering": patch_datas[pd_name].centerings,
            }
        )

    return tuple(pd_attrs)


def rename(hierarchy, names):
    return compute_hier_from(compute_rename, hierarchy, new_names=names)


def _compute_mul(patch_datas, **kwargs):
    names = kwargs["names"]

    # multiplication of a scalarField or VectorField by a scalar
    if "other" in kwargs:
        other = kwargs["other"]
        pd_attrs = []

        for name, pd_name in zip(names, patch_datas):
            pd_attrs.append(
                {
                    "name": name,
                    "data": other * patch_datas[pd_name].dataset[:],
                    "centering": patch_datas[pd_name].centerings,
                }
            )
    # multiplication of 2 scalarField
    elif "self_value" in patch_datas:
        dset_value = (
            patch_datas["self_value"].dataset[:] * patch_datas["other_value"].dataset[:]
        )
        pd_attrs = (
            {
                "name": "value",
                "data": dset_value,
                "centering": patch_datas["self_value"].centerings,
            },
        )

    return tuple(pd_attrs)


def _compute_add(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset_x = patch_datas["self_x"].dataset[:] + patch_datas["other_x"].dataset[:]
    dset_y = patch_datas["self_y"].dataset[:] + patch_datas["other_y"].dataset[:]
    dset_z = patch_datas["self_z"].dataset[:] + patch_datas["other_z"].dataset[:]

    return (
        {"name": "x", "data": dset_x, "centering": patch_datas[ref_name].centerings},
        {"name": "y", "data": dset_y, "centering": patch_datas[ref_name].centerings},
        {"name": "z", "data": dset_z, "centering": patch_datas[ref_name].centerings},
    )


def _compute_sub(patch_datas, **kwargs):
    ref_name = next(iter(patch_datas.keys()))

    dset_x = patch_datas["self_x"].dataset[:] - patch_datas["other_x"].dataset[:]
    dset_y = patch_datas["self_y"].dataset[:] - patch_datas["other_y"].dataset[:]
    dset_z = patch_datas["self_z"].dataset[:] - patch_datas["other_z"].dataset[:]

    return (
        {"name": "x", "data": dset_x, "centering": patch_datas[ref_name].centerings},
        {"name": "y", "data": dset_y, "centering": patch_datas[ref_name].centerings},
        {"name": "z", "data": dset_z, "centering": patch_datas[ref_name].centerings},
    )


def _compute_neg(patch_datas, **kwargs):
    names = kwargs["new_names"]
    pd_attrs = []

    for name in names:
        pd_attrs.append(
            {
                "name": name,
                "data": -patch_datas[name].dataset[:],
                "centering": patch_datas[name].centerings,
            }
        )

    return tuple(pd_attrs)


def _compute_truediv(patch_datas, **kwargs):
    names = kwargs["res_names"]
    pd_attrs = []

    # the denominator is a scalar field which name is "value"
    # hence the associated patchdata has to be removed from the list
    # of patchdata from the vectorField of the numerator
    left_ops = {k: v for k, v in patch_datas.items() if k != "value"}
    right_op = patch_datas["value"]
    for name, left_op in zip(names, left_ops.values()):
        pd_attrs.append(
            {
                "name": name,
                "data": left_op.dataset / right_op.dataset,
                "centering": patch_datas[name].centerings,
            }
        )

    return tuple(pd_attrs)


def _compute_scalardiv(patch_datas, **kwargs):
    names = kwargs["res_names"]
    scalar = kwargs["scalar"]
    pd_attrs = []

    left_ops = {k: v for k, v in patch_datas.items()}
    for name, left_op in zip(names, left_ops.values()):
        pd_attrs.append(
            {
                "name": name,
                "data": left_op.dataset / scalar,
                "centering": patch_datas[name].centerings,
            }
        )

    return tuple(pd_attrs)


@dataclass
class EqualityReport:
    failed: List[Tuple[str, Any, Any]] = field(default_factory=lambda: [])

    def __bool__(self):
        return not self.failed

    def __repr__(self):
        for msg, ref, cmp in self:
            print(msg)
            try:
                if type(ref) is FieldData:
                    phut.assert_fp_any_all_close(ref[:], cmp[:], atol=1e-16)
            except AssertionError as e:
                print(e)
        return self.failed[0][0]

    def __call__(self, reason, ref=None, cmp=None):
        self.failed.append((reason, ref, cmp))
        return self

    def __getitem__(self, idx):
        return (self.failed[idx][1], self.failed[idx][2])

    def __iter__(self):
        return self.failed.__iter__()

    def __reversed__(self):
        return reversed(self.failed)


def hierarchy_compare(this, that, atol=1e-16):
    eqr = EqualityReport()

    if not isinstance(this, PatchHierarchy) or not isinstance(that, PatchHierarchy):
        return eqr("class type mismatch")

    if this.ndim != that.ndim or this.domain_box != that.domain_box:
        return eqr("dimensional mismatch")

    if this.time_hier.keys() != that.time_hier.keys():
        return eqr("timesteps mismatch")

    for tidx in this.times():
        patch_levels_ref = this.time_hier[tidx]
        patch_levels_cmp = that.time_hier[tidx]

        if patch_levels_ref.keys() != patch_levels_cmp.keys():
            return eqr("levels mismatch")

        for level_idx in patch_levels_cmp.keys():
            patch_level_ref = patch_levels_ref[level_idx]
            patch_level_cmp = patch_levels_cmp[level_idx]

            for patch_idx in range(len(patch_level_cmp.patches)):
                patch_ref = patch_level_ref.patches[patch_idx]
                patch_cmp = patch_level_cmp.patches[patch_idx]

                if patch_ref.patch_datas.keys() != patch_cmp.patch_datas.keys():
                    return eqr("data keys mismatch")

                for patch_data_key in patch_ref.patch_datas.keys():
                    patch_data_ref = patch_ref.patch_datas[patch_data_key]
                    patch_data_cmp = patch_cmp.patch_datas[patch_data_key]

                    if not patch_data_cmp.compare(patch_data_ref, atol=atol):
                        msg = f"data mismatch: {type(patch_data_ref).__name__} {patch_data_key}"
                        eqr(msg, patch_data_cmp, patch_data_ref)

                if not eqr:
                    return eqr

    return eqr


def single_patch_for_LO(hier, qties=None, skip=None):
    def _skip(qty):
        return (qties is not None and qty not in qties) or (
            skip is not None and qty in skip
        )

    cier = deepcopy(hier)
    sim = hier.sim
    layout = GridLayout(
        Box(sim.origin, sim.cells), sim.origin, sim.dl, interp_order=sim.interp_order
    )
    p0 = Patch(patch_datas={}, patch_id="", layout=layout)
    for t in cier.times():
        cier.time_hier[format_timestamp(t)] = {0: cier.level(0, t)}
        cier.level(0, t).patches = [deepcopy(p0)]
        l0_pds = cier.level(0, t).patches[0].patch_datas
        for k, v in hier.level(0, t).patches[0].patch_datas.items():
            if _skip(k):
                continue
            if isinstance(v, FieldData):
                l0_pds[k] = FieldData(
                    layout, v.field_name, None, centering=v.centerings
                )
                l0_pds[k].dataset = np.zeros(l0_pds[k].size)
                patch_box = hier.level(0, t).patches[0].box
                l0_pds[k][patch_box] = v[patch_box]

            elif isinstance(v, ParticleData):
                l0_pds[k] = deepcopy(v)
            else:
                raise RuntimeError("unexpected state")

        for patch in hier.level(0, t).patches[1:]:
            for k, v in patch.patch_datas.items():
                if _skip(k):
                    continue
                if isinstance(v, FieldData):
                    l0_pds[k][patch.box] = v[patch.box]
                elif isinstance(v, ParticleData):
                    l0_pds[k].dataset.add(v.dataset)
                else:
                    raise RuntimeError("unexpected state")
    return cier
