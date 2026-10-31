from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pyphare.core import phare_utilities as phut

from ...core.box import Box
from ...core.gridlayout import GridLayout
from ...core.phare_utilities import listify
from .hierarchy import HierarchyAccessor, PatchHierarchy, format_timestamp
from .patch import Patch
from .patchdata import FieldData, ParticleData, align_ghosts
from .patchlevel import PatchLevel

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
    "charge_flux_x": "CFx",
    "charge_flux_y": "CFy",
    "charge_flux_z": "CFz",
    "momentum_tensor_xx": "Mxx",
    "momentum_tensor_xy": "Mxy",
    "momentum_tensor_xz": "Mxz",
    "momentum_tensor_yy": "Myy",
    "momentum_tensor_yz": "Myz",
    "momentum_tensor_zz": "Mzz",
    "density": "rho",
    "mass_density": "rho",
    "charge_density": "rho",
    "tags": "tags",
    "rho": "mhdRho",
    "V_x": "mhdVx",
    "V_y": "mhdVy",
    "V_z": "mhdVz",
    "P": "mhdP",
    "rhoV_x": "mhdRhoVx",
    "rhoV_y": "mhdRhoVy",
    "rhoV_z": "mhdRhoVz",
    "Etot": "mhdEtot",
}


def nbr_ranks(hier):
    """
    returns the number of mpi ranks used in the given hierarchy
    """
    max_rank = 0
    t0 = hier.times()[0]
    for lvl in hier.levels(t0).values():
        for patch in lvl.patches:
            rank = patch.attrs["mpi_rank"]
            max_rank = max(max_rank, rank)
    return max_rank


def patch_per_rank(hier):
    """
    returns the number of patch per mpi rank for each time step
    """
    nbranks = nbr_ranks(hier)
    ppr = {}
    for t in hier.times():
        ppr[t] = dict.fromkeys(np.arange(nbranks + 1), 0)
        for lvl in hier.levels(t).values():
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
    for patch_levels in hierarchy.time_hier.values():
        for plvl in patch_levels.values():
            for patch in plvl.patches:
                pdatas = patch.patch_datas
                domain_pdata = next(
                    (pdname, pd) for pdname, pd in pdatas.items() if "domain" in pdname
                )

                lghost_pdatas = [
                    (pdname, pd)
                    for pdname, pd in pdatas.items()
                    if "levelGhost" in pdname
                ]

                lghost_pdata = lghost_pdatas[0] if lghost_pdatas else None

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
    counts = dict.fromkeys(hier.levels(), 0)
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
    time_hier = {}
    for t in reference_hier.times():
        patch_levels = {}
        for ilvl in range(reference_hier.levelNbr()):
            patch_levels[ilvl] = PatchLevel(
                ilvl, new_patches_from(compute, hierarchies, ilvl, t, **kwargs)
            )
        time_hier[t] = patch_levels

    return PatchHierarchy(
        time_hier,
        domain_box,
        data_files=reference_hier.data_files,
        selection_box=reference_hier.selection_box,
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


def new_patchdatas_from(compute, patch, **kwargs):
    if not patch.patch_datas:  # empty patches are kept, nothing to compute
        return {}
    return {data.name: data for data in compute(patch, **kwargs)}


def new_patches_from(compute, hierarchies, ilvl, t, **kwargs):
    reference_hier = hierarchies[0]
    new_patches = []
    ref_patches = reference_hier.level(ilvl, time=t).patches
    for ip, ref_patch in enumerate(ref_patches):
        patch_datas = extract_patchdatas(hierarchies, ilvl, t, ip)
        if len(hierarchies) > 1:
            pdatas = list(patch_datas.values())
            aligned = align_ghosts(pdatas)
            if aligned is not pdatas:
                patch_datas = dict(zip(patch_datas.keys(), aligned))
        patch = ref_patch.copy_as(patch_datas=patch_datas)
        accessor = HierarchyAccessor(reference_hier, t, ilvl, ip)
        new_patches.append(
            ref_patch.copy_as(
                patch_datas=new_patchdatas_from(
                    compute, patch, accessor=accessor, **kwargs
                ),
            )
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




@dataclass
class EqualityReport:
    failed: list[tuple[str, Any, Any]] = field(default_factory=list)

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
        return self.failed[0][0] if self.failed else "=="

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

        for level_idx in patch_levels_cmp:
            patch_level_ref = patch_levels_ref[level_idx]
            patch_level_cmp = patch_levels_cmp[level_idx]

            for patch_idx in range(len(patch_level_cmp.patches)):
                patch_ref = patch_level_ref.patches[patch_idx]
                patch_cmp = patch_level_cmp.patches[patch_idx]

                if patch_ref.box != patch_cmp.box:
                    return eqr("patch box mismatch", patch_ref.box, patch_cmp.box)

                if patch_ref.patch_datas.keys() != patch_cmp.patch_datas.keys():
                    return eqr("data keys mismatch")

                for patch_data_key in patch_ref.patch_datas:
                    patch_data_ref = patch_ref.patch_datas[patch_data_key]
                    patch_data_cmp = patch_cmp.patch_datas[patch_data_key]

                    ret = patch_data_ref.compare(patch_data_cmp, atol=atol)
                    if not bool(ret):
                        msg = f"data mismatch: {type(patch_data_ref).__name__} {patch_data_key}"
                        if type(ret) is not bool:
                            msg += "\n" + str(ret)
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
    origin = [0] * len(sim.cells)
    box = Box(origin, sim.cells)

    def layout_of(patch_data):
        """whole-domain layout for that patch data: same as its own layout but for the
        box. Ghosts come from the data rather than from interp_order, which an MHD run
        does not have."""
        return GridLayout(
            box,
            origin,
            sim.dl,
            interp_order=patch_data.layout.interp_order,
            ghosts_nbr=patch_data.layout.ghosts_nbr,
        )

    p0 = Patch(patch_datas={}, patch_id="", box=box)
    for t in cier.times():
        cier.time_hier[format_timestamp(t)] = {0: cier.level(0, t)}
        cier.level(0, t).patches = [deepcopy(p0)]
        l0_pds = cier.level(0, t).patches[0].patch_datas
        for k, v in hier.level(0, t).patches[0].patch_datas.items():
            if _skip(k):
                continue
            if isinstance(v, FieldData):
                l0_pds[k] = FieldData(
                    layout_of(v),
                    v.name,
                    None,
                    centering=v.centerings,
                    ghosts_nbr=v.ghosts_nbr,
                )
                l0_pds[k].dataset = np.zeros(l0_pds[k].size)
                patch_box = hier.level(0, t).patches[0].box
                l0_pds[k][patch_box] = v[patch_box]

            elif isinstance(v, ParticleData):
                l0_pds[k] = deepcopy(v)
            else:
                raise TypeError("unexpected state")

        for patch in hier.level(0, t).patches[1:]:
            for k, v in patch.patch_datas.items():
                if _skip(k):
                    continue
                if isinstance(v, FieldData):
                    l0_pds[k][patch.box] = v[patch.box]
                elif isinstance(v, ParticleData):
                    l0_pds[k].dataset.add(v.dataset)
                else:
                    raise TypeError("unexpected state")
    return cier
