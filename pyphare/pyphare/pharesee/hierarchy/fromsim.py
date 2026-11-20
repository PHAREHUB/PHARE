import numpy as np

from ...core import gridlayout  #  import GridLayout, HybridGridLayoutFor
from ...core.box import Box
from ...core.phare_utilities import refinement_ratio
from ..particles import Particles
from .hierarchy import PatchHierarchy
from .hierarchy_utils import field_qties, isFieldQty, quantidic
from .patch import Patch
from .patchdata import FieldData, ParticleData
from .patchlevel import PatchLevel


def make_layout_for(simulator, patch, qty, dl):
    model = "mhd" if str(qty).startswith("mhd") else "hybrid"
    box = Box(patch.lower, patch.upper)
    origin = np.asarray(patch.origin.split(","), dtype=float)  # "x,y,z" from C++
    if model == "hybrid":
        return gridlayout.HybridGridLayoutFor(
            box,
            origin,
            dl,
            simulator.interp_order(),
            is_particle_layout=(qty == "particles"),
        )
    return gridlayout.MHDGridLayoutFor(
        box, origin, dl, simulator.simulation.reconstruction
    )


def hierarchy_from_sim(simulator, qty, pop=""):
    dw = simulator.data_wrangler()
    nbr_levels = dw.getNumberOfLevels()
    patch_levels = {}

    root_cell_width = np.asarray(simulator.cell_width())
    domain_box = Box([0] * len(root_cell_width), simulator.domain_box())
    assert domain_box.ndim == len(simulator.domain_box())

    for ilvl in range(nbr_levels):
        lvl_cell_width = root_cell_width / refinement_ratio**ilvl

        patches = {ilvl: [] for ilvl in range(nbr_levels)}
        getters = quantidic(ilvl, dw)

        if isFieldQty(qty):
            if qty.startswith("flux_"):
                if pop == "":
                    raise ValueError("must specify pop argument for flux")
                wpatches = getters[qty](pop)
            else:
                wpatches = getters[qty]()
            for patch in wpatches:
                layout = make_layout_for(simulator, patch, qty, lvl_cell_width)
                pdata = FieldData(layout, field_qties[qty], None)
                # C++ gives the field flat, with its ghost nodes
                pdata.dataset = np.asarray(patch.data).reshape(pdata.size)
                patches[ilvl].append(Patch({qty: pdata}))

        elif qty == "particles":
            if pop == "":
                raise ValueError("must specify pop argument for particles")
            # here the getter returns a dict like this
            # {'protons': {'patchGhost': [<pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78970>,
            # <pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78f70>],
            # 'domain': [<pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78d70>,
            # <pybindlibs.cpp.PatchDataContiguousParticles_1 at 0x119f78770>]}}

            # domain particles are assumed to always be here
            # but patchGhost and levelGhost may not be, depending on the level

            populationdict = getters[qty](pop)[pop]

            dom_dw_patches = populationdict["domain"]
            for patch in dom_dw_patches:
                patch_datas = {}

                layout = make_layout_for(simulator, patch, qty, lvl_cell_width)
                domain_particles = _particles_from(patch.data, lvl_cell_width)

                patch_datas[pop + "_particles"] = ParticleData(
                    layout, domain_particles, pop
                )
                patches[ilvl].append(Patch(patch_datas))

            # ok now let's add the patchGhost if present
            # note that patchGhost patches may not be the same list as the
            # domain patches... since not all patches may not have patchGhost while they do have
            # domain... while looping on the patchGhost items, we need to search in
            # the already created patches which one to which add the patchGhost particles

            for ghostParticles in ["levelGhost"]:
                if ghostParticles in populationdict:
                    for dwpatch in populationdict[ghostParticles]:
                        patchGhost_part = _particles_from(dwpatch.data, lvl_cell_width)

                        box = Box(dwpatch.lower, dwpatch.upper)

                        # now search which of the already created patches has the same box
                        # once found we add the new particles to the ones already present

                        patch = next(p for p in patches[ilvl] if p.box == box)
                        patch.patch_datas[pop + "_particles"].dataset.add(
                            patchGhost_part
                        )

        else:
            raise ValueError(f"{qty} is not a valid quantity")

        patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

    return PatchHierarchy({simulator.currentTime(): patch_levels}, domain_box)


def _particles_from(contiguous, dl):
    """Particles from C++ ContiguousParticles, whose members are flat arrays"""
    ndim = len(dl)
    weights = np.asarray(contiguous.weight)
    return Particles(
        icells=np.asarray(contiguous.iCell).reshape(-1, ndim),
        deltas=np.asarray(contiguous.delta).reshape(-1, ndim),
        v=np.asarray(contiguous.v).reshape(-1, 3),
        weights=weights,
        dl=np.zeros((weights.size, ndim)) + dl,
    )
