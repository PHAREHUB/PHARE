from .hierarchy_utils import isFieldQty, field_qties, quantidic, refinement_ratio
from .patchdata import FieldData, ParticleData
from .patch import Patch
from .patchlevel import PatchLevel
from .hierarchy import PatchHierarchy
from ..particles import Particles
from ...core import gridlayout  #  import GridLayout, HybridGridLayoutFor
from ...core.box import Box

import numpy as np


def origin_from(patch):
    """PatchData.origin is a comma separated string, built C++ side by Point::str()"""
    return [float(v) for v in patch.origin.split(",")]


def particles_from(patch, layout):
    """Build Particles from a data wrangler patch, whose buffers are all flat."""
    v = np.asarray(patch.data.v)
    v = v.reshape(int(v.size / 3), 3)
    nbr = v.shape[0]

    dl = np.zeros((nbr, layout.ndim))
    for i in range(layout.ndim):
        dl[:, i] = layout.dl[i]

    return Particles(
        icells=np.asarray(patch.data.iCell).reshape(nbr, layout.ndim),
        deltas=np.asarray(patch.data.delta).reshape(nbr, layout.ndim),
        v=v,
        weights=np.asarray(patch.data.weight),
        charges=np.asarray(patch.data.charge),
        dl=dl,
    )


def make_layout_for(simulator, patch, qty, dl):
    model = "mhd" if str(qty).startswith("mhd") else "hybrid"
    box = Box(patch.lower, patch.upper)
    origin = origin_from(patch)
    if model == "hybrid":
        # check for particle quantity?
        return gridlayout.HybridGridLayoutFor(box, origin, dl, simulator.interp_order())
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
            wpatches = getters[qty]()
            for patch in wpatches:
                patch_datas = {}
                layout = make_layout_for(simulator, patch, qty, lvl_cell_width)
                pdata = FieldData(layout, field_qties[qty], np.asarray(patch.data))
                # the wrangler hands back a flat buffer: give it the layout's shape,
                # as an h5 read would, so the dataset is indexable like any other
                pdata.dataset = pdata.dataset.reshape(pdata.size)
                patch_datas[qty] = pdata
                patches[ilvl].append(Patch(patch_datas))

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
                domain_particles = particles_from(patch, layout)

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
                        layout = make_layout_for(
                            simulator, dwpatch, qty, lvl_cell_width
                        )
                        patchGhost_part = particles_from(dwpatch, layout)

                        box = Box(dwpatch.lower, dwpatch.upper)

                        # now search which of the already created patches has the same box
                        # once found we add the new particles to the ones already present

                        patch = [p for p in patches[ilvl] if p.box == box][0]
                        patch.patch_datas[pop + "_particles"].dataset.add(
                            patchGhost_part
                        )

        else:
            raise ValueError("{} is not a valid quantity".format(qty))

        patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

    return PatchHierarchy(
        [patch_levels], domain_box, refinement_ratio, times=[simulator.currentTime()]
    )
