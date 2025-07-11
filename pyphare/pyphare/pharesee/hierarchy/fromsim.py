from .hierarchy_utils import isFieldQty, field_qties, quantidic, refinement_ratio
from .patchdata import FieldData, ParticleData
from .patch import Patch
from .patchlevel import PatchLevel
from .hierarchy import PatchHierarchy
from ..particles import Particles
from ...core.gridlayout import GridLayout
from ...core.box import Box

import numpy as np


def hierarchy_from_sim(simulator, qty, pop=""):
    dw = simulator.data_wrangler()
    nbr_levels = dw.getNumberOfLevels()
    patch_levels = {}

    root_cell_width = simulator.cell_width()
    domain_box = Box([0] * len(root_cell_width), simulator.domain_box())
    assert len(domain_box.ndim) == len(simulator.domain_box().ndim)

    for ilvl in range(nbr_levels):
        lvl_cell_width = root_cell_width / refinement_ratio**ilvl

        patches = {ilvl: [] for ilvl in range(nbr_levels)}
        getters = quantidic(ilvl, dw)

        if isFieldQty(qty):
            wpatches = getters[qty]()
            for patch in wpatches:
                patch_datas = {}
                lower = patch.lower
                upper = patch.upper
                origin = patch.origin
                layout = GridLayout(
                    Box(lower, upper),
                    origin,
                    lvl_cell_width,
                    interp_order=simulator.interporder(),
                )
                pdata = FieldData(layout, field_qties[qty], patch.data)
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

                lower = patch.lower
                upper = patch.upper
                origin = patch.origin
                layout = GridLayout(
                    Box(lower, upper),
                    origin,
                    lvl_cell_width,
                    interp_order=simulator.interp_order(),
                )
                v = np.asarray(patch.data.v).reshape(int(len(patch.data.v) / 3), 3)

                domain_particles = Particles(
                    icells=np.asarray(patch.data.iCell),
                    deltas=np.asarray(patch.data.delta),
                    v=v,
                    weights=np.asarray(patch.data.weight),
                    charges=np.asarray(patch.data.charge),
                )

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
                        v = np.asarray(dwpatch.data.v)
                        s = v.size
                        v = v[:].reshape(int(s / 3), 3)

                        patchGhost_part = Particles(
                            icells=np.asarray(dwpatch.data.iCell),
                            deltas=np.asarray(dwpatch.data.delta),
                            v=v,
                            weights=np.asarray(dwpatch.data.weight),
                            charges=np.asarray(dwpatch.data.charge),
                        )

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

    return PatchHierarchy(patch_levels, domain_box, time=simulator.currentTime())
