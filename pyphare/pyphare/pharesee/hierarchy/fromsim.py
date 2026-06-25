from pyphare import cpp
from pyphare.pharesee.hierarchy.hierarchy_utils import (
    isFieldQty,
    field_qties,
    quantidic,
    refinement_ratio,
)
from pyphare.pharesee.hierarchy.patchdata import FieldData, LiveParticleData
from pyphare.pharesee.particles import LiveParticles
from pyphare.pharesee.hierarchy.patch import Patch
from pyphare.pharesee.hierarchy.patchlevel import PatchLevel
from pyphare.pharesee.hierarchy.hierarchy import PatchHierarchy
from pyphare.core import gridlayout
from pyphare.core.box import Box

import numpy as np


def make_layout_for(simulator, patch_data, qty, dl):
    box = Box(patch_data.lower, patch_data.upper)
    origin = patch_data.origin
    sim = simulator.simulation
    if "MHDModel" in getattr(sim, "model_options", []):
        return gridlayout.MHDGridLayoutFor(box, origin, dl, sim.reconstruction)
    return gridlayout.HybridGridLayoutFor(
        box,
        origin,
        dl,
        simulator.interp_order(),
        is_particle_layout=(qty == "particles"),
    )


def hierarchy_from_sim(simulator, qty, pop="", hier=None, sync=False):
    """
    sync == False == live simulation data
    sync == True  == copy all data to rank 0 - eg. plotting
    """

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
            patch_datas = getters[qty]()

            if sync:
                patch_datas = dw.sync(patch_datas)
                if cpp.mpi_rank() > 0:
                    continue

            for patch_data in patch_datas:
                layout = make_layout_for(simulator, patch_data, qty, lvl_cell_width)

                patches[ilvl].append(
                    Patch({qty: FieldData(layout, field_qties[qty], patch_data.data)})
                )

        elif qty == "particles":  # domain only!
            if sync:
                raise ValueError("sync not supported for particles")
            if pop == "":
                raise ValueError("must specify pop argument for particles")

            for patch_data in getters[qty](pop):
                layout = make_layout_for(simulator, patch_data, qty, lvl_cell_width)
                live = LiveParticles(patch_data.data)

                patches[ilvl].append(Patch({qty: LiveParticleData(layout, live, pop)}))

        else:
            raise ValueError("{} is not a valid quantity".format(qty))

        patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

    if hier:
        for lvl_nbr, level in hier.levels(hier.times()[0]).items():
            new_level = patch_levels[lvl_nbr]
            for ip, patch in enumerate(level.patches):
                patch.patch_datas = {**patch.patch_datas, **new_level[ip].patch_datas}

        return hier
    return PatchHierarchy(patch_levels, domain_box, time=simulator.currentTime())
