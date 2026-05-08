from pyphare import cpp
from .hierarchy_utils import isFieldQty, field_qties, quantidic, refinement_ratio
from .patchdata import FieldData, ParticleData
from .patch import Patch
from .patchlevel import PatchLevel
from .hierarchy import PatchHierarchy
from ..particles import Particles
from ...core.gridlayout import GridLayout
from ...core.box import Box

import numpy as np


def patch_gridlayout(patch, dl, interp_order):
    return GridLayout(
        Box(patch.lower, patch.upper), patch.origin, dl, interp_order=interp_order
    )


def hierarchy_from_sim(simulator, qty, pop="", hier=None, sync=False):
    """
    sync==True will copy all data to rank 0!
    leaving other ranks with an empty hierarch!
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
                layout = patch_gridlayout(
                    patch_data, lvl_cell_width, simulator.interp_order()
                )
                patches[ilvl].append(
                    Patch({qty: FieldData(layout, field_qties[qty], patch_data.data)})
                )

        elif qty == "particles":  # domain only!
            if sync:
                raise ValueError("sync not supported for particles")
            if pop == "":
                raise ValueError("must specify pop argument for particles")

            for patch_data in getters[qty](pop):
                layout = patch_gridlayout(
                    patch_data, lvl_cell_width, simulator.interp_order()
                )

                patches[ilvl].append(  # ParticleData is SoA COPY!
                    Patch({qty: FieldData(layout, "tags", patch_data.data)})
                )

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
