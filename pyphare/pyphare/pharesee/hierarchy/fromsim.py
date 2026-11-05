import numpy as np

from ...core import gridlayout  #  import GridLayout, HybridGridLayoutFor
from ...core.box import Box
from ...core.phare_utilities import refinement_ratio
from ..particles import LiveParticles
from .hierarchy import PatchHierarchy
from .hierarchy_utils import field_qties, isFieldQty, quantidic
from .patch import Patch
from .patchdata import FieldData, LiveParticleData
from .patchlevel import PatchLevel


def make_layout_for(simulator, patch_data, qty, dl, model):
    box = Box(patch_data.lower, patch_data.upper)
    origin = patch_data.origin
    # levels can have different models, see DataWrangler.modelsPerLevel
    if model == "MHD":
        return gridlayout.MHDGridLayoutFor(
            box, origin, dl, simulator.simulation.reconstruction
        )
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

    from pyphare import cpp

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
                patch_datas = getters[qty](pop)
            else:
                patch_datas = getters[qty]()

            if sync:
                patch_datas = dw.sync(patch_datas)
                if cpp.mpi_rank() > 0:
                    continue

            for patch_data in patch_datas:
                layout = make_layout_for(
                    simulator, patch_data, qty, lvl_cell_width, dw.modelsPerLevel[ilvl]
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
                layout = make_layout_for(
                    simulator, patch_data, qty, lvl_cell_width, dw.modelsPerLevel[ilvl]
                )
                live = LiveParticles(patch_data.data)

                patches[ilvl].append(Patch({qty: LiveParticleData(layout, live, pop)}))

        else:
            raise ValueError(f"{qty} is not a valid quantity")

        patch_levels[ilvl] = PatchLevel(ilvl, patches[ilvl])

    if hier:
        for lvl_nbr, level in hier.levels(hier.times()[0]).items():
            new_level = patch_levels[lvl_nbr]
            for ip, patch in enumerate(level.patches):
                patch.patch_datas = {**patch.patch_datas, **new_level[ip].patch_datas}

        return hier
    return PatchHierarchy({simulator.currentTime(): patch_levels}, domain_box)
