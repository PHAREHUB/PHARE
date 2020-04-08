from phare.pp.diagnostics import Diagnostic, Patch, _ParticlePatchData, Particles
from .particle_overlap import ParticleOverlap, get_ghost_patch
from .particle_patch_overlap import DomainParticleOverlap
from .overlap import Overlap

from phare.pp.diagnostic.patch_data import particlesForDiags
from phare.pp.diagnostic.patch import REFINEMENT_RATIO


class LevelParticleOverlap(ParticleOverlap):
    def __init__(
        self,
        coarseDomain,
        fineLevelGhost,
        data_name,
        nGhosts,
        sizes,
        coarsePatchGhost,
    ):
        ParticleOverlap.__init__(
            self, coarseDomain, fineLevelGhost, data_name, nGhosts, sizes
        )
        self.coarseDomain = coarseDomain
        self.coarsePatchGhost = coarsePatchGhost
        self.fineLevelGhost = fineLevelGhost


def getLevelGhostOverlaps(diags):
    if isinstance(diags, dict) and _ParticlePatchData.__name__ in diags:
        diags = diags[_ParticlePatchData.__name__]

    # Merge diagnotics of same population into single Particles object
    if all([isinstance(diag, Diagnostic) for diag in diags]):
        diags = particlesForDiags(diags)

    assert all([isinstance(diag, Particles) for diag in diags])

    return _levelGhost(LevelParticleOverlap, diags)


def _levelGhost(OverlapType, diags):
    overlaps = []
    for particles in diags:
        for refineLevel in range(1, len(list(particles.domainDiag.levels.keys()))):
            coarseDomainPatches = list(
                particles.domainDiag.levels[refineLevel - 1].patches.values()
            )
            for coarseDomain in coarseDomainPatches:
                fineLevelGhostPatches = list(
                    particles.lGhostDiag.levels[refineLevel].patches.values()
                )
                for fineLevelGhost in fineLevelGhostPatches:
                    overlaps.extend(
                        _intersect_coarser_1d(
                            OverlapType,
                            coarseDomain,
                            get_ghost_patch(
                                particles, coarseDomain, DomainParticleOverlap
                            ),
                            fineLevelGhost,
                        )
                    )

    return overlaps


def _intersect_coarser_1d(OverlapType, coarseDomain, coarsePatchGhost, fineLevelGhost):
    assert all(
        [isinstance(x, Patch) for x in [coarseDomain, coarsePatchGhost, fineLevelGhost]]
    )

    return _calculate_coarser_1d(
        OverlapType, coarseDomain, coarsePatchGhost, fineLevelGhost, "particles"
    )


def _calculate_coarser_1d(
    OverlapType, coarseDomain, coarsePatchGhost, fineLevelGhost, data_name
):
    """Particle ghost width is never > 2, as such a fine/coarse overlap will never
        include more than 1 coarsePatch/cell"""

    assert issubclass(OverlapType, Overlap)
    assert coarseDomain.patch_level.lvlNbr == fineLevelGhost.patch_level.lvlNbr - 1
    import numpy as np

    direction = "x"

    coarse_minX = coarseDomain.min_coord(direction)
    coarse_maxX = coarseDomain.max_coord(direction)
    refine_minX = fineLevelGhost.min_coord(direction)
    refine_maxX = fineLevelGhost.max_coord(direction)

    uniq_fineLevelGhostICell = np.unique(
        fineLevelGhost.patch_data.data("iCell")
    ).tolist()

    nGhosts = coarseDomain.patch_data.nGhosts(data_name)

    def coarse_refine_overlap(refine_start_idx):
        return OverlapType(
            coarseDomain,
            fineLevelGhost,
            data_name,
            nGhosts,
            [refine_start_idx],
            coarsePatchGhost=coarsePatchGhost,
        )

    def check_lower():
        # iCell check if data exists, otherwise patchGhost overlap assumed
        fineICell = (
            fineLevelGhost.patch_level.position_to_amr_ndex(refine_minX, direction) - 1
        )
        if (
            coarse_minX <= refine_minX <= coarse_maxX
            and fineICell in uniq_fineLevelGhostICell
        ):
            return coarse_refine_overlap(
                fineLevelGhost.patch_level.position_to_amr_ndex(refine_minX, direction)
                - nGhosts
            )

    def check_upper():
        # iCell check if data exists, otherwise patchGhost overlap assumed
        fineICell = fineLevelGhost.patch_level.position_to_amr_ndex(
            refine_maxX, direction
        )
        if (
            coarse_maxX >= refine_maxX >= coarse_minX
            and fineICell in uniq_fineLevelGhostICell
        ):
            return coarse_refine_overlap(
                fineLevelGhost.patch_level.position_to_amr_ndex(refine_maxX, direction)
            )

    return [x for x in [check_lower(), check_upper()] if x is not None]
