from phare.pp.diagnostics import Diagnostic, Patch, _ParticlePatchData, Particles
from .particle_overlap import ParticleOverlap, get_ghost_patch
from .patch_ghost_particle_overlap import PatchGhostParticleOverlap
from .overlap import Overlap

from phare.pp.diagnostic.patch_data import particlesForDiags
from phare.pp.diagnostic.patch import REFINEMENT_RATIO


class LevelGhostParticleOverlap(ParticleOverlap):
    """
    Class containing information about an actual CoarseLevelGhost/FinePatchGhost overlap

    Parameters:
    -------------------

    coarseDomainPatch   : patch containing coarse domain data
    fineLevelGhostPatch : patch containing fine levelGhost data
    data_name           : Unused for particles
    nGhosts             : number of ghosts associated with the Particles datasets
    sizes               : size of the overlap box
    particles           : Particles object containing all particle types (domain, patchghost, lvlGhost) for the whole hierarchy

    """

    def __init__(
        self,
        coarseDomainPatch,
        fineLevelGhostPatch,
        data_name,
        nGhosts,
        firstCommonFineiCell,
        particles,
    ):
        ParticleOverlap.__init__(
            self,
            coarseDomainPatch,
            fineLevelGhostPatch,
            data_name,
            nGhosts,
            firstCommonFineiCell,
        )
        self.coarseDomainPatch = coarseDomainPatch
        self.coarsePatchGhostPatch = get_ghost_patch(
            particles, coarseDomainPatch, PatchGhostParticleOverlap
        )
        self.fineLevelGhostPatch = fineLevelGhostPatch
        self.firstCommonFineiCell = firstCommonFineiCell
        self.particles = particles


def getLevelGhostOverlaps(diags):
    assert isinstance(diags, dict) and _ParticlePatchData.__name__ in diags

    # Merge particle diagnostics of same population into single Particles object
    particlesList = particlesForDiags(diags[_ParticlePatchData.__name__])

    import itertools as it

    # Flatten list of lists into single list
    return list(it.chain(*[_levelGhost(particles) for particles in particlesList]))


def _levelGhost(particles):
    overlaps = []
    for fineLevel in list(particles.domainDiag.levels.keys())[1:]:
        coarseLevel = fineLevel - 1
        coarseDomainPatches = particles.domainDiag.levels[coarseLevel].patches_list()
        fineLevelGhostPatches = particles.lGhostDiag.levels[fineLevel].patches_list()
        for coarseDomain in coarseDomainPatches:
            for fineLevelGhost in fineLevelGhostPatches:
                overlaps += _intersect_coarser_1d(
                    coarseDomain, fineLevelGhost, particles
                )
    return overlaps


def _intersect_coarser_1d(coarseDomain, fineLevelGhost, particles):
    """Particle ghost width is never > 2, as such a fine/coarse overlap will never
        include more than 1 coarsePatch/cell"""
    assert coarseDomain.patch_level.diag.dim == 1
    assert all([isinstance(x, Patch) for x in [coarseDomain, fineLevelGhost]])
    assert coarseDomain.patch_level.lvlNbr == fineLevelGhost.patch_level.lvlNbr - 1

    import numpy as np

    direction = "x"
    nGhosts = coarseDomain.patch_data.nGhosts()
    coarse_minX, coarse_maxX = coarseDomain.min_max_coords(direction)
    fine_minX, fine_maxX = fineLevelGhost.min_max_coords(direction)
    uniq_fineLevelGhostICell = np.unique(fineLevelGhost.data("iCell")).tolist()

    def _check(fine_X, shift=0):
        fineICell = (
            fineLevelGhost.patch_level.position_to_amr_index(fine_X, direction) + shift
        )
        has_overlap = coarse_minX <= fine_X <= coarse_maxX
        levelOverlapDataExists = fineICell in uniq_fineLevelGhostICell
        if has_overlap and levelOverlapDataExists:
            return LevelGhostParticleOverlap(
                coarseDomain,
                fineLevelGhost,
                "particles",
                nGhosts,
                [fineICell],
                particles=particles,
            )

    return [
        x for x in [_check(fine_minX, -nGhosts), _check(fine_maxX)] if x is not None
    ]
