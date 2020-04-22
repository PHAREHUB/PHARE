from phare.pp.diagnostics import Diagnostic, Patch, _ParticlePatchData, Particles
from .overlap import Overlap, calculate_1d as overlap_calculate_1d
from .particle_overlap import ParticleOverlap, get_ghost_patch
from .periodic_overlap import intralevel as periodic_intralevel_1d

from phare.pp.diagnostic.patch_data import particlesForDiags


class PatchGhostParticleOverlap(ParticleOverlap):
    """
    Class containing information about an actual PatchGhost/Domain overlap on the same level, made from overlap_calculate_1d

    Parameters:
    -------------------

    lowerXPatch : patch with lower X origin
    upperXPatch : patch with upper X origin
    data_name   : Unused for particles
    nGhosts     : number of ghosts associated with the Particles datasets
    sizes       : size of the overlap box
    particles   : Particles object containing all particle types (domain, patchghost, lvlGhost) for the whole hierarchy
    periodic    : bool, true if overlap detected from periodic_overlap check - needed for iCell deduction/comparison

    """

    def __init__(
        self, lowerXPatch, upperXPatch, data_name, nGhosts, sizes, particles, periodic
    ):
        ParticleOverlap.__init__(
            self, lowerXPatch, upperXPatch, data_name, nGhosts, sizes
        )

        self.domainPatch = lowerXPatch
        self.ghostPatch = get_ghost_patch(
            particles, upperXPatch, PatchGhostParticleOverlap
        )
        self.particles = particles
        self.periodic = periodic
        self.is_mirror = False

        import copy  # copy and swap domain/ghost patches

        self.mirror = copy.copy(self)
        self.mirror.domainPatch = upperXPatch
        self.mirror.ghostPatch = get_ghost_patch(
            particles, lowerXPatch, PatchGhostParticleOverlap
        )
        self.mirror.is_mirror = True
        self.mirror.mirror = self


def getPatchGhostOverlaps(diags):
    assert isinstance(diags, dict) and _ParticlePatchData.__name__ in diags

    # Merge particle diagnostics of same population into single Particles object
    particlesList = particlesForDiags(diags[_ParticlePatchData.__name__])

    return _intralevel_1d(particlesList) + _periodic_1d(particlesList)


# no patchghost patches = no patchghost overlap
def _level_has_patchghost(particles, lvlNbr):
    return len(particles.pGhostDiag.levels[lvlNbr].patches.keys())


def _calculate_1d(particles, refPatch, cmpPatch, periodic=False):
    assert refPatch.patch_level().diag().dim == 1
    assert isinstance(particles, Particles)
    assert all([isinstance(x, Patch) for x in [refPatch, cmpPatch]])

    if _level_has_patchghost(particles, refPatch.patch_level().lvlNbr):
        return overlap_calculate_1d(
            PatchGhostParticleOverlap,
            refPatch,
            cmpPatch,
            "particles",
            particles=particles,
            periodic=periodic,
        )
    return []


def _periodic_1d(particlesList):
    overlaps = []
    for particles in particlesList:
        diag = particles.domainDiag
        for lvlNbr, patch_level in diag.levels.items():
            periodic_overlaps = periodic_intralevel_1d(particles, patch_level)
            data_names = list(periodic_overlaps.keys())
            if len(data_names):
                # we only need one out of the five, hence [0]
                minX, maxX = periodic_overlaps[data_names[0]]
                if len(minX) and len(maxX):
                    overlaps += _calculate_1d(
                        particles, minX[0], maxX[0], periodic=True
                    )
    return overlaps


def _intralevel_1d(particlesList):
    overlaps = []
    for particles in particlesList:
        for lvlNbr, patch_level in particles.domainDiag.levels.items():
            patches = patch_level.patches_list()
            for i, refPatch in enumerate(patches):
                for cmpPatch in patches[i + 1 :]:
                    overlaps += _calculate_1d(particles, refPatch, cmpPatch)
    return overlaps
