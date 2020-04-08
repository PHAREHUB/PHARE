from phare.pp.diagnostics import Diagnostic, Patch, _ParticlePatchData, Particles
from .overlap import Overlap, _calculate_1d as overlap_calculate_1d
from .particle_overlap import ParticleOverlap, get_ghost_patch
from .periodic_overlap import intralevel as periodic_intralevel

from phare.pp.diagnostic.patch_data import particlesForDiags


class DomainParticleOverlap(ParticleOverlap):
    def __init__(
        self, domainPatch, ghostPatch, data_name, nGhosts, sizes, particles, periodic
    ):
        ParticleOverlap.__init__(
            self, domainPatch, ghostPatch, data_name, nGhosts, sizes
        )
        self.domain = domainPatch
        self.ghost = ghostPatch
        self.particles = particles
        self.periodic = periodic

        import copy  # copy and swap domain/ghost patches

        self.mirror = copy.copy(self)
        self.mirror.domain = particles.domainDiag.levels[
            domainPatch.patch_level.lvlNbr
        ].patches[ghostPatch.id]
        self.mirror.ghost = get_ghost_patch(
            particles, self.domain, DomainParticleOverlap
        )
        self.mirror.mirror = self


def getPatchGhostOverlaps(diags):
    if isinstance(diags, dict) and _ParticlePatchData.__name__ in diags:
        diags = diags[_ParticlePatchData.__name__]

    # Merge diagnotics of same population into single Particles object
    if all([isinstance(diag, Diagnostic) for diag in diags]):
        diags = particlesForDiags(diags)

    assert all([isinstance(diag, Particles) for diag in diags])

    return _intralevel(DomainParticleOverlap, diags) + _periodic(
        DomainParticleOverlap, diags
    )


def _calculate_1d(OverlapType, particles, patch0, patch1, periodic=False):
    assert isinstance(particles, Particles)
    assert all([isinstance(x, Patch) for x in [patch0, patch1]])

    return overlap_calculate_1d(
        OverlapType,
        patch0,
        get_ghost_patch(particles, patch1, OverlapType),
        "particles",
        particles=particles,
        periodic=periodic,
    )


def _periodic(OverlapType, diags):
    overlaps = []
    for particles in diags:
        diag = particles.domainDiag
        for patch_level_ids, patch_level in diag.levels.items():
            periodic_overlaps = periodic_intralevel(particles, patch_level)

            # we only need one out of the five, v is arbitrary (and short ;))
            if "v" in periodic_overlaps:
                minX, maxX = periodic_overlaps["v"]
                if len(minX) and len(maxX):
                    overlaps += _calculate_1d(
                        OverlapType, particles, minX[0], maxX[0], periodic=True
                    )

    return overlaps


def _intralevel(OverlapType, diags):
    overlaps = []
    for particles in diags:
        for patch_level_ids, patch_level in particles.domainDiag.levels.items():
            patches = list(patch_level.patches.values())
            for i, patch0 in enumerate(patches):
                for j in range(i + 1, len(patches)):
                    overlaps += _calculate_1d(
                        OverlapType, particles, patch0, patches[j]
                    )

    return overlaps
