from phare.pp.diagnostics import Patch, _Particle, Particles
from .overlap import Overlap, _calculate_1d as overlap_calculate_1d
from .particle_overlap import ParticleOverlap, get_ghost_patch
from .periodic_overlap import intralevel as periodic_intralevel


class DomainParticleOverlap(ParticleOverlap):
    def __init__(self, patch0, patch1, dataset_key, nGhosts, sizes):
        ParticleOverlap.__init__(self, patch0, patch1, dataset_key, nGhosts, sizes)


def getPatchGhostOverlaps(diags):
    if isinstance(diags, dict):
        if _Particle.__name__ in diags:
            diags = diags[_Particle.__name__]

    return _intralevel(DomainParticleOverlap, diags) + _periodic(
        DomainParticleOverlap, diags
    )


def _calculate_1d(clazz, particles, patch0, patch1, coords=[]):
    assert isinstance(particles, Particles)
    assert all([isinstance(x, Patch) for x in [patch0, patch1]])

    direction = "x"
    overlaps = []
    for p in [(patch0, patch1), (patch1, patch0)]:
        overlaps += overlap_calculate_1d(
            clazz, p[0], get_ghost_patch(particles, p[1], clazz), "particles"
        )

    if len(coords) is 0:  # periodicity overlap icell matching
        _, upper = sorted([patch0, patch1], key=lambda x: x.min_coord(direction))
        coords = [
            upper.patch_level.position_to_index(idx, direction) for idx in upper.origin
        ]

    for o in overlaps:
        o.origin = coords
        if patch0 == o.p0:
            o.domain = o.p0
            o.__dict__["ghost"] = o.__dict__.pop("p1")
        else:
            o.ghost = o.p0
            o.__dict__["domain"] = o.__dict__.pop("p1")
    return overlaps


def _periodic(clazz, diags):
    overlaps = []
    for particles in diags:
        diag = particles.domain
        for patch_level_ids, patch_level in diag.levels.items():
            # we only need one out of the five, v is arbitrary (and short ;))
            minX, maxX = periodic_intralevel(clazz, particles, patch_level.patches)["v"]
            if len(minX) and len(maxX):
                overlaps += _calculate_1d(
                    clazz, particles, minX[0], maxX[0], diag.sim.origin
                )
    return overlaps


def _intralevel(clazz, diags):
    overlaps = []
    for particles in diags:
        for patch_level_ids, patch_level in particles.domain.levels.items():
            patches = patch_level.patches
            for i, patch0 in enumerate(patches):
                for j in range(i + 1, len(patches)):
                    overlaps += _calculate_1d(clazz, particles, patch0, patches[j])
    return overlaps
