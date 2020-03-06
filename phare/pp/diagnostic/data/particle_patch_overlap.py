from phare.pp.diagnostics import Patch, _Particle, Particles
from .periodic_overlap import PeriodicOverlap
from .particle_overlap import ParticleOverlap
from .overlap import Overlap


class DomainParticleOverlap(ParticleOverlap):
    def __init__(self, patch0, patch1, dataset_key, nGhosts, offsets, sizes):
        ParticleOverlap.__init__(
            self, patch0, patch1, dataset_key, nGhosts, offsets, sizes
        )

    @classmethod
    def _calculate_1d(clazz, particles, patch0, patch1, coords=[]):
        assert isinstance(particles, Particles)
        assert all([isinstance(x, Patch) for x in [patch0, patch1]])

        direction = "x"
        overlaps = []
        for p in [(patch0, patch1), (patch1, patch0)]:
            overlaps += Overlap._calculate_1d(
                p[0],
                ParticleOverlap.get_ghost_patch(particles, p[1], clazz),
                "particles",
                clazz,
            )

        if len(coords) is 0:  # periodicity overlap icell matching
            _, upper = sorted([patch0, patch1], key=lambda x: x.min_coord(direction))
            coords = [
                upper.patch_level.position_to_index(idx, direction)
                for idx in upper.origin
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

    @classmethod
    def periodic(clazz, diags):
        overlaps = []
        for particles in diags:
            diag = particles.domain
            for patch_level_ids, patch_level in diag.levels.items():
                # we only need one out of the five, v is arbitrary (and short ;))
                minX, maxX = PeriodicOverlap.intralevel(particles, patch_level.patches)[
                    "v"
                ]
                if len(minX) and len(maxX):
                    overlaps += clazz._calculate_1d(
                        particles, minX[0], maxX[0], diag.sim.origin
                    )
        return overlaps

    @classmethod
    def intralevel(clazz, diags):
        overlaps = []
        for particles in diags:
            for patch_level_ids, patch_level in particles.domain.levels.items():
                patches = patch_level.patches
                for i, patch0 in enumerate(patches):
                    for j in range(i + 1, len(patches)):
                        overlaps += clazz._calculate_1d(particles, patch0, patches[j])
        return overlaps

    @classmethod
    def get(clazz, diags):
        if isinstance(diags, dict):
            if _Particle.__name__ in diags:
                diags = diags[_Particle.__name__]

        return clazz.intralevel(diags) + clazz.periodic(diags)
