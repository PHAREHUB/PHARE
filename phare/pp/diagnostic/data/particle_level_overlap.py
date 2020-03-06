from phare.pp.diagnostics import Patch, _Particle, Particles
from .particle_overlap import ParticleOverlap
from .particle_patch_overlap import DomainParticleOverlap
from .overlap import Overlap


class LevelParticleOverlap(ParticleOverlap):
    def __init__(self, patch0, patch1, dataset_key, nGhosts, offsets, sizes):
        ParticleOverlap.__init__(
            self, patch0, patch1, dataset_key, nGhosts, offsets, sizes
        )

    @classmethod
    def get(clazz, diags):
        if isinstance(diags, dict) and _Particle.__name__ in diags:
            diags = diags[_Particle.__name__]

        return clazz.levelGhost(diags)

    @classmethod
    def levelGhost(clazz, diags):

        overlaps = []

        for particles in diags:
            for refineLevel in range(1, len(list(particles.domain.levels.keys()))):
                for domainPatch in particles.domain.levels[refineLevel - 1].patches:
                    for levelGhostPatch in particles.lGhost.levels[refineLevel].patches:
                        overlaps.extend(
                            clazz.intersect_coarser_1d(
                                domainPatch,
                                ParticleOverlap.get_ghost_patch(
                                    particles, domainPatch, DomainParticleOverlap
                                ),
                                levelGhostPatch,
                            )
                        )
        return overlaps

    @classmethod
    def intersect_coarser_1d(clazz, coarseDomain, coarsePatchGhost, fineLevelGhost):
        assert all(
            [
                isinstance(x, Patch)
                for x in [coarseDomain, coarsePatchGhost, fineLevelGhost]
            ]
        )

        overlaps = clazz.calculate_coarser_1d(
            coarseDomain,
            coarsePatchGhost,
            fineLevelGhost,
            "particles",
            LevelParticleOverlap,
        )

        for o in overlaps:
            o.__dict__["fineLevelGhost"] = o.__dict__.pop("p1")
            o.coarseDomain, o.coarsePatchGhost = o.p0

        return overlaps

    @staticmethod
    def calculate_coarser_1d(
        coarseDomain, coarsePatchGhost, fineLevelGhost, dataset_key: str, clazz
    ):
        """Particle ghost width is never > 2, as such a fine/coarse overlap will never
            include more than 1 coarsePatch"""

        assert issubclass(clazz, Overlap)
        assert coarseDomain.patch_level.idx == fineLevelGhost.patch_level.idx - 1

        direction = "x"
        course_minX, course_maxX = (
            coarseDomain.min_coord(direction),
            coarseDomain.max_coord(direction),
        )
        refine_minX, refine_maxX = (
            fineLevelGhost.min_coord(direction),
            fineLevelGhost.max_coord(direction),
        )
        nGhosts = coarseDomain.dtype.nGhosts(dataset_key)

        def coarse_refine_overlap(refine_start_idx):
            return clazz(
                [coarseDomain, coarsePatchGhost],
                fineLevelGhost,
                dataset_key,
                nGhosts,
                [refine_start_idx],
                [nGhosts],
            )

        def check_lower():
            if course_minX <= refine_minX <= course_maxX:
                return coarse_refine_overlap(
                    fineLevelGhost.patch_level.position_to_index(refine_minX, direction)
                    - nGhosts
                )

        def check_upper():
            if course_maxX >= refine_maxX >= course_minX:
                return coarse_refine_overlap(
                    fineLevelGhost.patch_level.position_to_index(refine_maxX, direction)
                )

        return [x for x in [check_lower(), check_upper()] if x is not None]
