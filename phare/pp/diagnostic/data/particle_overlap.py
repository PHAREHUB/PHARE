import numpy as np, math
from phare.pp.diagnostics import Patch, _Particle, Particles
from .overlap import Overlap


class ParticleOverlap(Overlap):
    def __init__(self, p0: Patch, p1: Patch, dataset_key, nGhosts, sizes):
        Overlap.__init__(self, p0, p1, dataset_key, nGhosts, sizes)


def getParticleOverlapsFrom(diags):
    from .particle_level_overlap import getLevelGhostOverlaps
    from .particle_patch_overlap import getPatchGhostOverlaps

    return getLevelGhostOverlaps(diags) + getPatchGhostOverlaps(diags)


def get_ghost_patch(particles: Particles, p0: Patch, gType):
    from .particle_level_overlap import LevelParticleOverlap
    from .particle_patch_overlap import DomainParticleOverlap

    assert gType == DomainParticleOverlap or gType == LevelParticleOverlap

    gDiag = particles.pGhost if gType is DomainParticleOverlap else particles.lGhost
    return gDiag.levels[p0.patch_level.idx].patchDict[p0.id]


class ParticleOverlapComparator:
    def __init__(self, p0, indices):
        assert isinstance(p0, Patch) or type(p0).__name__.startswith(
            "ContiguousParticles"
        )
        if isinstance(p0, Patch):
            self.dim = len(p0.origin)
        else:  # or SoA ContigousParticles
            self.dim = math.floor(len(p0.iCell) / len(p0.weight))
        self.p0 = p0
        self.indices = indices

    def sort(self):  # ONLY WORKS FOR 1D!
        icells = self._get("iCell")
        delta = [
            (v + icells[self.indices[i]], self.indices[i])
            for i, v in enumerate(  # select deltas from indices
                list(map(self._get("delta").__getitem__, self.indices))
            )
        ]
        return sorted(delta, key=lambda x: x[0])

    def cmp(self, that, s, d0, d1, dim):
        a0, a1 = [x._get(s) for x in [self, that]]
        return all(
            [
                np.array_equiv(
                    a0[d0[i][1] * dim : d0[i][1] * dim + dim],
                    a1[d1[i][1] * dim : d1[i][1] * dim + dim],
                )
                for i in range(len(d0))
            ]
        )

    def __eq__(self, that):
        assert type(self) is type(that)

        d0, d1 = self.sort(), that.sort()
        return len(d0) == len(d1) and all(
            [
                self.cmp(that, "v", d0, d1, 3),
                self.cmp(that, "delta", d0, d1, self.dim),
                [self.cmp(that, s, d0, d1, 1) for s in ["weight", "charge"]],
            ]
        )

    def _get(self, attr):
        if isinstance(self.p0, Patch):
            return self.p0.dtype.get()[attr]
        return getattr(self.p0, attr)
