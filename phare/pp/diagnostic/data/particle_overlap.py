import numpy as np, math
from phare.pp.diagnostics import Patch, _ParticlePatchData, Particles
from .overlap import Overlap


class ParticleOverlap(Overlap):
    def __init__(self, patch0, patch1, data_name, nGhosts, sizes):
        Overlap.__init__(self, patch0, patch1, data_name, nGhosts, sizes)


def getParticleOverlapsFrom(diags):
    from .particle_level_overlap import getLevelGhostOverlaps
    from .particle_patch_overlap import getPatchGhostOverlaps

    return getLevelGhostOverlaps(diags) + getPatchGhostOverlaps(diags)


def get_ghost_patch(particles, patch, ghostType):
    from .particle_level_overlap import LevelParticleOverlap
    from .particle_patch_overlap import DomainParticleOverlap

    assert ghostType == DomainParticleOverlap or ghostType == LevelParticleOverlap

    gDiag = particles.pGhostDiag if ghostType is DomainParticleOverlap else particles.lGhostDiag

    return gDiag.levels[patch.patch_level.lvlNbr].patches[patch.id]


class ParticleOverlapComparator:
    def __init__(self, patch, indices):
        assert isinstance(patch, Patch) or type(patch).__name__.startswith(
            "ContiguousParticles"
        )
        if isinstance(patch, Patch):
            self.dim = len(patch.origin)
        else:  # or SoA ContigousParticles
            self.dim = math.floor(len(patch.iCell) / len(patch.weight))
        self.patch = patch
        self.indices = indices

    def sort(self):  # ONLY WORKS FOR 1D!
        """returns list[tuple()] sorted on tuple[0])
            tuple[0] = particle["icell"] + particle["delta"]
            tuple[1] = position of icell/delta in contigous arrays/particle index
        """

        icells = self._dataset("iCell")
        delta = [
            (particleDelta + icells[self.indices[i]], self.indices[i])
            for i, particleDelta in enumerate(
                # select deltas from indices
                list(map(self._dataset("delta").__getitem__, self.indices))
            )
        ]
        return sorted(delta, key=lambda x: x[0])

    def cmp(self, that, data_name, sortedDeltaTuple0, sortedDeltaTuple1, dim):
        dataset0, dataset1 = [x._dataset(data_name) for x in [self, that]]

        return all(
            [
                np.array_equiv(
                    dataset0[
                        sortedDeltaTuple0[i][1] * dim : sortedDeltaTuple0[i][1] * dim
                        + dim
                    ],
                    dataset1[
                        sortedDeltaTuple1[i][1] * dim : sortedDeltaTuple1[i][1] * dim
                        + dim
                    ],
                )
                for i in range(len(sortedDeltaTuple0))
            ]
        )

    def __eq__(self, that):
        assert type(self) is type(that)

        sortedDeltaTuple0, sortedDeltaTuple1 = self.sort(), that.sort()
        return len(sortedDeltaTuple0) == len(sortedDeltaTuple1) and all(
            [
                self.cmp(that, "v", sortedDeltaTuple0, sortedDeltaTuple1, 3),
                self.cmp(that, "delta", sortedDeltaTuple0, sortedDeltaTuple1, self.dim),
                [
                    self.cmp(that, data_name, sortedDeltaTuple0, sortedDeltaTuple1, 1)
                    for data_name in ["weight", "charge"]
                ],
            ]
        )

    def _dataset(self, data_name):
        if isinstance(self.patch, Patch):
            return self.patch.patch_data.data(data_name)
        # else ContigousParticles
        return getattr(self.patch, data_name)
