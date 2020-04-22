import numpy as np, math
from phare.pp.diagnostics import Patch, _ParticlePatchData, Particles
from .overlap import Overlap
from phare.data.wrangler import safe_round


class ParticleOverlap(Overlap):
    def __init__(self, patch0, patch1, data_name, nGhosts, sizes):
        Overlap.__init__(self, patch0, patch1, data_name, nGhosts, sizes)


def getParticleOverlapsFrom(diags):
    from .level_ghost_particle_overlap import getLevelGhostOverlaps
    from .patch_ghost_particle_overlap import getPatchGhostOverlaps

    return getLevelGhostOverlaps(diags) + getPatchGhostOverlaps(diags)


def get_ghost_patch(particles, refPatch, ghostType):
    from .level_ghost_particle_overlap import LevelGhostParticleOverlap
    from .patch_ghost_particle_overlap import PatchGhostParticleOverlap

    assert (
        ghostType == PatchGhostParticleOverlap or ghostType == LevelGhostParticleOverlap
    )

    gDiag = particles.lGhostDiag
    if ghostType is PatchGhostParticleOverlap:
        gDiag = particles.pGhostDiag

    return gDiag.levels[refPatch.patch_level().lvlNbr].patches[refPatch.id]


class ParticleOverlapComparator:
    """
    Class to compare PatchGhost with overlapping same level PatchDomain
          or compare LevelGhost with overlapping coarser level PatchDomain

    Parameters:
    -------------------

    patch     : Patch object / or C++ ContiguousParticles class representing a Patch object
    indices   : Subset of Particle indices to use for Comparision with other "self.patch"
                  i.e. self.deltas[self.indices[0]] == that.deltas[that.indices[0]]

    """

    def __init__(self, patch, indices):

        assert isinstance(patch, Patch) or type(patch).__name__.startswith(
            "ContiguousParticles"
        )
        if isinstance(patch, Patch):
            self.dim = len(patch.origin)
        else:  # or SoA ContigousParticles
            self.dim = safe_round(len(patch.iCell), len(patch.weight))
        self.patch = patch
        self.indices = indices

    def sortedIndices(self):  # ONLY WORKS FOR 1D!
        """returns list[indices] sorted on iCell + delta"""
        sorted_zip = sorted(
            zip(
                self._dataset("iCell")[self.indices]
                + self._dataset("delta")[self.indices],
                self.indices,
            )
        )
        return [index for _, index in sorted_zip]

    def cmp(self, that, data_name, sortedIndices, dim):
        ref_datasets = [x._dataset(data_name) for x in [self, that]]
        cmp_datasets = []
        for i, dataset in enumerate(ref_datasets):
            cmp_datasets.append(
                dataset.reshape(safe_round(len(dataset), dim), dim)[sortedIndices[i]]
            )
        return np.array_equiv(*cmp_datasets)

    def __eq__(self, that):
        assert type(self) is type(that)
        assert len(self.indices) == len(that.indices)

        sortedIndices = [self.sortedIndices(), that.sortedIndices()]
        return len(sortedIndices[0]) == len(sortedIndices[1]) and all(
            [
                self.cmp(that, "v", sortedIndices, 3),
                self.cmp(that, "delta", sortedIndices, self.dim),
                [
                    self.cmp(that, data_name, sortedIndices, 1)
                    for data_name in ["weight", "charge"]
                ],
            ]
        )

    def _dataset(self, data_name):
        if isinstance(self.patch, Patch):
            return np.asarray(self.patch.data(data_name))
        # else ContigousParticles
        return np.asarray(getattr(self.patch, data_name))
