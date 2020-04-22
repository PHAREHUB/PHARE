#!/usr/bin/env python3
#
# formatted with black

from tests.simulator import test_simulator as tst
from tests.simulator.py import InitValueValidation

from ddt import ddt, data

import unittest, os, phare.pharein as ph, numpy as np

from tests.diagnostic import dump_all_diags

from phare.pp.diagnostic.data.overlap import Overlap
from phare.pp.diagnostic.data.em_overlap import EMOverlap, getEMOverlapsFrom
from phare.pp.diagnostic.data.fluid_overlap import FluidOverlap, getFluidOverlapsFrom
from phare.pp.diagnostic.data.particle_overlap import (
    ParticleOverlap,
    ParticleOverlapComparator,
    getParticleOverlapsFrom,
)
from phare.pp.diagnostic.data.level_ghost_particle_overlap import (
    LevelGhostParticleOverlap,
)
from phare.pp.diagnostic.data.patch_ghost_particle_overlap import (
    PatchGhostParticleOverlap,
)

from phare.pp.diagnostics import (
    Diagnostics,
    Diagnostic,
    Patch,
    _EMPatchData,
    _FluidPatchData,
    _ParticlePatchData,
)


diag_out_dir = "phare_outputs/initializer"


@ddt
class OverlapValueValidation(InitValueValidation):
    min_interp = 1
    max_interp = 1

    def add_to_dict(dic):
        dic.update(InitValueValidation.diag_options(diag_out_dir))
        dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
        return dic

    valid1D = [
        add_to_dict(
            {
                "refinement_boxes": {
                    "L0": {
                        "B0": [
                            (5,),
                            (29,),
                        ],  # not touching overlap gap of 1 cell, 2 on fine level
                        "B1": [(31,), (55,)],
                    }
                },
            }
        ),
        add_to_dict(
            {
                "refinement_boxes": {
                    "L0": {
                        "B0": [
                            (5,),
                            (29,),
                        ],  # not touching overlap gap of 2 cell, 4 on fine level
                        "B1": [(32,), (55,)],
                    }
                },
            }
        ),
        add_to_dict(
            {
                "refinement_boxes": {
                    "L0": {"B0": [(0,), (9,)], "B1": [(55,), (64,)]},
                    "L1": {"B0": [(0,), (9,)], "B1": [(120,), (129,)]},
                    "L2": {"B0": [(0,), (9,)], "B1": [(250,), (259,)]},
                    "L3": {"B0": [(0,), (9,)], "B1": [(510,), (519,)]},
                    "L4": {"B0": [(0,), (9,)], "B1": [(1030,), (1039,)]},
                    "L5": {"B0": [(0,), (9,)], "B1": [(2070,), (2079,)]},
                    "L6": {"B0": [(0,), (9,)], "B1": [(4150,), (4159,)]},
                    "L7": {"B0": [(0,), (9,)], "B1": [(8310,), (8319,)]},
                },
                "max_nbr_levels": 9,
            }
        ),
        add_to_dict(
            {
                "refinement_boxes": {
                    "L0": {"B0": [(5,), (55,)]},
                    "L1": {"B0": [(20,), (100,)]},
                },
                "max_nbr_levels": 3,
            }
        ),
        add_to_dict({"refinement_boxes": {"L0": {"B0": [(5,), (54,)]}},}),
    ]

    @data(*valid1D)
    def test_1d(self, input):
        dim = 1
        test_id = self.id().split("_")[-1]  # first id == "test_1d_1"
        min_interp = OverlapValueValidation.min_interp
        max_interp = OverlapValueValidation.max_interp + 1
        for interp in range(min_interp, max_interp):
            diags = self.runAndDump(dim, interp, input)
            self._checkFields(getEMOverlapsFrom(diags))
            self._checkFields(getFluidOverlapsFrom(diags))
            self._checkParticles(getParticleOverlapsFrom(diags), dim, interp)

    def _checkFields(self, overlaps):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            patch0_data, patch1_data = overlap.shared_data()
            self.assertTrue(len(patch0_data) == len(patch1_data))
            if len(patch0_data) > 1:  # fluid is 1 value and can be 0
                zeros = np.zeros(len(patch0_data))
                self.assertFalse(np.array_equal(patch0_data, zeros))
                self.assertFalse(np.array_equal(patch1_data, zeros))
            self.assertTrue(np.allclose(patch0_data, patch1_data, atol=1e-6))

    def _checkParticles(self, overlaps, dim, interp):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            if isinstance(overlap, PatchGhostParticleOverlap):
                self._checkPatchGhost(overlap)
                self._checkPatchGhost(overlap.mirror)
            else:
                self._checkLevelGhost(overlap, dim, interp)

    def _checkPatchGhost(self, overlap):
        assert isinstance(overlap, Overlap)

        domain, ghost = overlap.domainPatch, overlap.ghostPatch
        domain_iCell, ghost_iCell = [patch.data("iCell") for patch in [domain, ghost]]
        matching_iCells = []

        if overlap.periodic:
            """
              for periodic overlaps, the iCells do not actually match
              lower_x_patch iCell -1 == upper_x iCell domain_size - 1
            """
            from phare.pp.diagnostic.patch import max_amr_index_for

            max_x_amr_index = max_amr_index_for(domain.patch_level(), "x")
            for nGhost in range(1, overlap.nGhosts + 1):
                shift = nGhost - 1 if overlap.is_mirror else -nGhost
                matching_iCells.extend([shift, max_x_amr_index + shift])
        else:
            """
              otherwise we use the unique list of ghostPatch icells
            """
            uniq_domain_iCell = np.unique(domain_iCell).tolist()
            for i in np.unique(ghost_iCell).tolist():
                for index in np.where(np.array(uniq_domain_iCell) == i)[0]:
                    matching_iCells.append(uniq_domain_iCell[index])

        domain_indices, ghost_indices = [], []
        for i in np.unique(matching_iCells).tolist():
            domain_indices.extend(np.where(np.array(domain_iCell) == i)[0])
            ghost_indices.extend(np.where(np.array(ghost_iCell) == i)[0])

        self.assertEqual(
            ParticleOverlapComparator(domain, domain_indices),
            ParticleOverlapComparator(ghost, ghost_indices),
        )

    def _split_coarse(self, coarseDomain, coarsePatchGhost, dim, interp):
        """Merge patch arrays before passing off to splitter
            there's a chance the overlap covers domain and/or patchghost"""

        arrays = {dataset: [] for dataset in _ParticlePatchData.dataset_keys}
        for patch in [coarseDomain, coarsePatchGhost]:
            for dataset in _ParticlePatchData.dataset_keys:
                arrays[dataset].extend(patch.data(dataset)[:])

        contiguous_t = getattr(  # C++ type = ContiguousParticles<dim, interp>
            tst, "ContiguousParticles_" + str(dim) + "_" + str(interp)
        )
        contiguousParticles = contiguous_t(len(arrays["charge"]))  # charge is always 1d

        for dataset in _ParticlePatchData.dataset_keys:
            object.__setattr__(contiguousParticles, dataset, arrays[dataset][:])

        return contiguousParticles.split()

    def _checkLevelGhost(self, overlap, dim, interp):
        assert isinstance(overlap, Overlap)

        fineLevelGhost = overlap.fineLevelGhostPatch
        coarseSplitParticles = self._split_coarse(
            overlap.coarseDomainPatch, overlap.coarsePatchGhostPatch, dim, interp
        )

        """
        lower_x - overlap.firstCommonFineiCell[0] == finePatchLowerXAMRIndex - nGhosts
        upper_x - overlap.firstCommonFineiCell[0] == finePatchUpperXAMRIndex
        """
        icells = [overlap.firstCommonFineiCell[0] + g for g in range(overlap.nGhosts)]

        patchghost_indices, levelghost_indices = [], []
        for i in icells:
            patchghost_indices.extend(np.where(coarseSplitParticles.iCell == i)[0])
            levelghost_indices.extend(np.where(fineLevelGhost.data("iCell") == i)[0])

        self.assertEqual(len(patchghost_indices), len(levelghost_indices))
        self.assertEqual(
            ParticleOverlapComparator(coarseSplitParticles, patchghost_indices),
            ParticleOverlapComparator(fineLevelGhost, levelghost_indices),
        )


if __name__ == "__main__":
    unittest.main()
