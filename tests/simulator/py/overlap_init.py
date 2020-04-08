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
from phare.pp.diagnostic.data.particle_level_overlap import LevelParticleOverlap
from phare.pp.diagnostic.data.particle_patch_overlap import DomainParticleOverlap

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
    max_interp = 3

    def add_to_dict(dic):
        dic.update(InitValueValidation.diag_options(diag_out_dir))
        dic.update({"diags_fn": lambda model: dump_all_diags(model.populations)})
        return dic

    valid1D = [
        add_to_dict(
            {
                "id": 0,
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
                "id": 1,
                "refinement_boxes": {
                    "L0": {"B0": [(5,), (55,)]},
                    "L1": {"B0": [(20,), (100,)]},
                },
                "max_nbr_levels": 3,
            }
        ),
        add_to_dict({"id": 2, "refinement_boxes": {"L0": {"B0": [(5,), (54,)]}},}),
    ]

    @data(*valid1D)
    def test_1d(self, input):
        dim = 1
        min_interp = OverlapValueValidation.min_interp
        max_interp = OverlapValueValidation.max_interp + 1
        for interp in range(min_interp, max_interp):
            # blocks conflicts between runs - possible bug
            input["diag_options"]["options"]["dir"] = (
                diag_out_dir
                + "_"
                + str(dim)
                + "_"
                + str(interp)
                + "_"
                + str(input["id"])
            )
            diags = self.runAndDump(dim, interp, input)
            self._checkEM(getEMOverlapsFrom(diags))
            self._checkFluid(getFluidOverlapsFrom(diags))
            self._checkParticles(getParticleOverlapsFrom(diags), dim, interp)

    def _checkEM(self, overlaps):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            patch0_data, patch1_data = overlap.shared_data()
            self.assertTrue(np.allclose(patch0_data, patch1_data))

    def _checkFluid(self, overlaps):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            patch0_data, patch1_data = overlap.shared_data()
            self.assertTrue(np.allclose(patch0_data, patch1_data, atol=1e-07))

    def _checkParticles(self, overlaps, dim, interp):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            if isinstance(overlap, DomainParticleOverlap):
                self._checkPatchGhost(overlap)
                self._checkPatchGhost(overlap.mirror)
            else:
                self._checkLevelGhost(overlap, dim, interp)

    def _checkPatchGhost(self, overlap):
        assert isinstance(overlap, Overlap)
        from phare.pp.diagnostic.patch import max_amr_index_for

        domain, ghost = overlap.domainPatch, overlap.ghostPatch

        domain_iCell, ghost_iCell = [
            i.patch_data.data("iCell") for i in [domain, ghost]
        ]

        uniq_ghostiCell = np.unique(ghost_iCell).tolist()
        uniq_domain_iCell = np.unique(domain_iCell).tolist()

        max_x_amr_idx = max_amr_index_for(domain.patch_level, "x")

        icells = []

        if overlap.periodic:
            for index in [0, max_x_amr_idx]:
                icells.append(index)
                for nGhost in range(1, overlap.nGhosts + 1):
                    icells.append(index - nGhost)
        else:
            for i in uniq_ghostiCell:
                for index in np.where(np.array(uniq_domain_iCell) == i)[0]:
                    icells.append(uniq_domain_iCell[index])

        icells = np.unique(icells).tolist()

        pcells, gcells = [], []
        for i in icells:
            pcells.extend(np.where(np.array(domain_iCell) == i)[0])
            gcells.extend(np.where(np.array(ghost_iCell) == i)[0])

        self.assertTrue(len(pcells) > 0 and len(gcells) > 0)
        self.assertEqual(len(pcells), len(gcells))

        self.assertEqual(
            ParticleOverlapComparator(domain, pcells),
            ParticleOverlapComparator(ghost, gcells),
        )

    def _split_coarse(self, coarseDomain, coarsePatchGhost, dim, interp):
        """Merge patch arrays before passing off to splitter
            there's a chance the overlap covers domain and/or patchghost"""

        arrays = {dataset: [] for dataset in _ParticlePatchData.dataset_keys}
        for patch in [coarseDomain, coarsePatchGhost]:
            [
                arrays[dataset].extend(patch.patch_data.data(dataset)[:])
                for dataset in _ParticlePatchData.dataset_keys
            ]
        contiguous_t = getattr(  # C++ type = ContiguousParticles<dim, interp>
            tst, "ContiguousParticles_" + str(dim) + "_" + str(interp)
        )
        contiguousParticles = contiguous_t(len(arrays["charge"]))  # charge is always 1d
        [
            object.__setattr__(contiguousParticles, dataset, arrays[dataset][:])
            for dataset in _ParticlePatchData.dataset_keys
        ]
        return contiguousParticles.split()

    def _checkLevelGhost(self, overlap, dim, interp):
        assert isinstance(overlap, Overlap)

        fineLevelGhost = overlap.fineLevelGhostPatch
        coarseSplitParticles = self._split_coarse(
            overlap.coarseDomainPatch, overlap.coarsePatchGhostPatch, dim, interp
        )

        icells = [overlap.sizes[0] + g for g in range(overlap.nGhosts)]

        def where(cells):
            lst = []
            for i in icells:
                lst.extend(np.where(np.array(cells) == i)[0])
            return lst

        gcells = where(fineLevelGhost.patch_data.data("iCell"))
        pcells = where(coarseSplitParticles.iCell)

        self.assertEqual(len(pcells), len(gcells))
        self.assertEqual(
            ParticleOverlapComparator(coarseSplitParticles, pcells),
            ParticleOverlapComparator(fineLevelGhost, gcells),
        )


if __name__ == "__main__":
    unittest.main()
