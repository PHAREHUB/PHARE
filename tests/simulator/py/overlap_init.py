#!/usr/bin/env python3
#
# formatted with black

from tests.simulator import test_simulator as tst
from tests.simulator.py import InitValueValidation

from ddt import ddt, data

import unittest, os, phare.pharein as ph, numpy as np

from tests.diagnostic import dump_all_diags

from phare.pp.diagnostic.data.overlap import Overlap
from phare.pp.diagnostic.data.em_overlap import EMOverlap
from phare.pp.diagnostic.data.fluid_overlap import FluidOverlap
from phare.pp.diagnostic.data.particle_overlap import ParticleOverlap, ParticleOverlapComparator
from phare.pp.diagnostic.data.particle_level_overlap import LevelParticleOverlap
from phare.pp.diagnostic.data.particle_patch_overlap import DomainParticleOverlap

from phare.pp.diagnostics import Diagnostics, Diagnostic, Patch, _EM, _Fluid, _Particle


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
        add_to_dict({"id": 2, "refinement_boxes": {"L0": {"B0": [(10,), (14,)]}},}),
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
            diags = self._simulate_diagnostics(dim, interp, input)
            self._checkEM(EMOverlap.get(diags))
            self._checkFluid(FluidOverlap.get(diags))
            self._checkParticles(ParticleOverlap.get(diags), dim, interp)

    def _checkEM(self, overlaps):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            key, offsets = overlap.key, overlap.offsets
            for i in range(overlap.sizes[0]):
                x0, x1 = (
                    overlap.p0.dtype.get()[key][offsets[0][0] + i],
                    overlap.p1.dtype.get()[key][offsets[0][1] + i],
                )
                self.assertTrue(Overlap._float_equal(x0, x1))

    def _checkFluid(self, overlaps):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            key, offsets, nGhosts = overlap.key, overlap.offsets, overlap.nGhosts
            if overlap.sizes[0] == nGhosts * 2:
                x0, x1 = (
                    overlap.p0.dtype.get()[key][offsets[0][0] + nGhosts],
                    overlap.p1.dtype.get()[key][offsets[0][1] + nGhosts],
                )
                self.assertTrue(Overlap._float_equal(x0, x1))

    def _checkParticles(self, overlaps, dim, interp):
        self.assertTrue(len(overlaps))
        for overlap in overlaps:
            self._checkPatchGhost(overlap) if isinstance(
                overlap, DomainParticleOverlap
            ) else self._checkLevelGhost(overlap, dim, interp)

    def _checkPatchGhost(self, overlap):
        assert isinstance(overlap, Overlap)

        domain, ghost, coords = overlap.domain, overlap.ghost, overlap.origin
        domain_iCell, ghost_iCell = [i.dtype.get()["iCell"] for i in [domain, ghost]]

        uniq_gic = np.unique(ghost_iCell).tolist()
        max_idx = domain.patch_level.cells
        max_X_idx = max_idx[0]

        # upper periodicity
        for i in uniq_gic.copy():  # don't edit list being checked
            if i >= max_X_idx:
                uniq_gic.append(max_X_idx - i)

        icells = [  # filter icell values in ghost dataset not in overlap
            i
            for i in uniq_gic
            if coords[0] - overlap.nGhosts <= i <= coords[0] + overlap.nGhosts
        ]

        # lower periodicity
        for i in icells.copy():  # don't edit list being checked
            if i <= 0:
                icells.append(max_X_idx + i)

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

        arrays = {dataset: [] for dataset in _Particle.datasets}
        for patch in [coarseDomain, coarsePatchGhost]:
            [
                arrays[dataset].extend(patch.dtype.get()[dataset][:])
                for dataset in _Particle.datasets
            ]
        contiguous_t = getattr(  # C++ type = ContiguousParticles<dim, interp>
            tst, "ContiguousParticles_" + str(dim) + "_" + str(interp)
        )
        contiguousParticles = contiguous_t(len(arrays["charge"]))  # charge is always 1d
        [
            object.__setattr__(contiguousParticles, dataset, arrays[dataset][:])
            for dataset in _Particle.datasets
        ]
        return contiguousParticles.split()

    def _checkLevelGhost(self, overlap, dim, interp):
        assert isinstance(overlap, Overlap)

        if "iCell" not in overlap.fineLevelGhost.dtype.keys():
            """level ghost non existent - probably patchghost overlap on fine level
                can happen with three+ contiguous fine patches, the inner will lack any
                level ghost values
            """
            return

        fineLevelGhost = overlap.fineLevelGhost
        coarseSplitParticles = self._split_coarse(
            overlap.coarseDomain, overlap.coarsePatchGhost, dim, interp
        )

        icells = [overlap.offsets[0] + g for g in range(overlap.sizes[0])]

        def where(cells):
            lst = []
            for i in icells:
                lst.extend(np.where(np.array(cells) == i)[0])
            return lst

        gcells = where(fineLevelGhost.dtype.get()["iCell"])
        if len(gcells) == 0:
            return  # level ghost non existent - probably patchghost overlap on fine level

        pcells = where(coarseSplitParticles.iCell)
        self.assertEqual(len(pcells), len(gcells))
        self.assertEqual(
            ParticleOverlapComparator(coarseSplitParticles, pcells),
            ParticleOverlapComparator(fineLevelGhost, gcells),
        )


if __name__ == "__main__":
    unittest.main()
