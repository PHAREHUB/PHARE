#!/usr/bin/env python3
"""
With a seed, the initial particles of a cell depend on (seed, population, level, cell) only.
The same run with two different level-0 tilings must give the same particles in every cell.
Run it on several MPI ranks to also cover the rank count and ownership.
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

ndim = 2
interp = 1
cells = (40, 40)
dl = (0.2, 0.2)
ppc = 10
diag_dir = "phare_outputs/seed_tiling"


def config(largest_patch_size, diag_subdir):
    sim = ph.Simulation(
        time_step=0.001,
        time_step_nbr=1,
        cells=cells,
        dl=dl,
        interp_order=interp,
        smallest_patch_size=10,
        largest_patch_size=largest_patch_size,
        diag_options={
            "format": "phareh5",
            "options": {"dir": f"{diag_dir}/{diag_subdir}", "mode": "overwrite"},
        },
        strict=True,
    )

    def density(x, y):
        return 1.0 + 0.5 * np.cos(2 * np.pi * x / (cells[0] * dl[0]))

    def zero(x, y):
        return 0.0 * x

    def bx(x, y):
        return 1.0 + 0.0 * x

    def vth(x, y):
        return 0.3 + 0.1 * np.sin(2 * np.pi * y / (cells[1] * dl[1]))

    ph.MaxwellianFluidModel(
        bx=bx,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 2026},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.ParticleDiagnostics(
        quantity="domain", write_timestamps=[0.0], population_name="protons"
    )
    return sim


def particles_by_cell(diag_subdir):
    """level-0 domain particles at t=0 as {cell: sorted rows}, and the patch boxes"""
    hier = hierarchy_from(
        h5_filename=f"{diag_dir}/{diag_subdir}/ions_pop_protons_domain.h5",
        times=f"{0:.10f}",
    )
    boxes, rows = [], {}
    for patch in hier.level(0).patches:
        boxes.append((tuple(patch.box.lower), tuple(patch.box.upper)))
        particles = patch.patch_datas["protons_domain"].dataset
        data = np.column_stack(
            (
                particles.iCells,
                particles.deltas,
                particles.v,
                particles.weights,
            )
        )
        for row in data:
            rows.setdefault(tuple(row[:ndim].astype(int)), []).append(tuple(row))
    return sorted(boxes), {cell: sorted(r) for cell, r in rows.items()}


class SeedTilingTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_with(self, largest_patch_size, diag_subdir):
        ph.global_vars.sim = None
        Simulator(config(largest_patch_size, diag_subdir)).initialize().reset()

    def test_particles_do_not_depend_on_the_tiling(self):
        self.run_with(40, "large_patches")
        self.run_with(13, "small_patches")

        if cpp.mpi_rank() != 0:
            return

        boxes_a, cells_a = particles_by_cell("large_patches")
        boxes_b, cells_b = particles_by_cell("small_patches")

        self.assertNotEqual(boxes_a, boxes_b, "the two runs have the same tiling")
        self.assertEqual(len(cells_a), cells[0] * cells[1])
        self.assertEqual(cells_a.keys(), cells_b.keys())
        for cell, rows in cells_a.items():
            self.assertEqual(len(rows), ppc)
            self.assertEqual(rows, cells_b[cell], f"cell {cell}")


if __name__ == "__main__":
    unittest.main()
