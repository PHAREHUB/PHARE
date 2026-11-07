#!/usr/bin/env python3
"""
With a seed, the initial particles of a cell do not depend on the patches, so a run must not
depend on the patch layout: the same deck run with small patches and with one patch
per level must agree on every level.

Here: a uniform plasma, two levels (one refinement box), compared between 8-cell patches and one
patch per level, in 1D and 2D at every interpolation order. The number of domain particles per
level must be identical at every step. Level ghost particles are held by every patch whose ghost
box covers them; one entering the domain of a same-level patch belongs to that patch, the others
must drop their copy or it is duplicated.
"""

import unittest

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core.box import Box1D, Box2D
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = {1: 64, 2: 32}
refinement_box = {1: (16, 47), 2: (8, 23)}  # same bounds along every direction
ppc = {1: 128, 2: 64}
dl = 0.125
time_step = 0.00125
time_step_nbr = 5
small_patch_size = 8

cases = [(ndim, interp) for ndim in (1, 2) for interp in (1, 2, 3)]


def config(test, ndim, interp, patch_size):
    lower, upper = refinement_box[ndim]
    box = Box1D if ndim == 1 else Box2D
    sim = test.simulation(
        cells=[cells[ndim]] * ndim,
        dl=[dl] * ndim,
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        interp_order=interp,
        smallest_patch_size=patch_size,
        largest_patch_size=patch_size,
        refinement_boxes={"L0": {"B0": box(lower, upper)}},
        resistivity=0.0,
        hyper_resistivity=0.001,
        diag_options={
            "format": "phareh5",
            "options": {"dir": f"patch_layout/p{patch_size}", "mode": "overwrite"},
        },
    )

    def const(value):
        return lambda *xyz: value + np.zeros_like(xyz[0])

    ph.MaxwellianFluidModel(
        bx=const(0.0),
        by=const(0.0),
        bz=const(1.0),
        protons={
            "charge": 1,
            "mass": 1,
            "density": const(1.0),
            "vbulkx": const(0.0),
            "vbulky": const(0.0),
            "vbulkz": const(0.0),
            "vthx": const(0.3),
            "vthy": const(0.3),
            "vthz": const(0.3),
            "nbr_part_per_cell": ppc[ndim],
            "init": {"seed": 11},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.ParticleDiagnostics(
        quantity="domain",
        write_timestamps=timestamps(),
        population_name="protons",
    )
    return sim


def timestamps():
    return np.arange(time_step_nbr + 1) * time_step


def particle_counts(diag_dir):
    """number of domain particles per (time, level)"""
    counts = {}
    for time in timestamps():
        hier = hierarchy_from(
            h5_filename=f"{diag_dir}/ions_pop_protons_domain.h5", times=f"{time:.10f}"
        )
        for ilvl, level in hier.levels(time).items():
            counts[(time, ilvl)] = sum(
                patch.patch_datas["protons_domain"].dataset.size()
                for patch in level.patches
            )
    return counts


@ddt
class PatchLayoutIndependenceTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_layout(self, ndim, interp, patch_size):
        sim = config(self, ndim, interp, patch_size)
        Simulator(sim).run().reset()
        nbr_patches = len(
            hierarchy_from(
                h5_filename=f"{sim.diag_options['options']['dir']}/ions_pop_protons_domain.h5",
                times=f"{0:.10f}",
            ).level(0).patches
        )
        return sim.diag_options["options"]["dir"], nbr_patches

    @data(*cases)
    @unpack
    def test_particle_count_per_level(self, ndim, interp):
        print(f"{self._testMethodName}: {ndim}D interp {interp}")
        small, nbr_small = self.run_layout(ndim, interp, small_patch_size)
        whole, nbr_whole = self.run_layout(ndim, interp, cells[ndim])

        if cpp.mpi_rank() != 0:
            return

        self.assertEqual(nbr_whole, 1)
        self.assertGreater(nbr_small, 1)

        small_counts, whole_counts = particle_counts(small), particle_counts(whole)
        self.assertEqual(small_counts.keys(), whole_counts.keys())
        self.assertIn((timestamps()[-1], 1), small_counts)
        for key, count in whole_counts.items():
            self.assertEqual(small_counts[key], count, f"(time, level) {key}")


if __name__ == "__main__":
    unittest.main()
