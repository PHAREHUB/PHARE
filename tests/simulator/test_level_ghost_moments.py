#!/usr/bin/env python3
"""
Ohm's law reads the ion charge density and E-node flux one node beyond the patch. On a refined
level the level ghost particles cannot complete the nodes beyond the level, so their ghosts are
refined from the coarser level, in space and time (#1577).

Check: the fine ghost node that coincides with a coarse node outside the level must hold exactly
the coarse density there. At the end of a coarse step the last fill uses the coarse data at the
new time, and refining onto a node that coincides with a coarse node is a copy. Coarse nodes
outside the finer level are not changed by the synchronisation, so the equality is exact.
Without the refinement this node holds a partial deposit of level ghost particles.

The first ghost node, the one Ohm reads, sits between two coarse nodes, one of which is on the
level border and is overwritten by the synchronisation after the fill: it is checked loosely
against the average of its two coarse neighbours.

The same exact check holds for the E-node charge flux components that are primal along x (y and z
in 1D), with a drift along y so that the flux is not only noise.

The level border nodes of the E-node flux are completed by the level ghost particles, like the
density: there the flux over the density must be the drift velocity (about half of it if the level
ghost deposit of the E-node flux were missing), and the density must be the initial profile (about
half of it if the level ghost particles were not deposited at all).

The population flux the outputs compute (HybridLazyMoments, not the time loop) deposits the level
ghost particles with its own ownership map: on the level border it too must give the drift velocity
over the density the time loop deposited.
"""

import unittest

import numpy as np
from ddt import data, ddt

import pyphare.pharein as ph
from pyphare.core.box import Box1D
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

refined = (10, 25)  # level 0 cells -> level 1 primal nodes 20..52
ratio = 2
dl = 0.4


def density(x):
    return 1.0 + 0.3 * np.sin(2 * np.pi * x / 16.0)


def config(interp, vy=0.0):
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=0.005,
        time_step_nbr=2,
        cells=40,
        dl=dl,
        interp_order=interp,
        refinement_boxes={"L0": {"B0": Box1D(*refined)}},
        hyper_resistivity=0.001,
        strict=True,
    )

    def constant(value):
        return lambda x: value + 0.0 * x

    ph.MaxwellianFluidModel(
        bx=constant(1.0),
        by=constant(0.0),
        bz=constant(0.0),
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": constant(0.5),
            "vbulky": constant(vy),
            "vbulkz": constant(0.0),
            "vthx": constant(0.1),
            "vthy": constant(0.1),
            "vthz": constant(0.1),
            "nbr_part_per_cell": 50,
            "init": {"seed": 1},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.1)
    return sim


def node_values(patches):
    """AMR primal node index -> list of the values every patch holds there, ghosts included"""
    values = {}
    for patch in patches:
        g = patch.nGhosts
        data = np.asarray(patch.data)
        first = int(patch.lower[0]) - g
        for k, v in enumerate(data):
            values.setdefault(first + k, []).append(float(v))
    return values


@ddt
class LevelGhostMomentsTest(SimulatorTest):
    def tearDown(self):
        if hasattr(self, "simulator"):
            self.simulator.reset()
        super().tearDown()

    @data(1, 2, 3)
    def test_density_beyond_the_level_is_the_coarse_density(self, interp):
        self.simulator = Simulator(config(interp)).initialize()
        for _ in range(2):
            self.simulator.advance()

        wrangler = self.simulator.data_wrangler()
        coarse = node_values(wrangler.getPatchLevel(0).getDensity())
        fine = node_values(wrangler.getPatchLevel(1).getDensity())

        # second ghost node beyond each edge of level 1, on coarse nodes outside level 1
        lower_edge, upper_edge = refined[0] * ratio, (refined[1] + 1) * ratio
        for fine_node in (lower_edge - 2, upper_edge + 2):
            coarse_node = fine_node // ratio
            self.assertIn(fine_node, fine)
            expected = coarse[coarse_node][0]
            for value in fine[fine_node]:
                self.assertEqual(
                    value, expected, f"fine node {fine_node} interp {interp}"
                )

        # first ghost node beyond each edge, the one Ohm reads, between two coarse nodes
        # (with the fill: within 0.02 here; without it: off by 0.06 to 0.65)
        for fine_node in (lower_edge - 1, upper_edge + 1):
            left = fine_node // ratio
            expected = 0.5 * (coarse[left][0] + coarse[left + 1][0])
            for value in fine[fine_node]:
                self.assertLess(
                    abs(value - expected),
                    0.04,
                    f"fine node {fine_node} interp {interp}",
                )

    @data(1, 2, 3)
    def test_E_node_flux_beyond_the_level_is_the_coarse_flux(self, interp):
        self.simulator = Simulator(config(interp, vy=0.3)).initialize()
        for _ in range(2):
            self.simulator.advance()

        wrangler = self.simulator.data_wrangler()
        coarse = node_values(wrangler.getPatchLevel(0).getChargeFluxE("y"))
        fine = node_values(wrangler.getPatchLevel(1).getChargeFluxE("y"))

        lower_edge, upper_edge = refined[0] * ratio, (refined[1] + 1) * ratio
        for fine_node in (lower_edge - 2, upper_edge + 2):
            coarse_node = fine_node // ratio
            self.assertIn(fine_node, fine)
            expected = coarse[coarse_node][0]
            self.assertGreater(abs(expected), 0.1)  # the drift, not noise
            for value in fine[fine_node]:
                self.assertEqual(
                    value, expected, f"fine node {fine_node} interp {interp}"
                )

    @data(1, 2, 3)
    def test_E_node_flux_on_the_level_border_is_complete(self, interp):
        vy = 0.3
        self.simulator = Simulator(config(interp, vy=vy)).initialize()
        for _ in range(2):
            self.simulator.advance()

        level = self.simulator.data_wrangler().getPatchLevel(1)
        fluxE = node_values(level.getChargeFluxE("y"))
        n_level = node_values(level.getDensity())

        # measured: 0.28 to 0.30
        lower_edge, upper_edge = refined[0] * ratio, (refined[1] + 1) * ratio
        for node in (lower_edge, lower_edge + 1, upper_edge - 1, upper_edge):
            for f, n in zip(fluxE[node], n_level[node]):
                self.assertLess(abs(f / n - vy), 0.05, f"node {node} interp {interp}")

        # the ratio above cannot see a level ghost deposit missing from both: the density on the
        # border node alone (at order 1 the next node gets no level ghost contribution) against the
        # initial profile, 50 particles per cell (without the deposit: about 0.5)
        for node in (lower_edge, upper_edge):
            n0 = density(node * dl / ratio)
            for n in n_level[node]:
                self.assertLess(abs(n / n0 - 1), 0.25, f"node {node} interp {interp}")

    @data(1, 2, 3)
    def test_output_population_flux_on_the_level_border_is_complete(self, interp):
        vy = 0.3
        self.simulator = Simulator(config(interp, vy=vy)).initialize()
        for _ in range(2):
            self.simulator.advance()

        level = self.simulator.data_wrangler().getPatchLevel(1)
        flux = node_values(level.getFy("protons"))
        n_level = node_values(level.getPopDensities()["protons"])

        # without the level ghost deposit of the outputs' flux: about vy / 2
        lower_edge, upper_edge = refined[0] * ratio, (refined[1] + 1) * ratio
        for node in (lower_edge, upper_edge):
            for f, n in zip(flux[node], n_level[node]):
                self.assertLess(abs(f / n - vy), 0.05, f"node {node} interp {interp}")


if __name__ == "__main__":
    unittest.main()
