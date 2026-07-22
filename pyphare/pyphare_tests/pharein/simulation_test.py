import unittest
import numpy as np


import pyphare.pharein.global_vars as global_vars

from pyphare.core import phare_utilities
from pyphare.pharein import simulation
from pyphare.pharein import tagging


class TestSimulation(unittest.TestCase):
    def setUp(self):
        self.cells_array = [80, (80, 40), (80, 40, 12)]
        self.dl_array = [0.1, (0.1, 0.2), (0.1, 0.2, 0.3)]
        self.domain_size_array = [100.0, (100.0, 80.0), (100.0, 80.0, 20.0)]
        self.ndim = [1, 2]  # TODO https://github.com/PHAREHUB/PHARE/issues/232
        self.bcs = [
            "periodic",
            ("periodic", "periodic"),
            ("periodic", "periodic", "periodic"),
        ]
        self.layout = "yee"
        self.time_step = 0.001
        self.time_step_nbr = 1000
        self.final_time = 1.0
        global_vars.sim = None

    def test_dl(self):
        for cells, domain_size, dim, bc in zip(
            self.cells_array, self.domain_size_array, self.ndim, self.bcs
        ):
            j = simulation.Simulation(
                time_step_nbr=self.time_step_nbr,
                boundary_types=bc,
                cells=cells,
                domain_size=domain_size,
                final_time=self.final_time,
            )

            if phare_utilities.none_iterable(domain_size, cells):
                domain_size = phare_utilities.listify(domain_size)
                cells = phare_utilities.listify(cells)

            for d in np.arange(dim):
                self.assertEqual(j.dl[d], domain_size[d] / float(cells[d]))

            global_vars.sim = None

    def test_boundary_conditions(self):
        j = simulation.Simulation(
            time_step_nbr=1000,
            boundary_types="periodic",
            cells=80,
            domain_size=10,
            final_time=1.0,
        )

        for d in np.arange(j.ndim):
            self.assertEqual("periodic", j.boundary_types[d])

    def test_assert_boundary_condition(self):
        simulation.Simulation(
            time_step_nbr=1000,
            boundary_types="periodic",
            cells=80,
            domain_size=10,
            final_time=1000,
        )

    def test_time_step(self):
        s = simulation.Simulation(
            time_step_nbr=1000,
            boundary_types="periodic",
            cells=80,
            domain_size=10,
            final_time=10,
        )
        self.assertEqual(0.01, s.time_step)

    def test_tagging_unknown_top_level_key_raises(self):
        # a typo in a top-level tagging key must fail rather than silently fall back
        with self.assertRaises(ValueError):
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={"method": "lohner", "quantites": {"rho": 0.4}},
            )

    def test_tagging_without_refinement_raises(self):
        # a fully specified tagging= with refinement left at default must not be ignored
        with self.assertRaises(ValueError):
            tagging.resolve_tagging(
                max_nbr_levels=3,
                tagging={"method": "wavelet", "quantities": {"rho": 0.01}},
            )

    def test_tagging_valid_spec_normalizes(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "lohner",
                "quantities": {"B": 0.1},
                "params": {"reltol": 0.02},
            },
        )
        self.assertIsInstance(tag, tagging.LohnerTagging)
        self.assertEqual(tag.method, "lohner")
        self.assertEqual(tag.quantities, [("B", 0.1)])
        self.assertEqual(tag.reltol, 0.02)
        self.assertEqual(tag.abstol, 1e-30)  # unspecified -> field default

    def test_tagging_wavelet_level_scaling_is_bool(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={"method": "wavelet", "params": {"level_scaling": 0}},
        )
        self.assertIsInstance(tag, tagging.WaveletTagging)
        self.assertIs(tag.level_scaling, False)  # coerced 0 -> bool

    def test_tagging_none_when_not_tagging(self):
        self.assertIsNone(tagging.resolve_tagging())


if __name__ == "__main__":
    unittest.main()
