import unittest
import numpy as np


import pyphare.pharein.global_vars as global_vars

from pyphare.core import phare_utilities
from pyphare.pharein import simulation


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

    def make(self, **kwargs):
        global_vars.sim = None
        return simulation.Simulation(
            time_step_nbr=10, time_step=0.001, cells=80, dl=0.1, **kwargs
        )

    def test_collisions_defaults(self):
        self.assertIsNone(self.make().collisions)
        collisions = self.make(collisions={"nu0": 1e-3}).collisions
        self.assertEqual(
            collisions,
            {
                "nu0": 1e-3,
                "every": 1,
                "scheme": "nanbu",
                "exact_conservation": True,
                "seed": 0,
            },
        )

    def test_collisions_invalid(self):
        for bad in [
            {},
            {"nu0": -1.0},
            {"nu0": 1.0, "every": 0},
            {"nu0": 1.0, "scheme": "boltzmann"},
            {"nu0": 1.0, "seed": -1},
            {"nu0": 1.0, "frequency": 2},
        ]:
            with self.assertRaises(ValueError, msg=str(bad)):
                self.make(collisions=bad)

    def test_collision_nu0_for_the_solar_wind(self):
        # n = 5 /cm3, B = 5 nT: (c / VA)^4 / (4 pi n0 di^3) = 2.15e-8
        self.assertAlmostEqual(
            simulation.collision_nu0(5.0, 5.0, 1.0) / 2.1512e-8, 1.0, places=3
        )

    def test_time_step(self):
        s = simulation.Simulation(
            time_step_nbr=1000,
            boundary_types="periodic",
            cells=80,
            domain_size=10,
            final_time=10,
        )
        self.assertEqual(0.01, s.time_stepper.time_step)


if __name__ == "__main__":
    unittest.main()
