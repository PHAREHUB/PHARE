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

    def test_time_step(self):
        s = simulation.Simulation(
            time_step_nbr=1000,
            boundary_types="periodic",
            cells=80,
            domain_size=10,
            final_time=10,
        )
        self.assertEqual(0.01, s.time_step)

    def test_inner_boundary_sphere(self):
        s = simulation.Simulation(
            time_step_nbr=100,
            boundary_types=("periodic", "periodic"),
            cells=(80, 40),
            dl=(0.1, 0.2),
            inner_boundary={"name": "sphere_boundary", "shape": "sphere", "center": (1.0, 2.0), "radius": 0.5},
            final_time=1.0,
        )
        self.assertEqual("sphere", s.inner_boundary["shape"])
        self.assertEqual("sphere_boundary", s.inner_boundary["name"])
        self.assertEqual([1.0, 2.0], s.inner_boundary["center"])
        self.assertEqual(0.5, s.inner_boundary["radius"])
        global_vars.sim = None

    def test_inner_boundary_plane(self):
        s = simulation.Simulation(
            time_step_nbr=100,
            boundary_types=("periodic", "periodic", "periodic"),
            cells=(80, 40, 12),
            dl=(0.1, 0.2, 0.3),
            inner_boundary={
                "name": "plane_boundary",
                "shape": "plane",
                "point": (1.0, 2.0, 3.0),
                "normal": (0.0, 0.0, 2.0),
            },
            final_time=1.0,
        )
        self.assertEqual("plane", s.inner_boundary["shape"])
        self.assertEqual("plane_boundary", s.inner_boundary["name"])
        self.assertEqual([1.0, 2.0, 3.0], s.inner_boundary["point"])
        self.assertEqual([0.0, 0.0, 2.0], s.inner_boundary["normal"])
        global_vars.sim = None

    def test_inner_boundary_invalid(self):
        with self.assertRaises(ValueError):
            simulation.Simulation(
                time_step_nbr=100,
                boundary_types="periodic",
                cells=80,
                dl=0.1,
                inner_boundary={"shape": "sphere", "center": (0.0,), "radius": -1.0},
                final_time=1.0,
            )
        global_vars.sim = None


if __name__ == "__main__":
    unittest.main()
