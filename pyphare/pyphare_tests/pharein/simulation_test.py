import unittest
import numpy as np


import pyphare.pharein.global_vars as global_vars

from pyphare.core import phare_utilities
from pyphare.pharein import simulation
from pyphare.pharein.mhd_model import MHDModel


class TestSimulation(unittest.TestCase):
    def setUp(self):
        self.cells_array = [80, (80, 40), (80, 40, 12)]
        self.dl_array = [0.1, (0.1, 0.2), (0.1, 0.2, 0.3)]
        self.domain_size_array = [100.0, (100.0, 80.0), (100.0, 80.0, 20.0)]
        self.ndim = [1, 2]  # TODO https://github.com/PHAREHUB/PHARE/issues/232
        self.layout = "yee"
        self.time_step = 0.001
        self.time_step_nbr = 1000
        self.final_time = 1.0
        global_vars.sim = None

    def test_dl(self):
        for cells, domain_size, dim in zip(
            self.cells_array, self.domain_size_array, self.ndim
        ):
            j = simulation.Simulation(
                time_step_nbr=self.time_step_nbr,
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

    def test_boundaries_default_periodic(self):
        j = simulation.Simulation(
            time_step_nbr=1000,
            cells=80,
            domain_size=10,
            final_time=1.0,
        )

        for d in np.arange(j.ndim):
            self.assertTrue(j.periodicities[d])

    def test_boundary_types_kwarg_rejected(self):
        with self.assertRaises(ValueError):
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
            cells=80,
            domain_size=10,
            final_time=10,
        )
        self.assertEqual(0.01, s.time_stepper.time_step)

    # ---- physical outer boundary conditions ----------------------------------

    def _mhd_kwargs(self, **overrides):
        kwargs = dict(
            time_step_nbr=1,
            cells=80,
            domain_size=10,
            final_time=1.0,
            model_options=["MHDModel"],
            mhd_timestepper="SSPRK4_5",
            reconstruction="WENOZ",
            limiter="None",
            riemann="Rusanov",
        )
        kwargs.update(overrides)
        return kwargs

    def test_boundaries_default_none(self):
        # without a boundaries dict every location is periodic and defaults to 'none'
        global_vars.sim = None
        s = simulation.Simulation(**self._mhd_kwargs())
        self.assertEqual([True], s.periodicities)
        for loc in ("xlower", "xupper"):
            self.assertEqual("none", s.boundaries[loc].type)

    def test_physical_direction_requires_both_locations(self):
        # giving only one side of a direction must be rejected
        global_vars.sim = None
        with self.assertRaises(ValueError):
            simulation.Simulation(
                **self._mhd_kwargs(boundaries={"xlower": {"type": "open"}})
            )

    def test_inflow_velocity_scalar_normalized_signed(self):
        # a scalar inflow speed becomes the signed inward-normal component
        global_vars.sim = None
        s = simulation.Simulation(
            **self._mhd_kwargs(
                boundaries={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "velocity": 2.0,
                        "density": 1.0,
                        "pressure": 1.0,
                        "B": [0.5, 1.0, 0.0],
                    },
                    "xupper": {"type": "open"},
                },
            )
        )
        self.assertEqual([False], s.periodicities)
        vx, vy, vz = s.boundaries["xlower"].velocity
        self.assertEqual((2.0, 0.0, 0.0), (vx, vy, vz))  # +x inward at lower

    def test_inflow_scalar_B_rejected(self):
        # the magnetic field of an inflow must be a 3-vector, not a scalar
        global_vars.sim = None
        with self.assertRaises((TypeError, ValueError)):
            simulation.Simulation(
                **self._mhd_kwargs(
                    boundaries={
                        "xlower": {
                            "type": "super-magnetofast-inflow",
                            "velocity": 2.0,
                            "density": 1.0,
                            "pressure": 1.0,
                            "B": 0.5,
                        },
                        "xupper": {"type": "open"},
                    },
                )
            )

    def _inflow_sim(self, velocity=2.0, B=(0.5, 1.0, 0.0), **overrides):
        global_vars.sim = None
        return simulation.Simulation(
            **self._mhd_kwargs(
                boundaries={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "velocity": velocity,
                        "density": 1.0,
                        "pressure": 1.0,
                        "B": list(B),
                    },
                    "xupper": {"type": "open"},
                },
                **overrides,
            )
        )

    def test_inflow_consistent_normal_b_accepted(self):
        s = self._inflow_sim()
        MHDModel(bx=lambda x: 0.5 + 0.1 * np.sin(np.pi * x / 10.0), by=lambda x: 1.0)
        self.assertIs(s, global_vars.sim)
        self.assertIsNotNone(s.model)

    def test_inflow_inconsistent_normal_b_rejected(self):
        self._inflow_sim()
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x: 0.7 + x * 0, by=lambda x: 1.0)

    def test_inflow_oblique_inconsistent_normal_b_rejected(self):
        self._inflow_sim(velocity=[2.0, 0.5, 0.0])
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x: 0.7 + x * 0, by=lambda x: 1.0)

    def test_inflow_normal_b_checked_along_whole_face(self):
        self._inflow_sim(cells=(80, 40), domain_size=(10.0, 8.0))
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x, y: 0.5 + 0.01 * np.sin(2 * np.pi * y / 8.0))

    def test_inflow_normal_b_checked_at_upper_face(self):
        global_vars.sim = None
        simulation.Simulation(
            **self._mhd_kwargs(
                boundaries={
                    "xlower": {"type": "open"},
                    "xupper": {
                        "type": "super-magnetofast-inflow",
                        "velocity": 2.0,
                        "density": 1.0,
                        "pressure": 1.0,
                        "B": [0.5, 1.0, 0.0],
                    },
                },
            )
        )
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x: np.where(x > 9.0, 0.6, 0.5), by=lambda x: 1.0)


if __name__ == "__main__":
    unittest.main()
