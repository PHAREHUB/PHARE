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

    def test_boundary_conditions_default_none(self):
        # a physical boundary with no explicit dict defaults every location to 'none'
        global_vars.sim = None
        s = simulation.Simulation(**self._mhd_kwargs(boundary_types="periodic"))
        for loc in ("xlower", "xupper"):
            self.assertEqual("none", s.boundary_conditions[loc]["type"])
        self.assertEqual("ideal_gas", s.eos)

    def test_physical_boundary_requires_type(self):
        # a physical boundary left as 'none' must be rejected
        global_vars.sim = None
        with self.assertRaises(KeyError):
            simulation.Simulation(**self._mhd_kwargs(boundary_types="physical"))

    def test_inflow_velocity_scalar_normalized_signed(self):
        # a scalar inflow speed becomes the signed inward-normal component
        global_vars.sim = None
        s = simulation.Simulation(
            **self._mhd_kwargs(
                boundary_types="physical",
                boundary_conditions={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "data": {
                            "velocity": 2.0,
                            "density": 1.0,
                            "pressure": 1.0,
                            "B": [0.5, 1.0, 0.0],
                        },
                    },
                    "xupper": {"type": "super-magnetofast-outflow"},
                },
            )
        )
        vx, vy, vz = s.boundary_conditions["xlower"]["data"]["velocity"]
        self.assertEqual((2.0, 0.0, 0.0), (vx, vy, vz))  # +x inward at lower

    def test_inflow_scalar_B_rejected(self):
        # the magnetic field of an inflow must be a 3-vector, not a scalar
        global_vars.sim = None
        with self.assertRaises((TypeError, ValueError)):
            simulation.Simulation(
                **self._mhd_kwargs(
                    boundary_types="physical",
                    boundary_conditions={
                        "xlower": {
                            "type": "super-magnetofast-inflow",
                            "data": {
                                "velocity": 2.0,
                                "density": 1.0,
                                "pressure": 1.0,
                                "B": 0.5,
                            },
                        },
                        "xupper": {"type": "super-magnetofast-outflow"},
                    },
                )
            )

    def test_callable_B_component_accepted_for_both_inflow_types(self):
        # each component of B (and the other inflow quantities) may be a space-time
        # callable f(x[, y[, z]], t); this is accepted for both super-magnetofast-inflow
        # and free-pressure-inflow (no BC-type-specific restriction).
        global_vars.sim = None
        Bx = lambda x, t: 0.5
        simulation.Simulation(
            **self._mhd_kwargs(
                boundary_types="physical",
                boundary_conditions={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "data": {
                            "velocity": 2.0,
                            "density": 1.0,
                            "pressure": 1.0,
                            "B": [Bx, 1.0, 0.0],
                        },
                    },
                    "xupper": {"type": "super-magnetofast-outflow"},
                },
            )
        )
        global_vars.sim = None
        simulation.Simulation(
            **self._mhd_kwargs(
                boundary_types="physical",
                boundary_conditions={
                    "xlower": {
                        "type": "free-pressure-inflow",
                        "data": {"velocity": 2.0, "density": 1.0, "B": [Bx, 1.0, 0.0]},
                    },
                    "xupper": {"type": "fixed-pressure-outflow",
                               "data": {"pressure": 1.0}},
                },
            )
        )

    def test_inflow_accepts_callable_density(self):
        # density may also be a space-time callable f(x, t), not just a constant
        global_vars.sim = None
        rho = lambda x, t: 1.0
        s = simulation.Simulation(
            **self._mhd_kwargs(
                boundary_types="physical",
                boundary_conditions={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "data": {
                            "velocity": 2.0,
                            "density": rho,
                            "pressure": 1.0,
                            "B": [0.5, 1.0, 0.0],
                        },
                    },
                    "xupper": {"type": "super-magnetofast-outflow"},
                },
            )
        )
        self.assertIn("xlower", s.boundary_conditions)

    def test_inflow_accepts_callable_velocity_component(self):
        # a single velocity component may be a space-time callable f(x, t)
        global_vars.sim = None
        vx = lambda x, t: 2.0
        s = simulation.Simulation(
            **self._mhd_kwargs(
                boundary_types="physical",
                boundary_conditions={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "data": {
                            "velocity": [vx, 0.0, 0.0],
                            "density": 1.0,
                            "pressure": 1.0,
                            "B": [0.5, 1.0, 0.0],
                        },
                    },
                    "xupper": {"type": "super-magnetofast-outflow"},
                },
            )
        )
        self.assertIn("xlower", s.boundary_conditions)

    def test_inflow_accepts_callable_pressure(self):
        # pressure may also be a space-time callable f(x, t)
        global_vars.sim = None
        p = lambda x, t: 1.0
        s = simulation.Simulation(
            **self._mhd_kwargs(
                boundary_types="physical",
                boundary_conditions={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "data": {
                            "velocity": 2.0,
                            "density": 1.0,
                            "pressure": p,
                            "B": [0.5, 1.0, 0.0],
                        },
                    },
                    "xupper": {"type": "super-magnetofast-outflow"},
                },
            )
        )
        self.assertIn("xlower", s.boundary_conditions)

    def test_free_pressure_inflow_accepts_callable_density(self):
        # free-pressure-inflow has no pressure key; density may be callable there too
        global_vars.sim = None
        rho = lambda x, t: 1.0
        s = simulation.Simulation(
            **self._mhd_kwargs(
                boundary_types="physical",
                boundary_conditions={
                    "xlower": {
                        "type": "free-pressure-inflow",
                        "data": {"velocity": 2.0, "density": rho, "B": [0.5, 1.0, 0.0]},
                    },
                    "xupper": {"type": "fixed-pressure-outflow",
                               "data": {"pressure": 1.0}},
                },
            )
        )
        self.assertIn("xlower", s.boundary_conditions)


if __name__ == "__main__":
    unittest.main()
