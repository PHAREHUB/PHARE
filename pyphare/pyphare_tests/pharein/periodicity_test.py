"""
The periodicity check of the user functions (MaxwellianFluidModel.validate), shared by the
hybrid and fully kinetic models: a function that is not periodic along a periodic direction
makes the two copies of a border node differ from t = 0 (each patch evaluates it on its side).
"""

import unittest
from unittest.mock import patch

import numpy as np

from pyphare.pharein import global_vars, simulation
from pyphare.pharein.maxwellian_fluid_model import MaxwellianFluidModel

cells = (40, 20)
dl = (0.2, 0.2)
Lx, Ly = cells[0] * dl[0], cells[1] * dl[1]


def constant(value):
    return lambda *xyz: value + 0.0 * xyz[0]


def bump_at_y_border(x, y):
    """1 plus a bump on (Lx/2, 0), periodic in x, with no periodic image in y: B jumps across
    y = 0 near x = Lx/2 only (by 0.1), and by 4e-19 on the line x = 0"""
    periodic_in_x = np.exp(20.0 * (np.cos(2 * np.pi * (x - 0.5 * Lx) / Lx) - 1.0))
    return 1.0 + 0.1 * periodic_in_x * np.exp(-(y**2))


def periodic_wave(x, y):
    return 1.0 + 0.1 * np.cos(2 * np.pi * x / Lx) * np.sin(2 * np.pi * y / Ly)


def model(bx):
    population = {
        "charge": 1,
        "density": constant(1.0),
        **{v: constant(0.0) for v in ("vbulkx", "vbulky", "vbulkz")},
        **{v: constant(0.1) for v in ("vthx", "vthy", "vthz")},
    }
    return MaxwellianFluidModel(
        bx=bx, by=constant(0.0), bz=constant(0.0), protons=population
    )


class TestPeriodicityCheck(unittest.TestCase):
    def setUp(self):
        global_vars.sim = None

    def tearDown(self):
        global_vars.sim = None

    def sim(self, strict=False):
        return simulation.Simulation(
            time_step_nbr=10, time_step=0.001, cells=cells, dl=dl, strict=strict
        )

    def printed(self, bx):
        with patch("pyphare.cpp.print_rank0") as print_rank0:
            model(bx)
        printed = (str(a) for c in print_rank0.call_args_list for a in c.args)
        return " ".join(line for line in printed if "Warning" in line)

    def test_a_periodic_function_is_not_flagged(self):
        self.sim()
        self.assertEqual("", self.printed(periodic_wave))

    def test_a_jump_away_from_the_axes_is_flagged(self):
        self.sim()
        printed = self.printed(bump_at_y_border)
        self.assertIn("Bx along Y", printed)
        self.assertNotIn("Bx along X", printed)

    def test_strict_refuses_it(self):
        self.sim(strict=True)
        with (
            self.assertRaisesRegex(RuntimeError, "not periodic"),
            patch("pyphare.cpp.print_rank0"),
        ):
            model(bump_at_y_border)


if __name__ == "__main__":
    unittest.main()
