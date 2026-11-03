#!/usr/bin/env python3
#
# the electron closure parameters given in python must reach the C++ dictionary
# whatever their numeric type (int, float, numpy scalars), and anything else
# must be refused with a message naming the entry

import unittest

import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import startMPI

ph.NO_GUI()


def setup_simulation(density=lambda x: 1.0, sim_args=None, **electrons):
    ph.global_vars.sim = None
    sim = ph.Simulation(
        **(sim_args or {}),
        time_step_nbr=1,
        final_time=0.001,
        cells=20,
        dl=0.3,
    )

    def one(x):
        return 1.0

    def zero(x):
        return 0.0

    ph.MaxwellianFluidModel(
        bx=one,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": one,
            "vthy": one,
            "vthz": one,
        },
    )
    ph.ElectronModel(**electrons)
    return sim


class ElectronInputsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # as the Simulator does before filling the dictionary: the C++ side is then
        # torn down in order at exit
        startMPI()

    def tearDown(self):
        # the dictionary holds the python init functions: drop it while python
        # is still alive, as the Simulator does on reset
        ph.clearDict()
        ph.global_vars.sim = None

    def test_numeric_types_are_accepted(self):
        cases = [
            {"closure": "isothermal", "Te": 0},
            {"closure": "isothermal", "Te": 0.12},
            {"closure": "isothermal", "Te": np.int64(1)},
            {"closure": "isothermal", "Te": np.float32(0.12)},
        ]
        for electrons in cases:
            with self.subTest(**{k: repr(v) for k, v in electrons.items()}):
                setup_simulation(**electrons)
                ph.populateDict()

    def test_other_types_are_refused(self):
        for bad in ([0.12], True, None):
            with self.subTest(Te=repr(bad)):
                setup_simulation(closure="isothermal", Te=bad)
                with self.assertRaisesRegex(TypeError, "pressure_closure/Te"):
                    ph.populateDict()


    def test_a_negative_net_charge_is_refused(self):
        setup_simulation(
            density=lambda x: -1.0 + 0 * x, sim_args={"min_initial_density": 0.01}, closure="isothermal", Te=0.1
        )
        with self.assertRaisesRegex(ValueError, "below min_initial_density"):
            ph.populateDict()

    def test_the_density_is_not_checked_unless_asked_for(self):
        def hole(x):
            return 1.0 - 0.999 * np.exp(-(((x - 3.0) / 0.5) ** 2))

        # existing inputs with rarefied regions are unchanged
        setup_simulation(density=hole, closure="isothermal", Te=0.1)
        ph.populateDict()

        setup_simulation(
            density=hole, sim_args={"min_initial_density": 0.01}, closure="isothermal", Te=0.1
        )
        with self.assertRaisesRegex(ValueError, "below min_initial_density"):
            ph.populateDict()


    def test_min_initial_density_must_be_a_finite_non_negative_number(self):
        setup_simulation(sim_args={"min_initial_density": "0.01"}, closure="isothermal", Te=0.1)
        with self.assertRaises(TypeError):
            ph.populateDict()
        for bad in (-1.0, float("inf"), float("nan")):
            with self.subTest(min_initial_density=bad):
                setup_simulation(
                    sim_args={"min_initial_density": bad}, closure="isothermal", Te=0.1
                )
                with self.assertRaises(ValueError):
                    ph.populateDict()
        # 0 disables the check
        setup_simulation(
            density=lambda x: 0.0 * x,
            sim_args={"min_initial_density": 0},
            closure="isothermal",
            Te=0.1,
        )
        ph.populateDict()

    def test_infinite_initial_density_is_refused(self):
        setup_simulation(density=lambda x: np.inf + 0 * x, closure="isothermal", Te=0.1)
        with self.assertRaisesRegex(ValueError, "inf"):
            ph.populateDict()



if __name__ == "__main__":
    unittest.main()
