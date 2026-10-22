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


def setup_simulation(density=lambda x: 1.0, sim_args={}, **electrons):
    ph.global_vars.sim = None
    sim = ph.Simulation(
        **sim_args,
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
            dict(closure="isothermal", Te=0),
            dict(closure="isothermal", Te=0.12),
            dict(closure="isothermal", Te=np.int64(1)),
            dict(closure="isothermal", Te=np.float32(0.12)),
            dict(closure="polytropic", Gamma=2, Pe=lambda x: 0.1),
            dict(closure="polytropic", Gamma=np.float64(5 / 3), Pe=lambda x: 0.1),
        ]
        for electrons in cases:
            with self.subTest(**{k: repr(v) for k, v in electrons.items()}):
                setup_simulation(**electrons)
                ph.populateDict()

    def test_other_types_are_refused(self):
        for bad in ([0.12], True, None):
            with self.subTest(Te=repr(bad)):
                setup_simulation(closure="isothermal", Te=bad)
                with self.assertRaisesRegex(ValueError, "pressure_closure/Te"):
                    ph.populateDict()

    def test_too_low_initial_density_is_refused(self):
        def hole(x):
            return 1.0 - 0.999 * np.exp(-(((x - 3.0) / 0.5) ** 2))

        with self.assertRaisesRegex(ValueError, "below min_density"):
            setup_simulation(density=hole, closure="isothermal", Te=0.1)

        # a deliberately low density is accepted when min_density allows it
        setup_simulation(
            density=hole, sim_args={"min_density": 1e-4}, closure="isothermal", Te=0.1
        )

    def test_polytropic_default_gamma_is_five_thirds(self):
        sim = setup_simulation(closure="polytropic")
        self.assertEqual(sim.electrons.closure.Gamma, 5.0 / 3.0)


if __name__ == "__main__":
    unittest.main()
