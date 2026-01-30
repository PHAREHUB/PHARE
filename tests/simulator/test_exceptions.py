#!/usr/bin/env python3

import unittest

import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator import simulator


ph.NO_GUI()


def setup_model(ppc=100):
    def density(*xyz):
        return 1.0

    def by(*xyz):
        return 1

    def bz(*xyz):
        return 1

    def bx(*xyz):
        return 1.0

    def vx(*xyz):
        return 1000

    def vy(*xyz):
        return 1000

    def vz(*xyz):
        return 1000

    def vthx(*xyz):
        return 0.01

    def vthy(*xyz):
        return 0.01

    def vthz(*xyz):
        return 0.01

    vvv = {
        "vbulkx": vx,
        "vbulky": vy,
        "vbulkz": vz,
        "vthx": vthx,
        "vthy": vthy,
        "vthz": vthz,
    }

    model = ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "mass": 1,
            "charge": 1,
            "density": density,
            **vvv,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1337},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.12)

    ph.FluidDiagnostics(  # NO TIMESTAMPS!
        quantity="density", write_timestamps=[], population_name="protons"
    )
    return model


out = "phare_outputs/boris_move_test/"
simArgs = {
    "time_step_nbr": 1,
    "final_time": 0.001,
    "boundary_types": "periodic",
    "cells": 20,
    "dl": 0.3,
    "diag_options": {
        "format": "phareh5",
        "options": {"dir": out, "mode": "overwrite", "allow_emergency_dumps": True},
    },
}


class BorisTwoCellJumpTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(BorisTwoCellJumpTest, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None

    def test(self):
        simulation = ph.Simulation(**simArgs)
        setup_model()
        with self.assertRaises(Exception):
            simulator.exit_on_exception = False
            simulator.Simulator(simulation).run()

        # the following will fail if there is no emergency diag
        Run(out).GetN(0, "protons").plot(filename="exceptional_pop_density.png")


if __name__ == "__main__":
    simulator.startMPI()
    unittest.main()
