#!/usr/bin/env python3

import unittest

import numpy as np
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = 100
dl = 0.2
middle = cells * dl / 2
max_nbr_levels = 3
time_step_nbr = 3
thresholds = {"default": 0.05, "lohner": 0.3}
diag_root = "phare_outputs/tagging_quantities"


def step(low, high):
    def fn(x):
        return np.where(x < middle, low, high)

    return fn


def mhd_simulation(tagging, diag_dir):
    time_step = 0.01
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=dl,
        refinement="tagging",
        max_nbr_levels=max_nbr_levels,
        max_mhd_level=max_nbr_levels,
        tagging=tagging,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        eta=0.0,
        nu=0.0,
        gamma=5.0 / 3.0,
        reconstruction="Linear",
        limiter="VanLeer",
        riemann="Rusanov",
        mhd_timestepper="TVDRK2",
        hall=False,
        model_options=["MHDModel"],
    )
    ph.MHDModel(
        density=step(1.0, 2.0),
        vx=step(0.0, 1.0),
        vy=step(0.0, 0.0),
        vz=step(0.0, 0.0),
        bx=step(1.0, 1.0),
        by=step(1.0, 2.0),
        bz=step(0.0, 0.0),
        p=step(1.0, 2.0),
    )
    ph.ElectromagDiagnostics(
        quantity="B", write_timestamps=[0.0, time_step * time_step_nbr]
    )
    return sim


def hybrid_simulation(tagging, diag_dir):
    time_step = 0.001
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=dl,
        interp_order=1,
        refinement="tagging",
        max_nbr_levels=max_nbr_levels,
        tagging=tagging,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
        hyper_resistivity=0.002,
        resistivity=0.001,
    )
    ph.MaxwellianFluidModel(
        bx=step(1.0, 1.0),
        by=step(1.0, 2.0),
        bz=step(0.0, 0.0),
        protons={
            "charge": 1,
            "density": step(1.0, 1.0),
            "vbulkx": step(0.0, 0.0),
            "vbulky": step(0.0, 0.0),
            "vbulkz": step(0.0, 0.0),
            "vthx": step(0.3, 0.3),
            "vthy": step(0.3, 0.3),
            "vthz": step(0.3, 0.3),
            "nbr_part_per_cell": 20,
            "init": {"seed": 1},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    ph.ElectromagDiagnostics(
        quantity="B", write_timestamps=[0.0, time_step * time_step_nbr]
    )
    return sim


class TaggingQuantitiesTest(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None
        super().tearDown()

    def check_refines_at_startup(self, model, make_simulation, quantities):
        for method, threshold in thresholds.items():
            for quantity in quantities:
                with self.subTest(method=method, quantity=quantity):
                    ph.global_vars.sim = None
                    diag_dir = f"{diag_root}/{model}_{method}_{quantity}"
                    self.register_diag_dir_for_cleanup(diag_dir)
                    tagging = {"method": method, "quantities": {quantity: threshold}}
                    sim = make_simulation(tagging, diag_dir)
                    self.simulator = Simulator(sim)
                    self.simulator.run()
                    self.simulator.reset()
                    self.simulator = None

                    run = Run(diag_dir)
                    times = run.times("B")
                    self.assertEqual(len(times), 2)
                    for time in times:
                        self.assertEqual(
                            run.GetB(time).levelNbr(time), max_nbr_levels, time
                        )

    def test_mhd_allowed_quantities_refine_at_startup(self):
        self.check_refines_at_startup(
            "mhd", mhd_simulation, ("B", "density", "momentum", "Etot")
        )

    def test_hybrid_allowed_quantities_refine_at_startup(self):
        self.check_refines_at_startup("hybrid", hybrid_simulation, ("B", "By"))


if __name__ == "__main__":
    unittest.main()
