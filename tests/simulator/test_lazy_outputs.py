#!/usr/bin/env python3
"""
The time loop keeps the population densities and E-node fluxes up to date, not the population
fluxes on the moment nodes, the charge flux there, the bulk velocity or the momentum tensors: those
are computed when an output needs them (HybridLazyMoments, #1577). A writer that forgot to compute
them would write the values of the previous output silently.

Check: a population flux diagnostic triggered only by elapsed wall-clock time, never by a simulation
time, writes after an advance the flux the Python getters compute at that time.
"""

import os
import unittest
from time import sleep

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy.fromh5 import h5_filename_from, h5_time_grp_key
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

out = "phare_outputs/lazy_outputs"


def config():
    ph.global_vars.sim = None
    sim = ph.Simulation(
        time_step=0.002,
        time_step_nbr=2,
        cells=[16, 16],
        dl=[0.4, 0.4],
        interp_order=1,
        largest_patch_size=8,
        smallest_patch_size=8,
        refinement_boxes={"L0": {"B0": [(4, 4), (11, 11)]}},
        hyper_resistivity=0.001,
        strict=True,
        diag_options={
            "format": "phareh5",
            "options": {"dir": out, "mode": "overwrite"},
        },
    )

    def constant(value):
        return lambda x, y: value + 0.0 * x

    ph.MaxwellianFluidModel(
        bx=constant(1.0),
        by=constant(0.3),
        bz=constant(0.2),
        protons={
            "charge": 1,
            "density": lambda x, y: 1.0 + 0.3 * np.sin(2 * np.pi * x / 6.4),
            "vbulkx": constant(0.3),
            "vbulky": constant(0.2),
            "vbulkz": constant(-0.1),
            **{f"vth{k}": constant(0.2) for k in "xyz"},
            "nbr_part_per_cell": 50,
            "init": {"seed": 1},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.1)
    # wall-clock triggers only: the first write at the initial dump, the second at the first
    # dump call after an advance
    ph.FluidDiagnostics(
        quantity="flux", population_name="protons", elapsed_timestamps=[0, 0]
    )
    (diag,) = sim.diagnostics.values()
    return sim, diag


class LazyOutputsTest(SimulatorTest):
    def tearDown(self):
        if hasattr(self, "simulator"):
            self.simulator.reset()
        super().tearDown()

    def test_elapsed_time_write_computes_the_flux(self):
        import h5py  # see doc/conventions.md section 2.1.1

        sim, diag = config()
        self.assertEqual(len(diag.write_timestamps), 0)
        self.simulator = Simulator(sim).setup()
        sleep(1)  # so that the wall-clock triggers have elapsed
        self.simulator.initialize()
        self.simulator.advance()

        # the getters compute the flux at the current time themselves
        wrangler = self.simulator.data_wrangler()
        expected = {}
        for lvl in range(wrangler.getNumberOfLevels()):
            level = wrangler.getPatchLevel(lvl)
            for c, getter in (
                ("x", level.getFx),
                ("y", level.getFy),
                ("z", level.getFz),
            ):
                for patch in getter("protons"):
                    lower = tuple(int(v) for v in patch.lower)
                    expected[(lvl, lower, c)] = np.asarray(patch.data)

        h5 = h5py.File(os.path.join(out, h5_filename_from(diag)), "r")
        times = sorted(h5[h5_time_grp_key].keys(), key=float)
        self.assertEqual(len(times), 2, times)  # initial dump, and after the advance
        written = {}
        for lvlgrp, level in h5[h5_time_grp_key][times[-1]].items():
            if not lvlgrp.startswith("pl"):
                continue
            for patch in level.values():
                lower = tuple(int(v) for v in patch.attrs["lower"])
                for name, dataset in patch.items():
                    written[(int(lvlgrp[2:]), lower, name[-1])] = dataset[()]

        self.assertEqual(set(written), set(expected))
        for key, value in written.items():
            # assert_array_equal accepts NaN on both sides: a flux left unset must fail here
            self.assertTrue(np.isfinite(value).all(), str(key))
            np.testing.assert_array_equal(
                value.ravel(),
                expected[key].ravel().astype(value.dtype),
                err_msg=str(key),
            )


if __name__ == "__main__":
    unittest.main()
