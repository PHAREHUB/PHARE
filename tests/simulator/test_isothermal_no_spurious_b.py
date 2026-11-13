#!/usr/bin/env python3
#
# Regression test for #1527: with isothermal electrons, the pressure term of Ohm's law is a
# gradient, so it must not create magnetic field. Start from B = 0 and V = 0 with an uneven
# density (plus particle noise), advance, and check that B stays 0 to round-off.
# 2D: in 1D the curl of the pressure term vanishes whatever its discretisation.
#

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator import SimulatorTest

ph.NO_GUI()

diag_dir = "phare_outputs/isothermal_no_spurious_b"
time_step = 0.002
time_step_nbr = 20
final_time = time_step * time_step_nbr


def config():
    def density(x, y):
        return 1.0 + 0.5 * np.exp(-((x - 2.5) ** 2 + (y - 2.0) ** 2) / 1.0)

    def zero(x, y):
        return 0.0 * x

    def vth(x, y):
        return 0.0 * x + 0.3

    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=(20, 20),
        dl=(0.25, 0.25),
        resistivity=0.0,
        hyper_resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )
    ph.MaxwellianFluidModel(
        bx=zero,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": 50,
            "init": {"seed": 1527},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.2)
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=[final_time])
    return sim


def max_abs(h5_path, prefix):
    """largest |value| of each dataset whose name starts with prefix, over all patches;
    raises on any non-finite value, which max() would otherwise hide"""
    import h5py

    out = {}

    def visit(name, item):
        key = name.split("/")[-1]
        if isinstance(item, h5py.Dataset) and key.startswith(prefix):
            values = np.abs(item[()])
            if not np.isfinite(values).all():
                raise ValueError(f"{name}: non-finite values")
            out[key] = max(out.get(key, 0.0), float(values.max()))

    with h5py.File(h5_path, "r") as f:
        f.visititems(visit)
    return out


class IsothermalNoSpuriousB(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def test_b_stays_zero(self):
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config()).run().reset()

        cpp.mpi_barrier()  # before the assertions, so a failure on rank 0 cannot skip it
        if cpp.mpi_rank() == 0:
            E = max_abs(f"{diag_dir}/EM_E.h5", "EM_E")
            B = max_abs(f"{diag_dir}/EM_B.h5", "EM_B")
            maxE = max(E.values())
            self.assertGreater(maxE, 0.05)  # the density drives a pressure term
            # B = 0 up to round-off of E: anything above is field created by the pressure term
            self.assertLess(max(B.values()), 1e-12 * maxE, msg=f"B {B}, E {E}")


if __name__ == "__main__":
    startMPI()
    unittest.main()
