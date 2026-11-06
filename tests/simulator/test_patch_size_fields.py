#!/usr/bin/env python3
"""
With a seed, the initial particles of a cell do not depend on the patches (see test_seed_tiling.py),
so neither should the moments and E at t=0. Patches of width G+1, G the field ghost width, put
two patches G+1 cells apart, which corrupted the total ion flux and E at orders 2 and 3 (#1575).
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

# divisible by both patch sizes, G + 1 and G + 2, so that all patches have that width
cells = 30
diag_dir = "phare_outputs/patch_size_fields"
ghost_width = {2: 4, 3: 4}  # field ghost width per interp order (hybrid_options.hpp)
ppc = {1: 128, 2: 32, 3: 8}  # per dimension

FIELDS = {
    "E": [("EM_E.h5", f"EM_E_{c}") for c in "xyz"],
    "B": [("EM_B.h5", f"EM_B_{c}") for c in "xyz"],
    "n": [("ions_pop_protons_density.h5", "density")],
    "population flux": [("ions_pop_protons_flux.h5", f"flux_{c}") for c in "xyz"],
    "total flux": [("ions_charge_flux.h5", f"charge_flux_{c}") for c in "xyz"],
}


def config(ndim, interp, patch_size, out):
    sim = ph.Simulation(
        cells=[cells] * ndim,
        dl=[0.125] * ndim,
        time_step=0.00125,
        time_step_nbr=1,
        interp_order=interp,
        smallest_patch_size=patch_size,
        largest_patch_size=patch_size,
        resistivity=0.0,
        hyper_resistivity=0.0,
        diag_options={
            "format": "phareh5",
            "options": {"dir": out, "mode": "overwrite"},
        },
    )

    def const(value):
        return lambda *xyz: value + np.zeros_like(xyz[0])

    ph.MaxwellianFluidModel(
        bx=const(0.0),
        by=const(0.0),
        bz=const(1.0),
        protons={
            "charge": 1,
            "mass": 1,
            "density": const(1.0),
            "vbulkx": const(0.0),
            "vbulky": const(0.0),
            "vbulkz": const(0.0),
            "vthx": const(0.3),
            "vthy": const(0.3),
            "vthz": const(0.3),
            "nbr_part_per_cell": ppc[ndim],
            "init": {"seed": 11},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    times = [0.0]
    for quantity in ("E", "B"):
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=times)
    for quantity in ("density", "flux"):
        ph.FluidDiagnostics(
            quantity=quantity, population_name="protons", write_timestamps=times
        )
    ph.FluidDiagnostics(quantity="charge_flux", write_timestamps=times)
    return sim


def domain(out, ndim, filename, key):
    """level-0 field at t=0 on the domain, ghosts removed, each node taken once, and the
    patch lower bounds"""
    import h5py

    field, lowers = np.full([cells] * ndim, np.nan), []
    with h5py.File(f"{out}/{filename}") as h5:
        for patch in h5["t"]["0.0000000000"]["pl0"].values():
            data = patch[key]
            ghosts = int(data.attrs["ghosts"])
            lower = [int(v) for v in patch.attrs["lower"]]
            ncells = [int(v) for v in patch.attrs["nbrCells"]]
            dst = tuple(slice(lo, lo + n) for lo, n in zip(lower, ncells))
            src = tuple(slice(ghosts, ghosts + n) for n in ncells)
            field[dst] = data[...].astype(float)[src]
            lowers.append(tuple(lower))
    return field, sorted(lowers)


class PatchSizeFieldsTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_with(self, ndim, interp, patch_size):
        out = f"{diag_dir}/{ndim}d/interp{interp}/patch{patch_size}"
        ph.global_vars.sim = None
        Simulator(config(ndim, interp, patch_size, out)).initialize().reset()
        return out

    def check(self, ndim):
        for interp, G in ghost_width.items():
            reference = self.run_with(ndim, interp, cells)
            # G + 1 is the gap of #1575, G + 2 a control
            for patch_size in (G + 1, G + 2):
                with self.subTest(ndim=ndim, interp=interp, patch_size=patch_size):
                    out = self.run_with(ndim, interp, patch_size)
                    _, ref_lowers = domain(reference, ndim, *FIELDS["E"][0])
                    _, lowers = domain(out, ndim, *FIELDS["E"][0])
                    self.assertEqual(ref_lowers, [(0,) * ndim])
                    self.assertEqual(len(lowers), (cells // patch_size) ** ndim)
                    for name, fields in FIELDS.items():
                        for filename, key in fields:
                            np.testing.assert_allclose(
                                domain(out, ndim, filename, key)[0],
                                domain(reference, ndim, filename, key)[0],
                                rtol=0,
                                atol=1e-12,
                                err_msg=f"{name} {key}",
                            )

    def test_1d(self):
        self.check(1)

    def test_2d(self):
        self.check(2)

    def test_3d(self):
        self.check(3)


if __name__ == "__main__":
    unittest.main()
