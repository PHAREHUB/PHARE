#!/usr/bin/env python3
"""
Every level-0 patch holding a given E or B node, ghosts included, must hold the same value.

The initial condition below is a magnetic bump off the centre of a periodic box, so it is
not periodic: a ghost across the periodic boundary, evaluated at its own coordinate, differs
from the domain node it stands for, and so do the two copies of a face on the periodic
boundary. Hyper-resistivity makes Ohm's law read the current on the first ghost layer,
which is the curl of the magnetic ghosts.
"""

import unittest

import numpy as np
from ddt import data, ddt, unpack

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.core import box as boxm
from pyphare.pharesee.geometry import hierarchy_overlaps
from pyphare.pharesee.hierarchy import hierarchy_from
from pyphare.simulator.simulator import Simulator, startMPI
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = 32
dl = 0.2
time_step = 0.001
time_step_nbr = 4
center, width = 2.0, 0.8


def config(ndim, interp, hyper_mode, diag_dir):
    sim = ph.Simulation(
        cells=[cells] * ndim,
        dl=[dl] * ndim,
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        interp_order=interp,
        smallest_patch_size=8,
        largest_patch_size=16 if ndim > 1 else 8,
        resistivity=0.0,
        hyper_resistivity=0.001,
        hyper_mode=hyper_mode,
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )

    def bump(*xyz):
        return np.exp(-sum((x - center) ** 2 for x in xyz) / width**2)

    def const(value):
        return lambda *xyz: value + 0.0 * xyz[0]

    ph.MaxwellianFluidModel(
        bx=lambda *xyz: 1.0 + bump(*xyz),
        by=bump,
        bz=bump,
        protons={
            "charge": 1,
            "density": const(1.0),
            "vbulkx": const(2.0),
            "vbulky": const(1.0),
            "vbulkz": const(0.0),
            "vthx": const(0.3),
            "vthy": const(0.3),
            "vthz": const(0.3),
            "nbr_part_per_cell": 20,
            "init": {"seed": 1},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    timestamps = time_step * np.arange(time_step_nbr + 1)
    for quantity in ("E", "B"):
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    return sim


@ddt
class HybridSharedFieldsTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    @data(
        *[(ndim, interp, "constant") for ndim in (1, 2) for interp in (1, 2, 3)],
        (2, 3, "spatial"),
    )
    @unpack
    def test_level0_E_and_B_copies_agree(self, ndim, interp, hyper_mode):
        diag_dir = f"phare_outputs/hybrid_shared_fields/{ndim}_{interp}_{hyper_mode}"
        self.register_diag_dir_for_cleanup(diag_dir)
        Simulator(config(ndim, interp, hyper_mode, diag_dir)).run().reset()
        if cpp.mpi_rank() > 0:
            return

        for quantity in ("E", "B"):
            hier = hierarchy_from(h5_filename=f"{diag_dir}/EM_{quantity}.h5")
            for time in time_step * np.arange(time_step_nbr + 1):
                overlaps = hierarchy_overlaps(hier, time)[0]
                self.assertGreater(len(overlaps), 0)
                for overlap in overlaps:
                    box, offsets = overlap["box"], overlap["offset"]
                    values = [
                        boxm.select(
                            pd.dataset,
                            boxm.amr_to_local(box, boxm.shift(pd.ghost_box, offset)),
                        )
                        for pd, offset in zip(overlap["pdatas"], offsets)
                    ]
                    np.testing.assert_array_equal(
                        *values,
                        err_msg=f"{overlap['name']} at t={time} on {box}",
                    )


if __name__ == "__main__":
    startMPI()
    unittest.main()
