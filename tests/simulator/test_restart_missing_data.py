#!/usr/bin/env python3
"""A restart file that lacks data this version needs is refused with a clear error.

Restart files restore only the patch data they contain. A file written before a field was added
(e.g. the ion current on the E nodes, #1577) used to leave that field unallocated and the run
crashed in the first step. Here a restart is written, chargeFluxE is deleted from it, and the
restarted run must stop with the error naming the field, not a segfault.
"""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


def run(out, restart_time=None):
    import pyphare.pharein as ph
    from pyphare.simulator.simulator import Simulator

    ph.NO_GUI()
    ph.global_vars.sim = None
    restart_options = dict(dir=out, mode="overwrite")
    if restart_time is None:
        restart_options["timestamps"] = [0.05]
    else:
        restart_options["restart_time"] = restart_time
    sim = ph.Simulation(
        time_step=0.01, final_time=0.1, cells=(16,), dl=(0.5,), interp_order=1,
        largest_patch_size=8, smallest_patch_size=8, restart_options=restart_options,
    )

    def c(v):
        return lambda x: v + 0.0 * x

    ph.MaxwellianFluidModel(
        bx=c(1.0), by=c(0.0), bz=c(0.0),
        protons={"charge": 1, "mass": 1.0, "density": c(1.0),
                 **{f"vbulk{k}": c(0.0) for k in "xyz"},
                 **{f"vth{k}": c(0.3) for k in "xyz"},
                 "nbr_part_per_cell": 20, "init": {"seed": 1}},
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)
    Simulator(sim, log_to_file=False).run().reset()


def drop_from_restart(out, name):
    """make the file look written by a version without `name`: drop its patch data and remove
    it from each patch's SAMRAI patch_data_namelist"""
    import h5py  # here, not at the top: see doc/conventions.md (imports)

    key = f"{name}##default"
    dropped = 0
    for path in Path(out).rglob("proc.*"):
        with h5py.File(path, "a") as h5:
            patches = []
            h5.visit(lambda n: patches.append(n) if n.endswith("/" + key) else None)
            for data_path in patches:
                patch = h5[data_path.rsplit("/", 1)[0]]
                del patch[key]
                names = [n for n in patch["patch_data_namelist"][()] if n.decode() != key]
                dtype = patch["patch_data_namelist"].dtype
                attrs = dict(patch["patch_data_namelist"].attrs)  # SAMRAI's type marker
                del patch["patch_data_namelist"]
                patch.create_dataset("patch_data_namelist", data=np.array(names, dtype=dtype))
                patch["patch_data_namelist"].attrs.update(attrs)
                patch["patch_data_namelist_count"][0] = len(names)
                dropped += 1
    return dropped


class RestartMissingDataTest(unittest.TestCase):
    def test_a_restart_file_without_a_needed_field_is_refused(self):
        with tempfile.TemporaryDirectory() as out:
            subprocess.run([sys.executable, __file__, "write", out], check=True)
            self.assertGreater(drop_from_restart(out, "chargeFluxE"), 0)
            restarted = subprocess.run([sys.executable, __file__, "restart", out],
                                       capture_output=True, text=True)
            self.assertNotEqual(restarted.returncode, 0)
            self.assertNotEqual(restarted.returncode, -11)  # not a segfault
            self.assertIn("restart file does not contain data", restarted.stdout + restarted.stderr)
            self.assertIn("chargeFluxE", restarted.stdout + restarted.stderr)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "write":
        run(sys.argv[2])
    elif len(sys.argv) == 3 and sys.argv[1] == "restart":
        run(sys.argv[2], restart_time=0.05)
    else:
        unittest.main()
