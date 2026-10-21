#!/usr/bin/env python3

import numpy as np

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator, startMPI

from tests.simulator import SimulatorTest
from tests.functional.mhd_harris import harris

ph.NO_GUI()

reference_dir = "phare_outputs/mhd_harris_reference"
split_dir = "phare_outputs/mhd_harris_time_varying_b0"

b0 = (1.0, 0.5, 1.0)
omega = 1.0


def external_field():
    def az(x, y, t):
        return (b0[0] * y - b0[1] * x) * np.sin(omega * t)

    def ay(x, y, t):
        return b0[2] * x * np.sin(omega * t)

    def daz_dt(x, y, t):
        return omega * (b0[0] * y - b0[1] * x) * np.cos(omega * t)

    def day_dt(x, y, t):
        return omega * b0[2] * x * np.cos(omega * t)

    return {
        "type": "user-defined",
        "potential": (None, ay, az),
        "potential_time_derivative": (None, day_dt, daz_dt),
    }


def run_harris(diag_dir, external_field=None):
    sim = harris.config(diag_dir=diag_dir, external_field=external_field)
    if external_field is not None:
        ph.ElectromagDiagnostics(quantity="B0", write_timestamps=harris.timestamps)
    Simulator(sim).run().reset()
    ph.global_vars.sim = None
    return sim


def interior(patch_data):
    return patch_data.dataset[
        tuple(slice(g, -g if g else None) for g in patch_data.ghosts_nbr)
    ]


def level0(hierarchy, time):
    return {
        str(patch.box): patch.patch_datas for patch in hierarchy.level(0, time).patches
    }


def max_differences(reference, split, time, quantities):
    ref_patches, split_patches = level0(reference, time), level0(split, time)
    assert ref_patches.keys() == split_patches.keys(), "level 0 patches differ"

    diffs, scale = dict.fromkeys(quantities, 0.0), 0.0
    for box, ref_pds in ref_patches.items():
        for qty in quantities:
            ref = interior(ref_pds[qty])
            diff = np.max(np.abs(interior(split_patches[box][qty]) - ref))
            diffs[qty] = max(diffs[qty], diff)
            scale = max(scale, np.max(np.abs(ref)))
    return {qty: diff / scale if scale > 0 else diff for qty, diff in diffs.items()}


def max_b0_error(run, time):
    hierarchy = run.GetB0(time, all_primal=False)
    error = 0.0
    for patch_datas in level0(hierarchy, time).values():
        for c, b in zip("xyz", b0):
            expected = b * np.sin(omega * time)
            error = max(
                error, np.max(np.abs(interior(patch_datas[f"B0{c}"]) - expected))
            )
    return error


fields = {
    "GetMHDrho": ["mhdRho"],
    "GetMHDP": ["mhdP"],
    "GetMHDV": ["mhdVx", "mhdVy", "mhdVz"],
    "GetB": ["Bx", "By", "Bz"],
}

tolerances = {
    "mhdRho": 5e-2,
    "mhdP": 5e-2,
    **dict.fromkeys(["mhdVx", "mhdVy", "mhdVz"], 2.5e-1),
    **dict.fromkeys(["Bx", "By", "Bz"], 5e-2),
}


def compare(reference, split):
    differences = {}
    for time in harris.timestamps:
        for getter, quantities in fields.items():
            ref = getattr(reference, getter)(time, all_primal=False)
            spl = getattr(split, getter)(time, all_primal=False)
            for qty, diff in max_differences(ref, spl, time, quantities).items():
                differences[(time, qty)] = diff
    return differences


class HarrisTimeVaryingB0Test(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def test_run(self):
        self.register_diag_dir_for_cleanup(reference_dir)
        self.register_diag_dir_for_cleanup(split_dir)

        sim = run_harris(reference_dir)
        run_harris(split_dir, external_field())

        b0_errors, differences = {}, {}
        if not sim.dry_run and cpp.mpi_rank() == 0:
            reference, split = Run(reference_dir), Run(split_dir)
            b0_errors = {t: max_b0_error(split, t) for t in harris.timestamps}
            differences = compare(reference, split)
            for (time, qty), difference in differences.items():
                print(
                    f"t={time:5.1f} {qty:7s} max relative difference {difference:.3e}"
                )
        cpp.mpi_barrier()

        for time, error in b0_errors.items():
            self.assertLess(error, 1e-6, f"B0 at t={time}")
        for (time, qty), difference in differences.items():
            self.assertLess(difference, tolerances[qty], f"{qty} at t={time}")
        return self


if __name__ == "__main__":
    startMPI()
    HarrisTimeVaryingB0Test().test_run().tearDown()
