#!/usr/bin/env python3
"""
FullyKineticModel restarts (P5d, T17): T15's moving light pulses (two levels, level 1 regridded
every two coarse steps) in one run of STEPS coarse steps, against the same run restarted from its
restart file at step N, for N in RESTART_STEPS, boundary_cleaning_mask off and on. 1D, T15's
command line:
    1 1 --cleaning 0.5

A restart does not re-create the levels: SAMRAI reads the hierarchy and every patch data the
resources manager registered (fields, particles, the solver's and the messenger's scratch, their
time labels) from the file (TimeRefinementIntegrator::initializeHierarchy), no level initializer
runs. The first step rewrites the solver's and the messenger's other records before reading them,
except the cells the finer level synchronized at the end of the last step (the cleaning mask).

Check: both runs write a restart file at the last step, and every patch data array in them
(doubles, the particles included) is bitwise the same, on every rank. Two restart steps, one on
each side of a level 1 regrid: the mask must know the cells of the level 1 that synchronized,
which a regrid at the end of step N has already replaced.
"""

import glob
import os
import shutil
import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

# parses the command line
from tests.simulator import test_fully_kinetic as fk
from tests.simulator import test_fully_kinetic_amr_regrid as t15

ph.NO_GUI()

RESTART_STEPS = (t15.STEPS // 2 - 1, t15.STEPS // 2)
MASKS = [False, True] if fk.cleaning > 0 else [False]

fk.diag_dir = (
    f"phare_outputs/fully_kinetic_amr_restart_{t15.ndim}d_o{t15.interp_order}"
    f"_filter{fk.filter_passes}_cleaning{fk.cleaning}"
)


def restarts_dir(subdir):
    return f"{fk.diag_dir}/{subdir}/restarts"


def time_dir(subdir, time):
    """Hierarchy::restartFilePathForTime: five decimals, eleven characters"""
    return f"{restarts_dir(subdir)}/{time:011.5f}"


def restart_files(subdir, time):
    """{proc.NNNNNNN: path}, the SAMRAI restart file of every rank"""
    paths = glob.glob(f"{time_dir(subdir, time)}/restore.*/nodes.*/proc.*")
    return {os.path.basename(path): path for path in paths}


def patch_arrays(path):
    """{name: bytes} of every patch data dataset of a SAMRAI restart file"""
    import h5py

    arrays = {}

    def keep(name, obj):
        if isinstance(obj, h5py.Dataset) and "-patch_" in name:
            arrays[name] = np.asarray(obj[()]).tobytes()

    with h5py.File(path, "r") as f:
        f["PHARE_hierarchy"].visititems(keep)
    return arrays


class FullyKineticAMRRestartTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None

    def run_case(self, subdir, mask, restart_step=None):
        ph.global_vars.sim = None
        last = t15.times[-1]
        if restart_step is None:
            simulation = dict(
                restart_options={
                    "dir": restarts_dir(subdir),
                    "mode": "overwrite",
                    "timestamps": [n * fk.time_step for n in RESTART_STEPS] + [last],
                }
            )
            diag_times = None
        else:
            start = restart_step * fk.time_step
            simulation = dict(
                time_step_nbr=t15.STEPS - restart_step,
                restart_options={
                    "dir": restarts_dir(subdir),
                    "mode": "overwrite",
                    "timestamps": [last],
                    "restart_time": start,
                },
            )
            diag_times = [last]
        sim = t15.config(
            True,
            subdir=subdir,
            diag_times=diag_times,
            boundary_cleaning_mask=mask,
            **simulation,
        )
        Simulator(sim).run().reset()

    def test_restart(self):
        if t15.ndim != 1:
            self.skipTest("T17 is 1D")
        for mask in MASKS:
            reference = f"reference_mask{int(mask)}"
            self.run_case(reference, mask)
            for step in RESTART_STEPS:
                with self.subTest(mask=mask, restart_step=step):
                    self._restarted(reference, mask, step)

    def _restarted(self, reference, mask, step):
        restarted = f"restarted{step}_mask{int(mask)}"
        start = step * fk.time_step
        # the restarted run reads and writes its own directory
        if cpp.mpi_rank() == 0:
            shutil.rmtree(restarts_dir(restarted), ignore_errors=True)
            shutil.copytree(time_dir(reference, start), time_dir(restarted, start))
        cpp.mpi_barrier()
        self.run_case(restarted, mask, restart_step=step)
        if cpp.mpi_rank() != 0:
            return

        last = t15.times[-1]
        ref, new = restart_files(reference, last), restart_files(restarted, last)
        self.assertEqual(len(ref), cpp.mpi_size())
        self.assertEqual(ref.keys(), new.keys())
        for proc in sorted(ref):
            a, b = patch_arrays(ref[proc]), patch_arrays(new[proc])
            self.assertEqual(a.keys(), b.keys(), proc)
            differ = sorted(name for name in a if a[name] != b[name])
            print(
                f"mask {mask}, restart at step {step}, {proc}: {len(a)} patch data arrays,"
                f" {len(differ)} differ"
            )
            self.assertEqual(differ, [], f"{proc}: first {differ[:5]}")


if __name__ == "__main__":
    unittest.main()
