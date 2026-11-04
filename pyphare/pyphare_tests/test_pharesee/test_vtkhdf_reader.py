"""
reading vtkhdf files PHARE writes, including the bad ones seen on real runs: a time listed
but not written (NSteps), a time dumped again after a restart, no interpOrder (MHD),
offsets saturated by 32 bit types
"""

import os
import tempfile
import unittest

import numpy as np

from pyphare.core.box import Box
from pyphare.pharesee.hierarchy.for_vtk.patchdata import VtkFieldDatasetAccessor
from pyphare.pharesee.hierarchy.fromvtkhdf5 import (
    VtkFile,
    get_times_from_h5,
    hierarchy_fromvtkhdf,
)

cells = 10
rows_per_step = (cells + 1) * 4  # primal nodes, written 4 times in 1d


def write_file(path, times, nsteps=None, interp_order=1, data_offsets=None):
    """one 1d patch per step, every row of step i holds the value i"""
    import h5py  # see doc/conventions.md section 2.1.1

    nvalues = len(times)
    nsteps = nvalues if nsteps is None else nsteps
    with h5py.File(path, "w") as f:
        f.attrs["dimension"] = 1
        f.attrs["domain_box"] = [cells - 1]
        if interp_order is not None:
            f.attrs["interpOrder"] = interp_order
        root = f.create_group("VTKHDF")
        steps = root.create_group("Steps")
        steps.attrs["NSteps"] = nsteps
        steps.create_dataset("Values", data=np.asarray(times, dtype=np.float64))
        level = root.create_group("Level0")
        level.attrs["Spacing"] = np.asarray(
            [0.1], dtype=np.float32
        )  # one per dim, as written
        level.create_dataset(
            "AMRBox", data=np.asarray([[0, cells - 1, 0, 0, 0, 0]] * nvalues)
        )
        data = np.concatenate([np.full(rows_per_step, i) for i in range(nvalues)])
        level.create_group("PointData").create_dataset("data", data=data[:, None])
        step_level = steps.create_group("Level0")
        step_level.create_dataset("NumberOfAMRBox", data=np.ones(nvalues, dtype=int))
        step_level.create_dataset("AMRBoxOffset", data=np.arange(nvalues))
        if data_offsets is None:
            data_offsets = np.arange(nvalues) * rows_per_step
        step_level.create_group("PointDataOffset").create_dataset(
            "data", data=np.asarray(data_offsets)
        )


def write_two_patches(path):
    """one step, two patches of different sizes, rows of patch j hold j + 1"""
    import h5py  # see doc/conventions.md section 2.1.1

    boxes = [(0, 3), (4, 9)]  # 5 and 7 primal nodes
    with h5py.File(path, "w") as f:
        f.attrs["dimension"] = 1
        f.attrs["domain_box"] = [cells - 1]
        f.attrs["interpOrder"] = 1
        steps = f.create_group("VTKHDF").create_group("Steps")
        steps.attrs["NSteps"] = 1
        steps.create_dataset("Values", data=[0.0])
        level = f["VTKHDF"].create_group("Level0")
        level.attrs["Spacing"] = np.asarray([0.1], dtype=np.float32)
        level.create_dataset("AMRBox", data=[[lo, up, 0, 0, 0, 0] for lo, up in boxes])
        data = [np.full((up - lo + 2) * 4, j + 1) for j, (lo, up) in enumerate(boxes)]
        level.create_group("PointData").create_dataset(
            "data", data=np.concatenate(data)[:, None]
        )
        step_level = steps.create_group("Level0")
        step_level.create_dataset("NumberOfAMRBox", data=[len(boxes)])
        step_level.create_dataset("AMRBoxOffset", data=[0])
        step_level.create_group("PointDataOffset").create_dataset("data", data=[0])


class CountingDataset:
    """stands in for an h5py dataset, counts reads"""

    def __init__(self, array):
        self.array, self.reads = array, 0

    def __getitem__(self, item):
        self.reads += 1
        return self.array[item]


class VtkhdfReaderTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "ions_charge_density.vtkhdf")

    def tearDown(self):
        self.tmp.cleanup()

    def test_times_stop_at_nsteps(self):
        write_file(self.path, [0.0, 0.1, 0.2], nsteps=2)  # last dump never written
        vtk_file = VtkFile(self.path)
        np.testing.assert_array_equal(vtk_file.times(), [0.0, 0.1])
        self.assertFalse(vtk_file.has_time(0.2))
        np.testing.assert_array_equal(get_times_from_h5(self.path), [0.0, 0.1])

    def test_restart_duplicate_reads_last_copy(self):
        write_file(self.path, [0.0, 0.1, 0.2, 0.2, 0.3])  # restarted at 0.2
        vtk_file = VtkFile(self.path)
        np.testing.assert_array_equal(vtk_file.times(), [0.0, 0.1, 0.2, 0.3])
        self.assertEqual(vtk_file.time_idx(0.2), 3)

        hier = hierarchy_fromvtkhdf(self.path)  # all times
        self.assertEqual(len(hier.times()), 4)
        patch = hier.level(0, time=0.2).patches[0]
        np.testing.assert_array_equal(np.unique(patch["rho"][:]), [3])

    def test_no_interp_order(self):
        write_file(self.path, [0.0], interp_order=None)  # MHD files
        self.assertEqual(VtkFile(self.path).interp_order, 0)
        hierarchy_fromvtkhdf(self.path, 0.0)

    def test_saturated_offsets_raise(self):
        saturated = [0, rows_per_step, rows_per_step, rows_per_step]
        write_file(self.path, [0.0, 0.1, 0.2, 0.3], data_offsets=saturated)
        vtk_file = VtkFile(self.path)
        with self.assertRaises(ValueError) as error:
            vtk_file.patches(0, vtk_file.time_idx(0.0))
        self.assertIn("corrupt Steps/Level0 offsets", str(error.exception))

    def test_patch_offsets(self):
        write_two_patches(self.path)
        boxes, offsets = VtkFile(self.path).patches(0, 0)
        self.assertEqual([(b.lower[0], b.upper[0]) for b in boxes], [(0, 3), (4, 9)])
        np.testing.assert_array_equal(offsets, [0, 5 * 4])
        patches = hierarchy_fromvtkhdf(self.path, 0.0).level(0, time=0.0).patches
        for j, patch in enumerate(patches):
            np.testing.assert_array_equal(np.unique(patch["rho"][:]), [j + 1])

    def test_patch_read_once(self):
        dataset = CountingDataset(np.arange(10.0)[:, None])
        accessor = VtkFieldDatasetAccessor(dataset, 0, 2, Box([0], [3]))
        np.testing.assert_array_equal(accessor[:], [2, 3, 4, 5])
        np.testing.assert_array_equal(accessor[1:3], [3, 4])
        np.testing.assert_array_equal(np.asarray(accessor), [2, 3, 4, 5])
        self.assertEqual(dataset.reads, 1)
        accessor[:][0] = -1  # callers get copies
        self.assertEqual(accessor[0], 2)

    def test_unknown_time(self):
        write_file(self.path, [0.0])
        with self.assertRaises(IndexError):
            VtkFile(self.path).time_idx(1.0)


if __name__ == "__main__":
    unittest.main()
