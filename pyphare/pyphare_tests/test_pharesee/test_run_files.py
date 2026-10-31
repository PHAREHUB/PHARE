"""
Run reads times from the files it finds, keyed by file stem, without importing h5py
until a file is read, and hierarchy_compare reports what differs
"""

import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from pyphare.pharesee.hierarchy.hierarchy_utils import hierarchy_compare
from pyphare.pharesee.run import Run
from pyphare_tests.test_pharesee import build_hierarchy


def write_times(path, times):
    import h5py  # see doc/conventions.md section 2.1.1

    with h5py.File(path, "w") as f:
        for time in times:
            f.create_group(f"t/{time:.10f}")


def hierarchy(refinement_boxes):
    return build_hierarchy(
        nbr_cells=np.asarray([20, 20]),
        origin=np.zeros(2),
        interp_order=1,
        domain_size=np.ones(2),
        cell_width=np.ones(2) / 20,
        refinement_boxes=refinement_boxes,
        quantities=["Bx", "By"],
    )


class RunFilesTest(unittest.TestCase):
    def test_import_does_not_load_h5py(self):
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path))
        code = "import sys, pyphare.pharesee.run; sys.exit('h5py' in sys.modules)"
        subprocess.run([sys.executable, "-c", code], env=env, check=True)

    def test_times_keyed_by_file_stem(self):
        with tempfile.TemporaryDirectory() as path:
            write_times(os.path.join(path, "EM_B.h5"), [0.1, 0.0])
            write_times(os.path.join(path, "ions_charge_density.h5"), [0.2])

            run = Run(path)
            times = run.all_times()
            self.assertEqual(sorted(times), ["EM_B", "ions_charge_density"])
            np.testing.assert_array_equal(run.times("EM_B"), [0.0, 0.1])  # sorted
            np.testing.assert_array_equal(run.times("ions_charge_density"), [0.2])
            self.assertIn(run.default_time, (0.0, 0.2))  # first time of some file
            self.assertEqual(Run(path, default_time=0.1).default_time, 0.1)

    def test_no_files(self):
        with tempfile.TemporaryDirectory() as path:
            with self.assertRaises(RuntimeError) as error:
                Run(path)
            self.assertIn("No HDF5 files", str(error.exception))

    def test_hierarchy_compare_reports(self):
        boxes = {"L0": {"B0": [(5, 5), (9, 9)]}}
        ref = hierarchy(boxes)

        same = hierarchy_compare(ref, hierarchy(boxes))
        self.assertTrue(same)
        self.assertEqual(repr(same), "==")

        other_box = hierarchy_compare(ref, hierarchy({"L0": {"B0": [(4, 4), (9, 9)]}}))
        self.assertFalse(other_box)
        self.assertEqual(repr(other_box), "patch box mismatch")

        cmp = hierarchy(boxes)
        cmp.level(1).patches[0]["By"].dataset[3, 3] += 1
        other_data = hierarchy_compare(ref, cmp)
        self.assertFalse(other_data)
        self.assertIn("data mismatch: FieldData By", repr(other_data))


if __name__ == "__main__":
    unittest.main()
