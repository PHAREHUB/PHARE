import unittest

import numpy as np
from ddt import data, ddt

from pyphare.core.phare_utilities import is_scalar, on_step_grid


@ddt
class PhareUtilitiesTest(unittest.TestCase):
    @data(1, 1.0, np.float64(1), np.int32(1), True, None, "str")
    def test_is_scalar(self, arg):
        self.assertTrue(is_scalar(arg))

    @data([1], (1,), {"a": 1}, {1}, np.array([1]), np.array(1), range(1))
    def test_is_not_scalar(self, arg):
        self.assertFalse(is_scalar(arg))

    def test_on_step_grid_rejects_both_sides_of_a_step(self):
        self.assertFalse(on_step_grid([0.26], 0.0, 0.1))
        self.assertFalse(on_step_grid([0.34], 0.0, 0.1))
        self.assertTrue(on_step_grid([0.0, 0.1, 0.30000000000000004], 0.0, 0.1))

    def test_on_step_grid_is_anchored_at_start(self):
        start, step = 0.300004, 0.1
        grid = start + np.arange(11) * step
        self.assertTrue(on_step_grid(grid, start, step))
        self.assertFalse(on_step_grid([0.4, 0.5], start, step))

    def test_on_step_grid_absorbs_float64_rounding(self):
        start, step, nbr = 0.30000000000000004, 1e-9, 10000
        self.assertTrue(on_step_grid([start, start + nbr * step], start, step))

    def test_on_step_grid_absorbs_rounding_of_a_large_start(self):
        start, step = 1048576.0, 1e-5
        self.assertTrue(on_step_grid(start + np.arange(1, 4) * step, start, step))
        self.assertFalse(on_step_grid([start + 0.5 * step], start, step))


if __name__ == "__main__":
    unittest.main()
