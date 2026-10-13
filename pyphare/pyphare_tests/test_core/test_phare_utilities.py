import unittest

import numpy as np
from ddt import ddt, data

from pyphare.core.phare_utilities import is_scalar


@ddt
class PhareUtilitiesTest(unittest.TestCase):
    @data(1, 1.0, np.float64(1), np.int32(1), True, None, "str")
    def test_is_scalar(self, arg):
        self.assertTrue(is_scalar(arg))

    @data([1], (1,), {"a": 1}, {1}, np.array([1]), np.array(1), range(1))
    def test_is_not_scalar(self, arg):
        self.assertFalse(is_scalar(arg))


if __name__ == "__main__":
    unittest.main()
