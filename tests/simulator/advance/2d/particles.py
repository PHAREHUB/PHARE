import unittest
from ddt import ddt, data, unpack
from pyphare.core.box import Box, Box2D, nDBox
from tests.simulator.test_advance import AdvanceTestBase

import matplotlib

matplotlib.use("Agg")  # for systems without GUI

ndim = 2
interp_orders = [1, 2, 3]
ppc = 25

@ddt
class AdvanceTest(AdvanceTestBase):

    # @data(
    #   {"L0": [Box2D(10, 20)]},
    #   {"L0": [Box2D(2, 12), Box2D(13, 25)]},
    # )
    # def test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(self, refinement_boxes):
    #     for interp_order in [1, 2, 3]:
    #         self._test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(ndim, interp_order, refinement_boxes, ppc=ppc)


    # @data(
    #   {"L0": [Box2D(10, 20)]},
    #   {"L0": [Box2D(2, 12), Box2D(13, 25)]},
    # )
    # def test_overlapped_particledatas_have_identical_particles(self, refinement_boxes):
    #     for interp_order in [1, 2, 3]:
    #         self._test_overlapped_particledatas_have_identical_particles(ndim, interp_order, refinement_boxes, ppc=ppc)


    def test_L0_particle_number_conservation(self):
        self._test_L0_particle_number_conservation(ndim, ppc=ppc)


if __name__ == "__main__":
    unittest.main()

"""
======================================================================
FAIL: test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles_1 (__main__.AdvanceTest)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/opt/py/py/lib/python3.9/site-packages/ddt.py", line 182, in wrapper
    return func(self, *args, **kwargs)
  File "/home/philix/git/phare/master/tests/simulator/advance/2d/particles.py", line 23, in test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles
    self._test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles(ndim, interp_order, refinement_boxes, ppc=ppc)
  File "/home/philix/git/phare/master/tests/simulator/test_advance.py", line 317, in _test_patch_ghost_particle_are_clone_of_overlaped_patch_domain_particles
    np.testing.assert_allclose(refdomain.deltas[sort_refdomain_idx], cmpghost.deltas[sort_cmpghost_idx], atol=1e-12)
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 1528, in assert_allclose
    assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 761, in assert_array_compare
    raise AssertionError(msg)
AssertionError:
Not equal to tolerance rtol=1e-07, atol=1e-12

(shapes (500, 2, 2), (50, 2, 2) mismatch)
 x: array([[[0.283109, 0.646739],
        [0.46636 , 0.953576]],
...
 y: array([[[0.786481, 0.307409],
        [0.131742, 0.711136]],
...


======================================================================
FAIL: test_overlapped_particledatas_have_identical_particles_1 (__main__.AdvanceTest)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/opt/py/py/lib/python3.9/site-packages/ddt.py", line 182, in wrapper
    return func(self, *args, **kwargs)
  File "/home/philix/git/phare/master/build/tests/simulator/advance/2d/particles.py", line 32, in test_overlapped_particledatas_have_identical_particles
    self._test_overlapped_particledatas_have_identical_particles(ndim, interp_order, refinement_boxes, ppc=ppc)
  File "/home/philix/git/phare/master/build/tests/simulator/test_advance.py", line 387, in _test_overlapped_particledatas_have_identical_particles
    np.testing.assert_array_equal(part1.iCells[idx1]+offsets[0], part2.iCells[idx2]+offsets[1])
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 932, in assert_array_equal
    assert_array_compare(operator.__eq__, x, y, err_msg=err_msg,
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 842, in assert_array_compare
    raise AssertionError(msg)
AssertionError:
Arrays are not equal

Mismatched elements: 2200 / 4400 (50%)
Max absolute difference: 1
Max relative difference: 0.05
 x: array([[[ 0, 19],
        [ 0, 19]],
...
 y: array([[[ 0, 20],
        [ 0, 20]],
...



"""
