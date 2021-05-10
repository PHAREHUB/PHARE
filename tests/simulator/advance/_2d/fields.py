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

    @data(
      {"L0": [Box2D(10, 19)]},
      {"L0": [Box2D(8, 20)]},
    )
    def test_overlaped_fields_are_equal(self, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr=3
        time_step=0.001
        diag_outputs=f"phare_overlaped_fields_are_equal_{ndim}_{self.ddt_test_id()}"
        for interp_order in [1, 2, 3]:
            datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                      time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim, nbr_part_per_cell=ppc)
            self._test_overlaped_fields_are_equal(time_step, time_step_nbr, datahier)


    @data({"L0": [Box2D(10, 19)]})
    def test_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts(self, refinement_boxes):
        print(f"{self._testMethodName}_{ndim}d")
        time_step_nbr=3
        time_step=0.001
        from pyphare.pharein.simulation import check_patch_size
        diag_outputs=f"phare_overlaped_fields_are_equal_with_min_max_patch_size_of_max_ghosts_{ndim}_{self.ddt_test_id()}"
        for interp_order in [1, 2, 3]:
            largest_patch_size, smallest_patch_size = check_patch_size(ndim, interp_order=interp_order, cells=[30] * ndim)
            datahier = self.getHierarchy(interp_order, refinement_boxes, "eb", diag_outputs=diag_outputs,
                                      smallest_patch_size=smallest_patch_size, largest_patch_size=smallest_patch_size,
                                      time_step=time_step, time_step_nbr=time_step_nbr, ndim=ndim, nbr_part_per_cell=ppc)
            self._test_overlaped_fields_are_equal(time_step, time_step_nbr, datahier)


    # @data( # see below
    #    ({"L0": {"B0": Box2D(10, 19)}}),
    #    ({"L0": {"B0": Box2D(10, 14), "B1": Box2D(15, 19)}}),
    #    ({"L0": {"B0": Box2D(6, 23)}}),
    #    ({"L0": {"B0": Box2D( 2, 12), "B1": Box2D(13, 25)}}),
    #    ({"L0": {"B0": Box2D( 5, 20)}, "L1": {"B0": Box2D(15, 19)}}),
    #    # ({"L0": {"B0": Box2D( 5, 20)}, "L1": {"B0": Box2D(12, 38)}, "L2": {"B0": Box2D(30, 52)} }), # particle cell jump > 1
    # )
    # def test_field_coarsening_via_subcycles(self, refinement_boxes):
    #     print(f"{self._testMethodName}_{ndim}d")
    #     for interp_order in [1, 2, 3]:
    #         self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)


    ## needs update in test_refine_field.py::refine
    # @data( # only supports a hierarchy with 2 levels
    #    ({"L0": [Box2D(5, 9)]}),
    #    ({"L0": [Box2D(5, 24)]}),
    #    ({"L0": [Box2D(5, 9), Box2D(20, 24)]}),
    # )
    # def test_field_level_ghosts_via_subcycles_and_coarser_interpolation(self, refinement_boxes):
    #     for interp in [1, 2, 3]:
    #         self._test_field_level_ghosts_via_subcycles_and_coarser_interpolation(ndim, interp, refinement_boxes)


if __name__ == "__main__":
    try:
        from concurrencytest import ConcurrentTestSuite, fork_for_tests
        unittest.TextTestRunner().run(ConcurrentTestSuite(unittest.TestLoader().loadTestsFromTestCase(AdvanceTest), fork_for_tests(5)))
    except ImportError as err:
        unittest.main()



"""
======================================================================
FAIL: test_field_coarsening_via_subcycles_1 (__main__.AdvanceTest)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/opt/py/py/lib/python3.9/site-packages/ddt.py", line 182, in wrapper
    return func(self, *args, **kwargs)
  File "/home/philix/git/phare/master/build/tests/simulator/advance/2d/fields.py", line 58, in test_field_coarsening_via_subcycles
    self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)
  File "/home/philix/git/phare/master/build/tests/simulator/test_advance.py", line 496, in _test_field_coarsening_via_subcycles
    np.testing.assert_allclose(coarse_pdDataset, afterCoarse, atol=1e-6)
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 1528, in assert_allclose
    assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 842, in assert_array_compare
    raise AssertionError(msg)
AssertionError:
Not equal to tolerance rtol=1e-07, atol=1e-06

Mismatched elements: 10 / 1640 (0.61%)
Max absolute difference: 0.24276234
Max relative difference: 6.5977116e-07
 x: array([-0.005596,  0.034836, -0.052158, ..., -0.117555, -0.05284 ,
       -0.005419], dtype=float32)
 y: array([-0.005596,  0.034836, -0.052158, ..., -0.117555, -0.05284 ,
       -0.005419], dtype=float32)

======================================================================
FAIL: test_field_coarsening_via_subcycles_2 (__main__.AdvanceTest)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/opt/py/py/lib/python3.9/site-packages/ddt.py", line 182, in wrapper
    return func(self, *args, **kwargs)
  File "/home/philix/git/phare/master/build/tests/simulator/advance/2d/fields.py", line 58, in test_field_coarsening_via_subcycles
    self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)
  File "/home/philix/git/phare/master/build/tests/simulator/test_advance.py", line 496, in _test_field_coarsening_via_subcycles
    np.testing.assert_allclose(coarse_pdDataset, afterCoarse, atol=1e-6)
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 1528, in assert_allclose
    assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 842, in assert_array_compare
    raise AssertionError(msg)
AssertionError:
Not equal to tolerance rtol=1e-07, atol=1e-06

Mismatched elements: 5 / 1640 (0.305%)
Max absolute difference: 0.10364687
Max relative difference: 1.4190681e-07
 x: array([-0.016569,  0.060141, -0.079471, ..., -0.165187, -0.072189,
        0.050054], dtype=float32)
 y: array([-0.016569,  0.060141, -0.079471, ..., -0.165187, -0.072189,
        0.050054], dtype=float32)

======================================================================
FAIL: test_field_coarsening_via_subcycles_3 (__main__.AdvanceTest)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/opt/py/py/lib/python3.9/site-packages/ddt.py", line 182, in wrapper
    return func(self, *args, **kwargs)
  File "/home/philix/git/phare/master/build/tests/simulator/advance/2d/fields.py", line 58, in test_field_coarsening_via_subcycles
    self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)
  File "/home/philix/git/phare/master/build/tests/simulator/test_advance.py", line 496, in _test_field_coarsening_via_subcycles
    np.testing.assert_allclose(coarse_pdDataset, afterCoarse, atol=1e-6)
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 1528, in assert_allclose
    assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 842, in assert_array_compare
    raise AssertionError(msg)
AssertionError:
Not equal to tolerance rtol=1e-07, atol=1e-06

Mismatched elements: 9 / 1640 (0.549%)
Max absolute difference: 0.0918031
Max relative difference: 4.2922574e-07
 x: array([-0.033605, -0.144496, -0.003578, ...,  0.207218,  0.200839,
        0.067077], dtype=float32)
 y: array([-0.033605, -0.144496, -0.003578, ...,  0.207218,  0.200839,
        0.067077], dtype=float32)

======================================================================
FAIL: test_field_coarsening_via_subcycles_4 (__main__.AdvanceTest)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/opt/py/py/lib/python3.9/site-packages/ddt.py", line 182, in wrapper
    return func(self, *args, **kwargs)
  File "/home/philix/git/phare/master/build/tests/simulator/advance/2d/fields.py", line 58, in test_field_coarsening_via_subcycles
    self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)
  File "/home/philix/git/phare/master/build/tests/simulator/test_advance.py", line 496, in _test_field_coarsening_via_subcycles
    np.testing.assert_allclose(coarse_pdDataset, afterCoarse, atol=1e-6)
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 1528, in assert_allclose
    assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 842, in assert_array_compare
    raise AssertionError(msg)
AssertionError:
Not equal to tolerance rtol=1e-07, atol=1e-06

Mismatched elements: 11 / 1640 (0.671%)
Max absolute difference: 0.18894267
Max relative difference: 3.259992e-07
 x: array([ 0.038159, -0.013643,  0.096166, ..., -0.060588, -0.059212,
       -0.091487], dtype=float32)
 y: array([ 0.038159, -0.013643,  0.096166, ..., -0.060588, -0.059212,
       -0.091487], dtype=float32)

======================================================================
FAIL: test_field_coarsening_via_subcycles_5 (__main__.AdvanceTest)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/opt/py/py/lib/python3.9/site-packages/ddt.py", line 182, in wrapper
    return func(self, *args, **kwargs)
  File "/home/philix/git/phare/master/build/tests/simulator/advance/2d/fields.py", line 58, in test_field_coarsening_via_subcycles
    self._test_field_coarsening_via_subcycles(ndim, interp_order, refinement_boxes)
  File "/home/philix/git/phare/master/build/tests/simulator/test_advance.py", line 496, in _test_field_coarsening_via_subcycles
    np.testing.assert_allclose(coarse_pdDataset, afterCoarse, atol=1e-6)
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 1528, in assert_allclose
    assert_array_compare(compare, actual, desired, err_msg=str(err_msg),
  File "/opt/py/py/lib/python3.9/site-packages/numpy/testing/_private/utils.py", line 842, in assert_array_compare
    raise AssertionError(msg)
AssertionError:
Not equal to tolerance rtol=1e-07, atol=1e-06

Mismatched elements: 5 / 702 (0.712%)
Max absolute difference: 0.12036321
Max relative difference: 5.8217904e-07
 x: array([-2.742830e-02, -3.085866e-02,  2.238256e-02,  7.562380e-02,
       -2.041059e-02, -1.164450e-01, -7.100931e-02, -2.557364e-02,
       -8.084903e-02, -1.361244e-01, -1.095806e-01, -8.303684e-02,...
 y: array([-2.742830e-02, -3.085866e-02,  2.238256e-02,  7.562380e-02,
       -2.041059e-02, -1.164450e-01, -7.100931e-02, -2.557364e-02,
       -8.084903e-02, -1.361244e-01,  0.000000e+00,  0.000000e+00,...

----------------------------------------------------------------------
Ran 10 tests in 302.621s
"""