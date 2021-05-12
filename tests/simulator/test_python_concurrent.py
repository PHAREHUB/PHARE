
import unittest
from concurrencytest import ConcurrentTestSuite, fork_for_tests

from tests.simulator.refinement_boxes import SimulatorRefineBoxInputs as SimulatorValidation

from tests.simulator.initialize.test_fields_init_1d import InitializationTest as InitField1d
from tests.simulator.initialize.test_particles_init_1d import InitializationTest as InitParticles1d

from tests.simulator.advance.test_fields_advance_1d import AdvanceTest as AdvanceField1d
from tests.simulator.advance.test_particles_advance_1d import AdvanceTest as AdvanceParticles1d

from tests.simulator.initialize.test_fields_init_2d import InitializationTest as InitField2d
from tests.simulator.initialize.test_particles_init_2d import InitializationTest as InitParticles2d

from tests.simulator.advance.test_fields_advance_2d import AdvanceTest as AdvanceField2d
from tests.simulator.advance.test_particles_advance_2d import AdvanceTest as AdvanceParticles2d

if __name__ == "__main__":

    test_classes_to_run = [
      SimulatorValidation,
      InitField1d,
      InitParticles1d,
      AdvanceField1d,
      AdvanceParticles1d,
      InitField2d,
      InitParticles2d,
      AdvanceField2d,
      AdvanceParticles2d
    ]

    loader = unittest.TestLoader()

    suite = unittest.TestSuite()
    for test_class in test_classes_to_run:
        suite.addTest(loader.loadTestsFromTestCase(test_class))

    import multiprocessing

    unittest.TextTestRunner().run(ConcurrentTestSuite(suite, fork_for_tests(multiprocessing.cpu_count())))


"""
import os
import unittest

ignore_paths = ["/functional/", "test_diagnostic_timestamps.py"]
tests_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

class CustomTestLoader(unittest.TestLoader):
    def _find_test_path(self, full_path, pattern="*.py", namespace=False):
        for ignore_path in ignore_paths:
            path = f"{ignore_path}"
            if path in full_path:
                return (None, False)

        original_isfile = os.path.isfile

        def patched_isfile(path):
            return str(path).endswith('__init__.py') or original_isfile(path)

        os.path.isfile = patched_isfile
        result = super()._find_test_path(full_path=full_path, pattern=pattern,
                                         namespace=namespace)
        os.path.isfile = original_isfile
        return result

if __name__ == "__main__":
    from concurrencytest import ConcurrentTestSuite, fork_for_tests
    all_tests = CustomTestLoader().discover(tests_dir)
    unittest.TextTestRunner().run(ConcurrentTestSuite(all_tests, fork_for_tests(22)))
"""