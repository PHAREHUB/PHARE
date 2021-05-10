
import unittest
from concurrencytest import ConcurrentTestSuite, fork_for_tests

from tests.simulator.refinement_boxes import SimulatorRefineBoxInputs as SimulatorValidation

from tests.simulator.initialize._1d.fields import InitializationTest as InitField1d
from tests.simulator.initialize._1d.particles import InitializationTest as InitParticles1d

from tests.simulator.advance._1d.fields import AdvanceTest as AdvanceField1d
from tests.simulator.advance._1d.particles import AdvanceTest as AdvanceParticles1d

from tests.simulator.initialize._2d.fields import InitializationTest as InitField2d
from tests.simulator.initialize._2d.particles import InitializationTest as InitParticles2d

from tests.simulator.advance._2d.fields import AdvanceTest as AdvanceField2d
from tests.simulator.advance._2d.particles import AdvanceTest as AdvanceParticles2d

if __name__ == "__main__":

    test_classes_to_run = [
      # SimulatorValidation,
      InitField1d, InitParticles1d,
      AdvanceField1d, AdvanceParticles1d,
      InitField2d, InitParticles2d,
      AdvanceField2d, AdvanceParticles2d
    ]

    loader = unittest.TestLoader()

    suite = unittest.TestSuite()
    for test_class in test_classes_to_run:
        suite.addTest(loader.loadTestsFromTestCase(test_class))

    unittest.TextTestRunner().run(ConcurrentTestSuite(suite, fork_for_tests(30)))
