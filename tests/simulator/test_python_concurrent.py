"""
  This script exists to minimize testing time by running all simulation/tests
    concurrently without needing to wait for any particular file or set of tests
"""

import os
import unittest
import multiprocessing

from tests.simulator.test_validation import SimulatorValidation

from tests.simulator.initialize.test_fields_init_1d import InitializationTest as InitField1d
from tests.simulator.initialize.test_particles_init_1d import InitializationTest as InitParticles1d

from tests.simulator.advance.test_fields_advance_1d import AdvanceTest as AdvanceField1d
from tests.simulator.advance.test_particles_advance_1d import AdvanceTest as AdvanceParticles1d

from tests.simulator.initialize.test_fields_init_2d import InitializationTest as InitField2d
from tests.simulator.initialize.test_particles_init_2d import InitializationTest as InitParticles2d

from tests.simulator.advance.test_fields_advance_2d import AdvanceTest as AdvanceField2d
from tests.simulator.advance.test_particles_advance_2d import AdvanceTest as AdvanceParticles2d


N_CORES = int(os.environ["N_CORES"]) if "N_CORES" in os.environ else multiprocessing.cpu_count()
MPI_RUN = int(os.environ["MPI_RUN"]) if "MPI_RUN" in os.environ else 1

def test_cmd(clazz, test_id):
    return f"mpirun -n {MPI_RUN} python3 -m {clazz.__module__} {clazz.__name__}.{test_id}"

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

    tests = []
    loader = unittest.TestLoader()
    for test_class in test_classes_to_run:
        for suite in loader.loadTestsFromTestCase(test_class):
            tests += [test_cmd(type(suite), suite._testMethodName)]

    from tools.python3 import run_mp
    run_mp(tests, N_CORES, check=True)
