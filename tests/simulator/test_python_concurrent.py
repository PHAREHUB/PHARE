"""
  This script exists to minimize testing time by running all simulation/tests
    concurrently without needing to wait for any particular file or set of tests
"""

import os
import time
import unittest
from multiprocessing import cpu_count

from tests.simulator.test_validation import SimulatorValidation
from tests.simulator.test_diagnostics import DiagnosticsTest  # mpirun -n 1/2/3/4

import tests.simulator.initialize.test_fields_init_1d as test_fields_init_1d
# from tests.simulator.initialize.test_fields_init_1d import (
#     InitializationTest as InitField1d,
# )
# from tests.simulator.initialize.test_particles_init_1d import (
#     InitializationTest as InitParticles1d,
# )
# from tests.simulator.advance.test_fields_advance_1d import AdvanceTest as AdvanceField1d
# from tests.simulator.advance.test_particles_advance_1d import (
#     AdvanceTest as AdvanceParticles1d,
# )
# from tests.simulator.initialize.test_fields_init_2d import (
#     InitializationTest as InitField2d,
# )
# from tests.simulator.initialize.test_particles_init_2d import (
#     InitializationTest as InitParticles2d,
# )
# from tests.simulator.advance.test_fields_advance_2d import AdvanceTest as AdvanceField2d
# from tests.simulator.advance.test_particles_advance_2d import (
#     AdvanceTest as AdvanceParticles2d,
# )

from phlop.testing.parallel_processor import load_test_cases_in, process

N_CORES = int(os.environ["N_CORES"]) if "N_CORES" in os.environ else cpu_count()
PRINT = int(os.environ["PRINT"]) if "PRINT" in os.environ else 0
MPI_RUN = os.environ.get("MPI_RUN", None)
MPI_RUN_EXTRA = os.environ.get("MPI_RUN_EXTRA", "")

def test_cmd(clazz, test_id, mpi_run):
    return f"mpirun {MPI_RUN_EXTRA} -n {mpi_run} python3 -m {clazz.__module__} {clazz.__name__}.{test_id}"



test_classes_to_run = [
    SimulatorValidation,
    DiagnosticsTest,
    test_fields_init_1d.InitializationTest,
    # InitField1d,
    # InitParticles1d,
    # AdvanceField1d,
    # AdvanceParticles1d,
    # InitField2d,
    # InitParticles2d,
    # AdvanceField2d,
    # AdvanceParticles2d,
]


def load_test_cases(classes, cores):
    return load_test_cases_in(classes, cores, test_cmd_fn=test_cmd)

def build_batches():
    batches = []
    if not MPI_RUN:
        batches += [load_test_cases(test_classes_to_run, 1)]
        batches += [load_test_cases(test_classes_to_run, 2)]
        batches += [load_test_cases([DiagnosticsTest], 3)]
        batches += [load_test_cases([DiagnosticsTest], 4)]
    else:
        batches += [load_test_cases(test_classes_to_run, int(MPI_RUN))]
    return batches


def main():
    process(build_batches())


if __name__ == "__main__":
    main()
