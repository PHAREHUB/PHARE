"""
  This script exists to minimize testing time by running all simulation/tests
    concurrently without needing to wait for any particular file or set of tests
"""

import os
import time
import unittest
from multiprocessing import Process, Queue, cpu_count

from tests.simulator.test_validation import SimulatorValidation
from tests.simulator.test_diagnostics import DiagnosticsTest  # mpirun -n 1/2/3/4
from tests.simulator.initialize.test_fields_init_1d import InitializationTest as InitField1d
from tests.simulator.initialize.test_particles_init_1d import InitializationTest as InitParticles1d
from tests.simulator.advance.test_fields_advance_1d import AdvanceTest as AdvanceField1d
from tests.simulator.advance.test_particles_advance_1d import AdvanceTest as AdvanceParticles1d
from tests.simulator.initialize.test_fields_init_2d import InitializationTest as InitField2d
from tests.simulator.initialize.test_particles_init_2d import InitializationTest as InitParticles2d
from tests.simulator.advance.test_fields_advance_2d import AdvanceTest as AdvanceField2d
from tests.simulator.advance.test_particles_advance_2d import AdvanceTest as AdvanceParticles2d

from tools.python3 import run, find_on_path

N_CORES = int(os.environ["N_CORES"]) if "N_CORES" in os.environ else cpu_count()
MPI_RUN = os.environ["MPI_RUN"] if "MPI_RUN" in os.environ else 1
PRINT   = int(os.environ["PRINT"]) if "PRINT" in os.environ else 0

MPI_RUN_EXTRA = os.environ["MPI_RUN_EXTRA"] if "MPI_RUN_EXTRA" in os.environ else ""

def test_cmd(clazz, test_id, mpi_run):
    return f"mpirun {MPI_RUN_EXTRA} -n {mpi_run} python3 -m {clazz.__module__} {clazz.__name__}.{test_id}"

test_classes_to_run = [
  # SimulatorValidation, DiagnosticsTest,
  # InitField1d,         InitParticles1d,
  # AdvanceField1d,      AdvanceParticles1d,
  # InitField2d,         InitParticles2d,
  # AdvanceField2d,
  AdvanceParticles2d
]

class TestBatch:
    def __init__(self, tests, mpi_run = 1):
        self.tests = tests
        self.mpi_run = mpi_run

def load_test_cases_in(classes, mpi_run = 1):
    tests, loader = [], unittest.TestLoader()
    for test_class in classes:
        for suite in loader.loadTestsFromTestCase(test_class):
            tests += [test_cmd(type(suite), suite._testMethodName, mpi_run)]
    return TestBatch(tests, mpi_run)

def build_batches():
    batches = []
    if MPI_RUN=="cmake":
        batches += [load_test_cases_in(test_classes_to_run, 1)]
        batches += [load_test_cases_in(test_classes_to_run, 2)]
        batches += [load_test_cases_in([DiagnosticsTest], 3)]
        batches += [load_test_cases_in([DiagnosticsTest], 4)]
    else:
        batches += [load_test_cases_in(test_classes_to_run, int(MPI_RUN))]
    return batches

class CallableTest:
    def __init__(self, bi, ti, cmd):
        self.bi = bi
        self.ti = ti
        self.cmd = cmd
        self.run = None

    def __call__(self, **kwargs):
        self.run = run(self.cmd.split(), shell=False, capture_output=True, check=True, print_cmd=False)
        return self

class CoreCount:
  def __init__(self, cores_avail):
      self.cores_avail = cores_avail
      self.proces = []
      self.fin = []

def runner(runnable, queue):
    queue.put(runnable())

def print_tests(batches):
    for batch in batches:
        for test in batch.tests:
            print(test)

if __name__ == "__main__":
    batches = build_batches()
    # batches = [TestBatch([batches[0].tests[0]])]
    if PRINT:
        print_tests(batches)

    else:
        cc = CoreCount(N_CORES)
        assert cc.cores_avail >= max([batch.mpi_run for batch in batches])
        cc.procs = [[] for batch in batches]
        cc.fin = [0 for batch in batches]
        pqueue = Queue()

        def launch_tests():
            for bi, batch in enumerate(batches):
                offset = len(cc.procs[bi])
                for ti, test in enumerate(batch.tests[offset:]):
                    ti += offset
                    if batch.mpi_run <= cc.cores_avail:
                        test = CallableTest(bi, ti, batches[bi].tests[ti])
                        cc.cores_avail -= batch.mpi_run
                        cc.procs[bi] += [Process(target=runner, args=(test, (pqueue),))]
                        cc.procs[bi][-1].daemon = True
                        cc.procs[bi][-1].start()

        def finished():
            b = True
            for bi, batch in enumerate(batches):
                b &= cc.fin[bi] == len(batch.tests)
            return b

        def waiter(queue):
            while True:
                proc = queue.get()
                time.sleep(.01) # don't throttle!
                if (isinstance(proc, CallableTest)):
                    print(proc.cmd, f"finished in {proc.run.t:.2f} seconds")
                    cc.cores_avail += batches[proc.bi].mpi_run
                    cc.fin[proc.bi] += 1
                    launch_tests()
                    if finished(): break

        launch_tests()
        waiter(pqueue)

