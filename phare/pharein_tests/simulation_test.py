import unittest
import numpy as np
import configparser
import pharein.globals as globals

from pharein import phare_utilities, uniform_model, diagnostics, simulation


class TestSimulation(unittest.TestCase):

    def setUp(self):
        self.cells_array = [80, (80, 40), (80, 40, 12)]
        self.dl_array = [0.1, (0.1, 0.2), (0.1, 0.2, 0.3)]
        self.domain_size_array = [100., (100., 80.), (100., 80., 20.)]
        self.dims = [1, 2, 3]
        self.bcs = ["periodic", ("periodic", "periodic"), ("periodic", "periodic", "periodic")]
        self.layout = "yee"
        self.time_step = 0.001
        self.time_step_nbr = 1000
        self.final_time = 1.
        globals.sim = None

    def test_dl(self):
        for cells, domain_size, dim, bc in zip(self.cells_array,
                                               self.domain_size_array,
                                               self.dims,
                                               self.bcs):

            j = simulation.Simulation(time_step_nbr=self.time_step_nbr,
                                      boundary_types=bc, cells=cells,
                                      domain_size=domain_size, final_time=self.final_time)

            if phare_utilities.none_iterable(domain_size, cells):
                domain_size = phare_utilities.listify(domain_size)
                cells = phare_utilities.listify(cells)

            for d in np.arange(dim):
                self.assertEqual(j.dl[d], domain_size[d]/float(cells[d]))

            globals.sim = None

    def test_boundary_conditions(self):
        j = simulation.Simulation(time_step_nbr=1000, boundary_types="periodic",
                                  cells=80, domain_size=10, final_time=1.)

        for d in np.arange(j.dims):
            self.assertEqual("periodic", j.boundary_types[d])

    def test_assert_boundary_condition(self):
        simulation.Simulation(time_step_nbr=1000,
                              boundary_types="periodic",
                              cells=80, domain_size=10, final_time=1000)

    def test_time_step(self):
        s = simulation.Simulation(time_step_nbr=1000,
                                  boundary_types="periodic",
                                  cells=80, domain_size=10, final_time=10)
        self.assertEqual(0.01, s.time_step)




###############################################################################
#
#                   usage example and run tests
#
#       the script: - creates a Simulation
#                   - registers 2 fluid diagnostics to the Simulation
#                   - registers 2 electromag diagnostics to the Simulation
#                   - creates a model 'UniformModel'
#                   - adds default electromagnetic fields to the model
#                   - adds 2 proton species, the second has a density of 2.
#                   - registers the model to the Simulation
#                   - prepare the job directories and input file
#
#                   - executs unit tests of this module
###############################################################################
if __name__ == '__main__':



    unittest.main()
