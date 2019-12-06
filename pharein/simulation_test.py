import unittest
import numpy as np
import configparser

import simulation
import uniform_model
import phare_utilities
import diagnostics


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

    def test_dl(self):
        for cells, domain_size, dim, bc in zip(self.cells_array,
                                               self.domain_size_array,
                                               self.dims,
                                               self.bcs):

            j = simulation.Simulation(time_step_nbr=self.time_step_nbr,
                                      boundary_types=bc, cells=cells ,
                                      domain_size=domain_size, final_time=self.final_time)

            if phare_utilities.none_iterable(domain_size, cells):
                domain_size = phare_utilities.listify(domain_size)
                cells = phare_utilities.listify(cells)

            for d in np.arange(dim):
                self.assertEqual(j.dl[d], domain_size[d]/float(cells[d]))

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

    def test_ini_file(self):
        """ serialize a simulation in an INI file and read it to check it is OK"""

        simu = simulation.Simulation(time_step_nbr=1000,
                                     boundary_types="periodic",
                                     cells=80,
                                     dl=0.1,
                                     final_time=1.)

        fd1 = diagnostics.FluidDiagnostics(name="FluidDiagnostics1",
                                           diag_type="rho_s",
                                           write_every=10,
                                           compute_every=5,
                                           start_iteration=0,
                                           last_iteration=1000,
                                           species_name="proton1")

        fd2 = diagnostics.FluidDiagnostics(name="FluidDiagnostics2",
                                           diag_type="flux_s",
                                           write_every=10,
                                           compute_every=5,
                                           start_iteration=0,
                                           last_iteration=1000,
                                           species_name="proton1")

        simu.add_diagnostics(fd1)
        simu.add_diagnostics(fd2)

        ed1 = diagnostics.ElectromagDiagnostics(name="ElectromagDiagnostics1",
                                                diag_type="E",
                                                write_every=10,
                                                compute_evert=5,
                                                start_teration=0,
                                                last_iteration=1000)

        ed2 = diagnostics.ElectromagDiagnostics(name="ElectromagDiagnostics2",
                                                diag_type="B",
                                                write_every=10,
                                                compute_evert=5,
                                                start_teration=0,
                                                last_iteration=1000)

        simu.add_diagnostics(ed1)
        simu.add_diagnostics(ed2)

        model = uniform_model.UniformModel()
        model.add_fields()
        model.add_species("proton1")
        model.add_species("proton2", density=2.)

        simu.set_model(model)

        simu.write_ini_file()

        config=configparser.ConfigParser()
        config.read("phare.ini")
        sections = config.sections()
        expected_sections=('Simulation','model','FluidDiagnostics1',
                           'FluidDiagnostics2',
                           'ElectromagDiagnostics1',
                           'ElectromagDiagnostics2')

        for section in sections:
            self.assertIn(section, expected_sections)

        numOfCellx = config.getint('Simulation', 'nbr_cells_x')
        numOfCelly = config.getint('Simulation', 'nbr_cells_y', fallback=0)
        numOfCellz = config.getint('Simulation', 'nbr_cells_z', fallback=0)

        cells = -12
        if numOfCellz != 0 and numOfCelly != 0 and numOfCellx != 0:
            cells = [numOfCellx, numOfCelly, numOfCellz]

        elif numOfCellx != 0 and numOfCelly != 0 and numOfCellz == 0:
            cells = [numOfCellx, numOfCelly]

        elif numOfCellx != 0 and numOfCelly == 0 and numOfCellz == 0:
            cells = numOfCellx

        self.assertEqual(cells, 80)

        dx = config.getfloat('Simulation', 'dx')
        dy = config.getfloat('Simulation', 'dy', fallback=0)
        dz = config.getfloat('Simulation', 'dz', fallback=0)

        dl = -12
        if dz != 0 and dy != 0 and dx != 0:
            dl = [dx, dy, dz]

        elif dx != 0 and dy != 0 and dz == 0:
            dl = [dx, dy]

        elif dx != 0 and dy == 0 and dz == 0:
            dl = dx

        self.assertEqual(dl, 0.1)

        Ox = config.getfloat('Simulation', 'origin_x', fallback=-12)
        Oy = config.getfloat('Simulation', 'origin_y', fallback=-12)
        Oz = config.getfloat('Simulation', 'origin_z', fallback=-12)

        origin = Ox

        self.assertEqual(origin, 0.)

        bc = config.get('Simulation', 'boundary_condition_x')
        self.assertEqual('periodic', bc)

        pusher = config.get('Simulation', 'particle_pusher')
        self.assertEqual('modified_boris', pusher)

        refined_particle_nbr = config.get('Simulation', 'refined_particle_nbr')
        self.assertEqual(2, int(refined_particle_nbr))




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
