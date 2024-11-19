#!/usr/bin/env python3
#
# formatted with black


import os
import sys
import yaml
import unittest
import numpy as np

from pyphare.cpp import cpp_lib
from pyphare.cpp import splitter_type
from pyphare.simulator.simulator import Simulator

from tests.simulator import NoOverwriteDict, populate_simulation
from tests.simulator.config import project_root

cpp = cpp_lib()


class SimulatorRefinedParticleNbr(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(SimulatorRefinedParticleNbr, self).__init__(*args, **kwargs)
        self.simulator = None
        with open(os.path.join(project_root, "res/amr/splitting.yml"), "r") as stream:
            try:
                self.yaml_root = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
                sys.exit(1)

    # needed in case exception is raised in test and Simulator not reset properly
    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()

    # while splitting particles may leave the domain area
    #  so we remove the particles from the border cells of each patch
    def _less_per_dim(self, dim, refined_particle_nbr, patch):
        if dim == 1:
            return refined_particle_nbr * 2
        cellNbr = patch.upper - patch.lower + 1
        dim2 = refined_particle_nbr * ((cellNbr[0] * 2 + (cellNbr[1] * 2)))
        if dim == 2:
            return dim2
        # TODO3D if dim==3?
        return dim2 * (cellNbr[2] * 2)

    def _do_dim(self, dim, min_diff, max_diff):
        from pyphare.pharein.simulation import valid_refined_particle_nbr

        for interp in range(1, 4):
            prev_split_particle_max = 0
            for refined_particle_nbr in valid_refined_particle_nbr[dim][interp]:
                simInput = {"refined_particle_nbr": refined_particle_nbr}
                self.simulator = Simulator(populate_simulation(dim, interp, **simInput))
                self.simulator.initialize()
                dw = self.simulator.data_wrangler()
                max_per_pop = 0
                leaving_particles = 0
                for pop, particles in dw.getPatchLevel(1).getParticles().items():
                    per_pop = 0
                    for key, patches in particles.items():
                        for patch in patches:
                            leaving_particles += self._less_per_dim(
                                dim, refined_particle_nbr, patch
                            )
                            per_pop += patch.data.size()
                    max_per_pop = max(max_per_pop, per_pop)
                prev_min_diff = prev_split_particle_max * min_diff
                self.assertTrue(max_per_pop > prev_min_diff - (leaving_particles))
                if prev_split_particle_max > 0:
                    prev_max_diff = prev_min_diff * dim * max_diff
                    self.assertTrue(max_per_pop < prev_max_diff)
                prev_split_particle_max = max_per_pop
                self.simulator = None

    #  1d
    #  refine 10 cells in 1d, ppc 100
    #    10 * 2 * ppc = 200
    #    10 * 3 * ppc = 300 300 / 200 = 1.5
    #    10 * 4 * ppc = 400 500 / 400 = 1.33
    #    10 * 5 * ppc = 500 500 / 400 = 1.25
    #  taking the minimul diff across permutations
    #  current to previous should be at least this times bigger
    PRIOR_MIN_DIFF_1d = 1.25
    PRIOR_MAX_DIFF_1d = 1.50

    def test_1d(self):
        This = type(self)
        self._do_dim(1, This.PRIOR_MIN_DIFF_1d, This.PRIOR_MAX_DIFF_1d)

    #  2d
    #  refine 10x10 cells in 2d, ppc 100
    #    10 * 10 * 4 * ppc = 400
    #    10 * 10 * 8 * ppc = 800 800 / 400 = 1.5
    #    10 * 10 * 9 * ppc = 900 900 / 800 = 1.125
    PRIOR_MIN_DIFF_2d = 1.125
    PRIOR_MAX_DIFF_2d = 1.50

    def test_2d(self):
        This = type(self)
        self._do_dim(2, This.PRIOR_MIN_DIFF_2d, This.PRIOR_MAX_DIFF_2d)

    # 3d
    #  refine 10x10x10 cells in 3d, ppc 100
    #    10 * 10 * 10 * 6 *  ppc = 6000
    #    10 * 10 * 10 * 12 * ppc = 12000 - 12000 / 6000  = 2
    #    10 * 10 * 10 * 27 * ppc = 27000 - 27000 / 12000 = 2.25
    PRIOR_MIN_DIFF_3d = 2
    PRIOR_MAX_DIFF_3d = 2.25

    def test_3d(self):
        This = type(self)
        self._do_dim(3, This.PRIOR_MIN_DIFF_3d, This.PRIOR_MAX_DIFF_3d)

    def _check_deltas_and_weights(self, dim, interp, refined_particle_nbr):
        yaml_dim = self.yaml_root["dimension_" + str(dim)]
        yaml_interp = yaml_dim["interp_" + str(interp)]
        yaml_n_particles = yaml_interp["N_particles_" + str(refined_particle_nbr)]

        yaml_delta = [float(s) for s in str(yaml_n_particles["delta"]).split(" ")]
        yaml_weight = [float(s) for s in str(yaml_n_particles["weight"]).split(" ")]

        splitter_t = splitter_type(dim, interp, refined_particle_nbr)
        np.testing.assert_allclose(yaml_delta, splitter_t.delta)
        np.testing.assert_allclose(yaml_weight, splitter_t.weight)

    def test_values(self):
        from pyphare.pharein.simulation import valid_refined_particle_nbr

        for dim in range(1, 4):
            for interp in range(1, 4):
                for refined_particle_nbr in valid_refined_particle_nbr[dim][interp]:
                    self._check_deltas_and_weights(dim, interp, refined_particle_nbr)


if __name__ == "__main__":
    # the following ensures the order of execution so delta/weights are verified first
    suite = unittest.TestSuite()
    suite.addTest(SimulatorRefinedParticleNbr("test_values"))
    for dim in range(1, 4):
        suite.addTest(SimulatorRefinedParticleNbr("test_" + str(dim) + "d"))
    unittest.TextTestRunner(failfast=True).run(suite)
