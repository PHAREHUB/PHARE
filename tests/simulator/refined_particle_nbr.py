#!/usr/bin/env python3
#
# formatted with black

from pyphare.cpp import cpp_lib
cpp = cpp_lib()

import os, sys, unittest, yaml
import numpy as np
import pyphare.pharein as ph
from pyphare.simulator.simulator import Simulator
from tests.simulator import NoOverwriteDict
from tests.simulator import populate_simulation
from pyphare.cpp import splitter_type

from tests.simulator.config import project_root


class SimulatorRefinedParticleNbr(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(SimulatorRefinedParticleNbr, self).__init__(*args, **kwargs)
        with open(os.path.join(project_root, "res/amr/splitting.yml"), 'r') as stream:
            try:
                self.yaml_root = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
                sys.exit(1)

    # while splitting particles may leave the domain area
    #  so we remove the particles from the border cells of each patch
    def _less_per_dim(self, dim, refined_particle_nbr, patch):
        if dim == 1:
            return refined_particle_nbr * 2
        cellNbr = patch.upper - patch.lower + 1
        if dim == 2:
            return refined_particle_nbr * ((cellNbr[0] * 2 + (cellNbr[1] * 2)))
        raise ValueError("Unhandled dimension for function")


    def _check_deltas_and_weights(self, dim, interp, refined_particle_nbr):
        yaml_dim = self.yaml_root["dimension_" + str(dim)]
        yaml_interp = yaml_dim["interp_" + str(interp)]
        yaml_n_particles = yaml_interp["N_particles_" + str(refined_particle_nbr)]

        yaml_delta = [float(s) for s in str(yaml_n_particles["delta"]).split(" ")]
        yaml_weight = [float(s) for s in str(yaml_n_particles["weight"]).split(" ")]

        splitter_t = splitter_type(dim, interp, refined_particle_nbr)
        np.testing.assert_allclose(yaml_delta, splitter_t.delta)
        np.testing.assert_allclose(yaml_weight, splitter_t.weight)


    def _do_dim(self, dim, min_diff, max_diff):
        from pyphare.pharein.simulation import valid_refined_particle_nbr

        for interp in range(1, 4):

            prev_split_particle_max = 0
            for refined_particle_nbr in valid_refined_particle_nbr[dim][interp]:

                self._check_deltas_and_weights(dim, interp, refined_particle_nbr)

                simInput = NoOverwriteDict({"refined_particle_nbr": refined_particle_nbr})
                self.simulator = Simulator(populate_simulation(dim, interp, **simInput))
                self.simulator.initialize()
                dw = self.simulator.data_wrangler()
                max_per_pop = 0
                leaving_particles = 0
                for pop, particles in dw.getPatchLevel(1).getParticles().items():
                    per_pop = 0
                    for key, patches in particles.items():
                        for patch in patches:
                            leaving_particles += self._less_per_dim(dim, refined_particle_nbr, patch)
                            per_pop += patch.data.size()
                    max_per_pop = max(max_per_pop, per_pop)
                prev_min_diff = prev_split_particle_max * min_diff

                # while splitting particles may leave the domain area
                #  so we remove the particles from the border cells of each patch
                self.assertTrue(
                    max_per_pop > prev_min_diff - (leaving_particles)
                )
                if prev_split_particle_max > 0:
                    prev_max_diff = prev_min_diff * dim * max_diff
                    self.assertTrue(max_per_pop < prev_max_diff)
                prev_split_particle_max = max_per_pop
                self.simulator = None

    """ 1d
      refine 10 cells in 1d, ppc 100
        10 * 2 * ppc = 200
        10 * 3 * ppc = 300 300 / 200 = 1.5
        10 * 4 * ppc = 400 500 / 400 = 1.33
        10 * 5 * ppc = 500 500 / 400 = 1.25
      taking the minimul diff across permutations
      current to previous should be at least this times bigger
    """
    PREVIOUS_ITERATION_MIN_DIFF_1d = 1.25
    PREVIOUS_ITERATION_MAX_DIFF_1d = 1.50

    def test_1d(self):
        This = type(self)
        self._do_dim(1, This.PREVIOUS_ITERATION_MIN_DIFF_1d, This.PREVIOUS_ITERATION_MAX_DIFF_1d)

    """ 2d
      refine 10x10 cells in 2d, ppc 100
      10 * 10 * 4 * ppc = 400
      10 * 10 * 8 * ppc = 800 800 / 400 = 1.5
      10 * 10 * 9 * ppc = 900 900 / 800 = 1.125
    """
    PREVIOUS_ITERATION_MIN_DIFF_2d = 1.125
    PREVIOUS_ITERATION_MAX_DIFF_2d = 1.50

    def test_2d(self):
        This = type(self)
        self._do_dim(2, This.PREVIOUS_ITERATION_MIN_DIFF_2d, This.PREVIOUS_ITERATION_MAX_DIFF_2d)

    def tearDown(self):
        # needed in case exception is raised in test and Simulator
        # not reset properly
        if self.simulator is not None:
            self.simulator.reset()



if __name__ == "__main__":
    unittest.main()
