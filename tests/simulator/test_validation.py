#!/usr/bin/env python3
#
# formatted with black

from pyphare.cpp import cpp_lib

cpp = cpp_lib()

import unittest

from ddt import data, ddt
from pyphare.core.box import Box, Box2D
from pyphare.simulator.simulator import Simulator

from tests.simulator import NoOverwriteDict, populate_simulation
from tests.simulator import SimulatorTest

out = "phare_outputs/valid/refinement_boxes/"
diags = {
    "diag_options": {"format": "phareh5", "options": {"dir": out, "mode": "overwrite"}}
}
restarts = {"restart_options": {"dir": out, "mode": "overwrite"}}


def dup(dic):
    dic = NoOverwriteDict(dic)
    dic.update(diags.copy())
    dic.update(restarts.copy())
    return dic


@ddt
class SimulatorValidation(SimulatorTest):
    def __init__(self, *args, **kwargs):
        super(SimulatorValidation, self).__init__(*args, **kwargs)
        self.simulator = None

    def tearDown(self):
        super(SimulatorValidation, self).tearDown()
        if self.simulator is not None:
            self.simulator.reset()

    def _do_dim(self, dim, input, valid: bool = False):
        for interp in range(1, 4):
            try:
                self.simulator = Simulator(populate_simulation(dim, interp, **input))
                self.simulator.setup()
                self.assertTrue(valid)
                self.simulator = None
            except ValueError as e:
                self.assertTrue(not valid)

    """
      The first set of boxes "B0": [(10,), (14,)]
      Are configured to force there to be a single patch on L0
      This creates a case with MPI that there are an unequal number of
      Patches across MPI domains. This case must be handled and not hang due
      to collective calls not being handled properly.
    """
    valid1D = [
        dup({"cells": [65], "refinement_boxes": {"L0": {"B0": [(10,), (14,)]}}}),
        dup({"cells": [65], "refinement_boxes": {"L0": {"B0": [(5,), (55,)]}}}),
        dup({"cells": [65], "refinement_boxes": {"L0": {"B0": Box(5, 55)}}}),
        dup({"cells": [65], "refinement_boxes": {"L0": [Box(5, 55)]}}),
        dup({"cells": [65], "refinement_boxes": {0: [Box(5, 55)]}}),
        dup({"cells": [65], "refinement_boxes": {0: [Box(0, 55)]}}),
        dup({"cells": [65], "refinement_boxes": {"L0": [Box(5, 14), Box(15, 25)]}}),
        dup(
            {
                "cells": [65],
                "refinement_boxes": {
                    "L0": [Box(5, 25)],
                    "L1": [Box(12, 48)],
                    "L2": [Box(60, 64)],
                },
            }
        ),
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(12, 48)]},
            }
        ),
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(20, 30)]},
            }
        ),
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(11, 49)]},
                "nesting_buffer": 1,
            }
        ),
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(10, 50)]},
            }
        ),
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(15, 49)]},
            }
        ),
        dup(
            {
                "cells": [65],
                "refinement_boxes": None,
                "smallest_patch_size": 20,
                "largest_patch_size": 20,
                "nesting_buffer": 10,
            }
        ),
        # finer box is within set of coarser boxes
        dup(
            {
                "cells": [65],
                "refinement_boxes": {
                    "L0": [Box(5, 9), Box(10, 15)],
                    "L1": [Box(11, 29)],
                },
            }
        ),
    ]

    invalid1D = [
        # finer box outside lower
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 24)], "L1": [Box(9, 30)]},
            }
        ),
        # finer box outside upper
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 24)], "L1": [Box(15, 50)]},
            }
        ),
        # overlapping boxes
        dup({"cells": [65], "refinement_boxes": {"L0": [Box(5, 15), Box(15, 25)]}}),
        # box.upper outside domain
        dup({"cells": [55], "refinement_boxes": {"L0": {"B0": [(5,), (65,)]}}}),
        # largest_patch_size > smallest_patch_size
        dup(
            {
                "smallest_patch_size": 100,
                "largest_patch_size": 64,
            }
        ),
        # refined_particle_nbr doesn't exist
        dup({"refined_particle_nbr": 1}),
        # L2 box incompatible with L1 box due to nesting buffer
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(11, 49)]},
                "nesting_buffer": 2,
            }
        ),
        # negative nesting buffer
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(11, 49)]},
                "nesting_buffer": -1,
            }
        ),
        # too large nesting buffer
        dup(
            {
                "cells": [65],
                "refinement_boxes": {"L0": [Box(5, 25)], "L1": [Box(11, 49)]},
                "nesting_buffer": 33,
            }
        ),
        dup(
            {
                "cells": [65],
                "refinement_boxes": None,
                "largest_patch_size": 20,
                "nesting_buffer": 46,
            }
        ),
        # finer box is not within set of coarser boxes
        dup(
            {
                "cells": [65],
                "refinement_boxes": {
                    "L0": [Box(5, 9), Box(11, 15)],
                    "L1": [Box(11, 29)],
                },
            }
        ),
    ]

    @data(*valid1D)
    def test_1d_valid(self, input):
        self._do_dim(1, input, True)

    @data(*invalid1D)
    def test_1d_invalid(self, input):
        self._do_dim(1, input)

    valid2D = [
        dup({"cells": [65, 65], "refinement_boxes": {"L0": [Box2D(5, 55)]}}),
        dup({"smallest_patch_size": None, "largest_patch_size": None}),
        dup({"smallest_patch_size": (10, 10), "largest_patch_size": (20, 20)}),
        dup({"smallest_patch_size": [10, 10], "largest_patch_size": [20, 20]}),
        dup({"smallest_patch_size": (10, 10), "largest_patch_size": None}),
        dup({"smallest_patch_size": None, "largest_patch_size": (20, 20)}),
        dup({"smallest_patch_size": (10, 10)}),
        dup({"largest_patch_size": (20, 20)}),
        dup({"smallest_patch_size": 10, "largest_patch_size": (20, 20)}),
        dup({"smallest_patch_size": (10, 10), "largest_patch_size": 20}),
        dup({"smallest_patch_size": [10, 10], "largest_patch_size": (20, 20)}),
        dup({"smallest_patch_size": (10, 10), "largest_patch_size": [20, 20]}),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": None,
                "smallest_patch_size": 20,
                "largest_patch_size": 20,
                "nesting_buffer": 10,
            }
        ),
        dup({"cells": [65, 65], "refinement_boxes": {"L0": {"B0": Box2D(5, 55)}}}),
        dup({"cells": [65, 65], "refinement_boxes": {"L0": [Box2D(5, 55)]}}),
        dup({"cells": [65, 65], "refinement_boxes": {0: [Box2D(5, 55)]}}),
        dup({"cells": [65, 65], "refinement_boxes": {0: [Box2D(0, 55)]}}),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 14), Box2D(15, 25)]},
            }
        ),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {
                    "L0": [Box2D(5, 25)],
                    "L1": [Box2D(12, 48)],
                    "L2": [Box2D(60, 64)],
                },
            }
        ),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(12, 48)]},
            }
        ),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(20, 30)]},
            }
        ),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(11, 49)]},
                "nesting_buffer": 1,
            }
        ),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(10, 50)]},
            }
        ),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(15, 49)]},
            }
        ),
    ]

    invalid2D = [
        # finer box outside lower
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 24)], "L1": [Box2D(9, 30)]},
            }
        ),
        # finer box outside lower
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 24)], "L1": [Box2D(9, 30)]},
            }
        ),
        # finer box outside upper
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 24)], "L1": [Box2D(15, 50)]},
            }
        ),
        # overlapping boxes
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 15), Box2D(15, 25)]},
            }
        ),
        # box.upper outside domain
        dup(
            {
                "cells": [55, 55],
                "refinement_boxes": {
                    "L0": {
                        "B0": Box2D(
                            5,
                            65,
                        )
                    }
                },
            }
        ),
        # largest_patch_size > smallest_patch_size
        dup(
            {
                "smallest_patch_size": 100,
                "largest_patch_size": 64,
            }
        ),
        # refined_particle_nbr doesn't exist
        dup({"refined_particle_nbr": 1}),
        # L2 box incompatible with L1 box due to nesting buffer
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(11, 49)]},
                "nesting_buffer": 2,
            }
        ),
        # negative nesting buffer
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(11, 49)]},
                "nesting_buffer": -1,
            }
        ),
        # too large nesting buffer
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {"L0": [Box2D(5, 25)], "L1": [Box2D(11, 49)]},
                "nesting_buffer": 33,
            }
        ),
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": None,
                "largest_patch_size": 20,
                "nesting_buffer": 46,
            }
        ),
        # finer box is not within set of coarser boxes
        dup(
            {
                "cells": [65, 65],
                "refinement_boxes": {
                    "L0": [Box2D(5, 9), Box2D(11, 15)],
                    "L1": [Box2D(11, 29)],
                },
            }
        ),
    ]

    @data(*valid2D)
    def test_2d_valid(self, input):
        self._do_dim(2, input, True)

    @data(*invalid2D)
    def test_2d_invalid(self, input):
        self._do_dim(2, input)


if __name__ == "__main__":
    unittest.main()
