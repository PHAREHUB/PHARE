"""
finest data of a deposited moment must not dip at coarse/fine borders, see
https://github.com/PHAREHUB/PHARE/issues/1558

At a coarse/fine border, level ghost particles only fill the first ghost cell(s), so the
moment ghost nodes there are only partially deposited. At interp order 1 the kept ghost node
was half filled, plotting a dip to about half the density on a uniform density.
"""

import unittest

import numpy as np

from pyphare.pharesee.run import Run
from pyphare.pharesee.run.utils import flat_finest_field
from pyphare.simulator.simulator import Simulator, startMPI

diag_dir = "phare_outputs/test_finest_moments_1d"
density = 2.2
time_step = 0.005
final_time = 10 * time_step
timestamps = [0.0, final_time]


def config():
    import pyphare.pharein as ph

    sim = ph.Simulation(
        smallest_patch_size=20,
        largest_patch_size=60,
        time_step=time_step,
        final_time=final_time,
        cells=500,
        dl=0.2,
        nesting_buffer=1,
        interp_order=1,
        refinement_boxes={"L0": {"B0": [(80,), (120,)], "B1": [(380,), (420,)]}},
        diag_options={
            "format": "phareh5",
            "options": {"dir": diag_dir, "mode": "overwrite"},
        },
    )

    def zero(x):
        return np.zeros_like(x)

    def vth(x):
        return 0.7 * np.ones_like(x)

    ph.MaxwellianFluidModel(
        bx=zero,
        by=zero,
        bz=zero,
        protons={
            "charge": 1,
            "density": lambda x: density * np.ones_like(x),
            "vbulkx": zero,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": 100,
            "init": {"seed": 1337},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.5)
    ph.FluidDiagnostics(quantity="charge_density", write_timestamps=timestamps)
    return sim


class FinestMoments1dTest(unittest.TestCase):
    def test_finest_charge_density_is_flat(self):
        import pyphare.pharein as ph

        ph.global_vars.sim = None
        Simulator(config()).run().reset()
        ph.global_vars.sim = None

        run = Run(diag_dir)
        for time in timestamps:
            rho = run.GetNi(time)
            self.assertEqual(rho.finest_level(time), 1)

            data, x = flat_finest_field(rho, "value", time=time)
            # 100 particles per cell: about 4% noise per node, the dips were about 45%
            deviation = np.abs(data / density - 1)
            worst = np.argmax(deviation)
            self.assertLess(
                deviation[worst],
                0.25,
                f"t={time}: {data[worst]} at x={x[worst]}, expected about {density}",
            )


if __name__ == "__main__":
    startMPI()
    unittest.main()
