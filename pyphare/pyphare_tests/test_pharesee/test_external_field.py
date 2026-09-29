#!/usr/bin/env python3
"""
Reading the external magnetic field B0 back out, in memory and from a dump.

These are post-processing tests: that B0 is *computed* at the right moments is
the simulator's business and is covered by tests/simulator/test_external_field.py.
What is checked here is the read surface B0 added - the getB0{x,y,z} bindings and
their hierarchy_utils entries, B0x/B0y/B0z in gridlayout's yee_centering, Run.GetB0
and the plotting dispatch.

B0 is a good probe for a read path because its exact value is known analytically
everywhere. The updater builds it as the discrete curl of a vector potential, so

    A0 = (0, 0, b0x*y - b0y*x)   =>   B0 = (b0x, b0y, 0)

is constant over the whole domain, ghost nodes included, at any resolution and on
any level. Anything that misreads shapes, centerings or patch association breaks
that identity.
"""

import unittest

import numpy as np

import pyphare.pharein as ph
from pyphare.pharesee.hierarchy import hierarchy_from, VectorField
from pyphare.pharesee.plotting import finest_field_plot
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

diag_outputs = "phare_outputs/test_pharesee_external_field"
time_step = 0.001
cells = (20, 20)
dl = (0.2, 0.2)

b0x_const, b0y_const = 0.3, -0.2


def az_const(x, y):
    """A0 = (0, 0, b0x*y - b0y*x)  =>  B0 = (b0x, b0y, 0), constant everywhere."""
    return b0x_const * y - b0y_const * x


class ExternalFieldReadbackTest(unittest.TestCase):
    def diag_dir(self):
        return f"{diag_outputs}/{self._testMethodName}"

    def setUp(self):
        ph.global_vars.sim = None

        sim = ph.Simulation(
            time_step=time_step,
            time_step_nbr=1,
            cells=cells,
            dl=dl,
            refinement_boxes={"L0": {"B0": [(5, 5), (14, 14)]}},
            smallest_patch_size=10,
            largest_patch_size=20,
            diag_options={
                "format": "phareh5",
                "options": {"dir": self.diag_dir(), "mode": "overwrite"},
            },
            external_field={
                "type": "user-defined",
                "potential": (None, None, az_const),
            },
        )

        def density(x, y):
            return 1.0

        def bx(x, y):
            return 1.0

        def zero(x, y):
            return 0.0

        def vth(x, y):
            return 0.3

        ph.MaxwellianFluidModel(
            bx=bx,
            by=zero,
            bz=zero,
            protons={
                "charge": 1,
                "density": density,
                **{f"vbulk{c}": zero for c in "xyz"},
                **{f"vth{c}": vth for c in "xyz"},
                "nbr_part_per_cell": 10,
                "init": {"seed": 13337},
            },
        )
        ph.ElectronModel(closure="isothermal", Te=0.12)
        ph.ElectromagDiagnostics(quantity="B0", write_timestamps=np.array([0.0]))

        self.sim = sim
        self.simulator = Simulator(sim)
        self.simulator.initialize()

    def tearDown(self):
        if self.simulator is not None:
            self.simulator.reset()
        self.simulator = None
        ph.global_vars.sim = None

    def dump(self):
        """write the diagnostics of the live simulator, then return a Run on them"""
        self.simulator.dump(self.simulator.currentTime(), self.sim.time_step)
        return Run(self.diag_dir())

    def test_b0_read_in_memory_is_the_analytic_field(self):
        """the getB0{x,y,z} bindings, and their field_qties/quantidic entries"""
        for qty, value in {
            "EM_B0_x": b0x_const,
            "EM_B0_y": b0y_const,
            "EM_B0_z": 0.0,
        }.items():
            hier = hierarchy_from(simulator=self.simulator, qty=qty)
            levels = hier.levels()
            self.assertGreaterEqual(len(levels), 2, f"{qty}: {len(levels)} level(s)")
            for ilvl, level in levels.items():
                self.assertTrue(level.patches, f"{qty} lvl={ilvl}: no patch")
                for patch in level.patches:
                    data = patch.patch_datas[qty].dataset[:]
                    self.assertGreater(data.size, 0)
                    np.testing.assert_allclose(
                        data, value, atol=1e-12, rtol=0, err_msg=f"{qty} lvl={ilvl}"
                    )

    def test_get_b0_all_primal(self):
        """GetB0 with all_primal needs B0x/B0y/B0z in gridlayout's yee_centering.

        B0 is constant here, so averaging a dual direction onto primal nodes
        leaves it unchanged and the identity survives the conversion exactly.
        """
        b0 = self.dump().GetB0(0.0)  # all_primal=True -> a VectorField
        self.assertIsInstance(b0, VectorField)
        self.assertEqual(b0.quantities(), ["x", "y", "z"])

        expected = {"x": b0x_const, "y": b0y_const, "z": 0.0}
        checked = 0
        for ilvl, level in b0.levels(0.0).items():
            for patch in level.patches:
                for name, pd in patch.patch_datas.items():
                    data = pd.dataset[:]
                    # _compute_to_primal leaves NaN in the ghost layers by design
                    valid = np.isfinite(data)
                    self.assertTrue(valid.any(), f"lvl={ilvl} {name}: all NaN")
                    np.testing.assert_allclose(
                        data[valid],
                        expected[name],
                        atol=1e-12,
                        rtol=0,
                        err_msg=f"lvl={ilvl} {name}",
                    )
                    checked += valid.sum()
        self.assertGreater(checked, 0)

    def test_finest_field_plot_accepts_the_b0_components(self):
        """finest_field_plot dispatches B0x/B0y/B0z to EM_B0.h5 and GetB0"""
        import matplotlib.pyplot as plt

        run = self.dump()
        for qty in ("B0x", "B0y", "B0z"):
            fig, ax = finest_field_plot(run.path, qty, time=0.0)
            self.assertIsNotNone(fig)
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
