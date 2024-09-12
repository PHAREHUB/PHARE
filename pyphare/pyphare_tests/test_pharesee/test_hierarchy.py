import unittest
from ddt import ddt
import numpy as np

from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.run import Run
from pyphare.pharesee.hierarchy import PatchHierarchy
from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.core.operators import dot, cross, sqrt, modulus, grad
from pyphare.core.box import Box

diag_outputs = "phare_outputs/"
time_step_nbr = 20
time_step = 0.005
final_time = time_step * time_step_nbr
dt = 10 * time_step
nt = int(final_time / dt) + 1
timestamps = dt * np.arange(nt)


@ddt
class PatchHierarchyTest(unittest.TestCase):
    def diag_dir(self):
        return diag_outputs + self._testMethodName

    def setUp(self):
        import pyphare.pharein as ph

        ph.global_vars.sim = None

        def config():
            sim = ph.Simulation(
                time_step_nbr=time_step_nbr,
                time_step=time_step,
                cells=(86, 86),
                dl=(0.3, 0.3),
                refinement="tagging",
                max_nbr_levels=2,
                hyper_resistivity=0.005,
                resistivity=0.001,
                diag_options={
                    "format": "phareh5",
                    "options": {"dir": self.diag_dir(), "mode": "overwrite"},
                },
            )

            def density(x, y):
                L = sim.simulation_domain()[1]
                return (
                    0.2
                    + 1.0 / np.cosh((y - L * 0.3) / 0.5) ** 2
                    + 1.0 / np.cosh((y - L * 0.7) / 0.5) ** 2
                )

            def S(y, y0, l):
                return 0.5 * (1.0 + np.tanh((y - y0) / l))

            def by(x, y):
                Lx = sim.simulation_domain()[0]
                Ly = sim.simulation_domain()[1]
                w1 = 0.2
                w2 = 1.0
                x0 = x - 0.5 * Lx
                y1 = y - 0.3 * Ly
                y2 = y - 0.7 * Ly
                w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
                w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
                w5 = 2.0 * w1 / w2
                return (w5 * x0 * w3) + (-w5 * x0 * w4)

            def bx(x, y):
                Lx = sim.simulation_domain()[0]
                Ly = sim.simulation_domain()[1]
                w1 = 0.2
                w2 = 1.0
                x0 = x - 0.5 * Lx
                y1 = y - 0.3 * Ly
                y2 = y - 0.7 * Ly
                w3 = np.exp(-(x0 * x0 + y1 * y1) / (w2 * w2))
                w4 = np.exp(-(x0 * x0 + y2 * y2) / (w2 * w2))
                w5 = 2.0 * w1 / w2
                v1 = -1
                v2 = 1.0
                return (
                    v1
                    + (v2 - v1) * (S(y, Ly * 0.3, 0.5) - S(y, Ly * 0.7, 0.5))
                    + (-w5 * y1 * w3)
                    + (+w5 * y2 * w4)
                )

            def bz(x, y):
                return 0.0

            def b2(x, y):
                return bx(x, y) ** 2 + by(x, y) ** 2 + bz(x, y) ** 2

            def T(x, y):
                K = 1
                temp = 1.0 / density(x, y) * (K - b2(x, y) * 0.5)
                assert np.all(temp > 0)
                return temp

            def vx(x, y):
                return 0.0

            def vy(x, y):
                return 0.0

            def vz(x, y):
                return 0.0

            def vthx(x, y):
                return np.sqrt(T(x, y))

            def vthy(x, y):
                return np.sqrt(T(x, y))

            def vthz(x, y):
                return np.sqrt(T(x, y))

            vvv = {
                "vbulkx": vx,
                "vbulky": vy,
                "vbulkz": vz,
                "vthx": vthx,
                "vthy": vthy,
                "vthz": vthz,
                "nbr_part_per_cell": 100,
            }

            ph.MaxwellianFluidModel(
                bx=bx,
                by=by,
                bz=bz,
                protons={
                    "charge": 1,
                    "density": density,
                    **vvv,
                    "init": {"seed": 12334},
                },
            )

            ph.ElectronModel(closure="isothermal", Te=0.1)

            for quantity in ["E", "B"]:
                ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)

            for quantity in ["charge_density", "bulkVelocity"]:
                ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)

            return sim

        Simulator(config()).run()

    def test_data_is_a_hierarchy(self):
        r = Run(self.diag_dir())
        B = r.GetB(0.0)
        self.assertTrue(isinstance(B, PatchHierarchy))

    def test_can_read_multiple_times(self):
        r = Run(self.diag_dir())
        times = (0.0, 0.1)
        B = r.GetB(times)
        E = r.GetE(times)
        Ni = r.GetNi(times)
        Vi = r.GetVi(times)
        for hier in (B, E, Ni, Vi):
            self.assertEqual(len(hier.times()), 2)
            self.assertTrue(np.allclose(hier.times().astype(np.float32), times))

    def test_hierarchy_is_refined(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertEqual(len(B.levels()), B.levelNbr())
        self.assertEqual(len(B.levels()), 2)
        self.assertEqual(len(B.levels(time)), 2)
        self.assertEqual(len(B.levels(time)), B.levelNbr(time))

    def test_can_get_nbytes(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertGreater(B.nbytes(), 0)

    def test_hierarchy_has_patches(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertGreater(B.nbrPatches(), 0)

    def test_access_patchdatas_as_hierarchies(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertTrue(isinstance(B.x, PatchHierarchy))
        self.assertTrue(isinstance(B.y, PatchHierarchy))
        self.assertTrue(isinstance(B.z, PatchHierarchy))

    def test_partial_domain_hierarchies(self):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle

        r = Run(self.diag_dir())
        time = 0.0
        box = Box((10, 5), (18, 12.5))
        B = r.GetB(time)
        Bpartial = r.GetB(time, selection_box=box)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        B.x.plot(plot_patches=True, ax=ax1)
        Bpartial.x.plot(plot_patches=True, ax=ax2)
        ax2.add_patch(
            Rectangle(
                box.lower,
                box.shape[0],
                box.shape[1],
                ec="w",
                fc="none",
                lw=3,
            )
        )
        ax2.set_title(f"{self.id()}")
        fig.savefig(f"{self.id()}.png")

        self.assertTrue(isinstance(Bpartial, PatchHierarchy))
        self.assertLess(Bpartial.nbrPatches(), B.nbrPatches())

    def test_scalarfield_quantities(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Vi = r.GetVi(time)
        self.assertTrue(isinstance(Ni, ScalarField))
        self.assertTrue(isinstance(Vi, VectorField))

    def test_sum_two_scalarfields(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Pe = r.GetPe(time)
        s = Ni + Pe
        self.assertTrue(isinstance(s, ScalarField))
        self.assertEqual(s.quantities(), ["value"])

    def test_sum_with_scalar(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        s1 = Ni + 0.1
        s2 = 0.1 + Ni
        for s in (s1, s2):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def test_scalarfield_difference(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Pe = r.GetPe(time)
        s1 = Ni - Pe
        s2 = -Ni
        s3 = Ni - 0.1
        for s in (s1, s2, s3):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def test_scalarfield_product(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Pe = r.GetPe(time)
        s1 = Ni * Pe
        s2 = Ni * 0.1
        for s in (s1, s2):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def test_scalarfield_division(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Pe = r.GetPe(time)
        s1 = Ni / Pe
        s2 = Ni / 0.1
        s3 = Ni / Ni
        s4 = 0.1 / Ni
        for s in (s1, s2, s3, s4):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def test_scalarfield_sqrt(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        nisqrt = sqrt(Ni)
        self.assertTrue(isinstance(nisqrt, ScalarField))

    def test_vectorfield_dot_product(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        s1 = dot(B, B)
        s2 = modulus(B)
        for s in (s1, s2):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def test_vectorfield_binary_ops(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        V = r.GetVi(time)
        s1 = B + V
        s2 = B - V
        # s3 = -B  # TODO
        s4 = 2 * B
        s5 = B * 2
        s6 = B / 10
        for s in (s1, s2, s4, s5, s6):
            self.assertTrue(isinstance(s, VectorField))
            self.assertEqual(s.quantities(), ["x", "y", "z"])

    def test_scalarfield_to_vectorfield_ops(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        s1 = grad(Ni)
        for s in (s1,):
            self.assertTrue(isinstance(s, VectorField))
            self.assertEqual(s.quantities(), ["x", "y", "z"])

    def test_vectorfield_to_vectorfield_ops(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        E = r.GetE(time)
        s1 = cross(E, B)
        for s in (s1,):
            self.assertTrue(isinstance(s, VectorField))
            self.assertEqual(s.quantities(), ["x", "y", "z"])


if __name__ == "__main__":
    unittest.main()
