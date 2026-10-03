import unittest
from ddt import ddt
import numpy as np

from pyphare.core.box import Box
from pyphare.pharesee.run import Run
from pyphare.simulator.simulator import Simulator
from pyphare.pharesee.hierarchy import PatchHierarchy
from pyphare.pharesee.hierarchy import ScalarField, VectorField
from pyphare.core.operators import dot, cross, sqrt, modulus, grad

diag_outputs = "phare_outputs/"
time_step_nbr = 20
time_step = 0.005
final_time = time_step * time_step_nbr
dt = 10 * time_step
nt = int(final_time / dt) + 1
timestamps = dt * np.arange(nt)


def patch_box_data(patch, qty):
    pd = patch[qty]
    return pd[patch.box]


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

            for quantity in ["density", "flux"]:
                ph.FluidDiagnostics(
                    quantity=quantity,
                    write_timestamps=timestamps,
                    population_name="protons",
                )

            return sim

        Simulator(config()).run()

    def _test_data_is_a_hierarchy(self):
        r = Run(self.diag_dir())
        B = r.GetB(0.0)
        self.assertTrue(isinstance(B, PatchHierarchy))

    def _test_can_read_multiple_times(self):
        r = Run(self.diag_dir())
        times = (0.0, 0.1)
        B = r.GetB(times)
        E = r.GetE(times)
        Ni = r.GetNi(times)
        Vi = r.GetVi(times)
        for hier in (B, E, Ni, Vi):
            self.assertEqual(len(hier.times()), 2)
            self.assertTrue(np.allclose(hier.times().astype(np.float32), times))

    def _test_hierarchy_is_refined(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertEqual(len(B.levels()), B.levelNbr())
        self.assertEqual(len(B.levels()), 2)
        self.assertEqual(len(B.levels(time)), 2)
        self.assertEqual(len(B.levels(time)), B.levelNbr(time))

    def _test_can_get_nbytes(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertGreater(B.nbytes(), 0)

    def _test_hierarchy_has_patches(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertGreater(B.nbrPatches(), 0)

    def _test_access_patchdatas_as_hierarchies(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        self.assertTrue(isinstance(B.x, PatchHierarchy))
        self.assertTrue(isinstance(B.y, PatchHierarchy))
        self.assertTrue(isinstance(B.z, PatchHierarchy))

    def _test_partial_domain_hierarchies(self):
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

    def _test_scalarfield_quantities(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Vi = r.GetVi(time)
        self.assertTrue(isinstance(Ni, ScalarField))
        self.assertTrue(isinstance(Vi, VectorField))

    def _test_sum_two_scalarfields(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Pe = r.GetPe(time)
        s = Ni + Pe
        self.assertTrue(isinstance(s, ScalarField))
        self.assertEqual(s.quantities(), ["value"])

    def _test_sum_with_scalar(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        s1 = Ni + 0.1
        s2 = 0.1 + Ni
        for s in (s1, s2):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def _test_scalarfield_difference(self):
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

    def _test_scalarfield_product(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        Pe = r.GetPe(time)
        s1 = Ni * Pe
        s2 = Ni * 0.1
        for s in (s1, s2):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def _test_scalarfield_division(self):
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

    def _test_scalarfield_sqrt(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        nisqrt = sqrt(Ni)
        self.assertTrue(isinstance(nisqrt, ScalarField))

    def _test_vectorfield_dot_product(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        s1 = dot(B, B)
        s2 = modulus(B)
        for s in (s1, s2):
            self.assertTrue(isinstance(s, ScalarField))
            self.assertEqual(s.quantities(), ["value"])

    def _test_vectorfield_binary_ops(self):
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

    def _test_scalarfield_to_vectorfield_ops(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        s1 = grad(Ni)
        for s in (s1,):
            self.assertTrue(isinstance(s, VectorField))
            self.assertEqual(s.quantities(), ["x", "y", "z"])

    def _test_vectorfield_to_vectorfield_ops(self):
        r = Run(self.diag_dir())
        time = 0.0
        B = r.GetB(time)
        E = r.GetE(time)
        s1 = cross(E, B)
        for s in (s1,):
            self.assertTrue(isinstance(s, VectorField))
            self.assertEqual(s.quantities(), ["x", "y", "z"])

    def _test_mixed_ghost_getters(self):
        """
        single population of charge 1, so Ni == N and Vi == Flux / N
        some getters drop ghosts and others do not, they must still be combinable
        """
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        N = r.GetN(time, "protons")
        Vi = r.GetVi(time)
        Flux = r.GetFlux(time, "protons")
        B = r.GetB(time)

        diff = Ni - N
        V = Flux / N
        self.assertTrue(isinstance(np.add(Ni, N), ScalarField))
        self.assertTrue(isinstance(dot(Flux, B), ScalarField))

        # Ez is primal in 2d, E keeps its ghosts with all_primal=False
        Ez = ScalarField(r.GetE(time, all_primal=False).Ez)
        self.assertTrue(isinstance(Ni * Ez, ScalarField))
        self.assertTrue(isinstance(np.multiply(Ni, Ez), ScalarField))

        for ilvl, lvl in diff.levels(time).items():
            for ip, patch in enumerate(lvl.patches):
                pd = patch["value"]
                np.testing.assert_allclose(pd[patch.box], 0, atol=1e-12)

                # Vi on patch borders can be refilled after it is computed from
                # the pop moments (fillIonBorders), so only compare inner nodes
                inner = (slice(1, -1),) * 2
                vi_patch = Vi.level(ilvl, time).patches[ip]
                for c in ["x", "y", "z"]:
                    np.testing.assert_allclose(
                        patch_box_data(V.level(ilvl, time).patches[ip], c)[inner],
                        patch_box_data(vi_patch, c)[inner],
                        atol=1e-12,
                    )

    def _test_grad_values(self):
        r = Run(self.diag_dir())
        time = 0.0
        for scalar in (r.GetNi(time), r.GetPe(time)):
            g = grad(scalar)
            for ilvl, lvl in scalar.levels(time).items():
                for ip, patch in enumerate(lvl.patches):
                    data = patch_box_data(patch, "value")
                    expected = np.gradient(data)
                    gpatch = g.level(ilvl, time).patches[ip]
                    for i, c in enumerate(["x", "y"]):
                        computed = patch_box_data(gpatch, c)
                        self.assertFalse(np.isnan(computed).any())
                        # np.gradient is one sided on the edges of the patch data
                        inner = (slice(1, -1),) * 2
                        np.testing.assert_allclose(
                            computed[inner], expected[i][inner], atol=1e-12
                        )

    def _test_patch_cut(self):
        r = Run(self.diag_dir())
        time = 0.0
        for hier, qty in ((r.GetNi(time), "value"), (r.GetB(time), "x")):
            for patch in hier.level(0, time).patches:
                pd = patch[qty]
                data = pd[patch.box]
                idx = data.shape[0] // 2
                cut = pd.origin[0] + (idx + 0.1) * pd.layout.dl[0]
                np.testing.assert_array_equal(patch(qty, x=cut), data[idx, :])

    def _test_hierarchy_cut(self):
        r = Run(self.diag_dir())
        time = 0.0
        Ni = r.GetNi(time)
        dl = Ni.level(0, time).patches[0]["value"].layout.dl
        L = 86 * dl[1]

        coords, cut = Ni(x=12.95)  # not on a patch border
        self.assertEqual(coords.shape, cut.shape)
        self.assertFalse(np.isnan(cut).any())
        self.assertTrue((np.diff(coords) > 0).all())  # no overlaps across patches
        self.assertLessEqual(np.diff(coords).max(), dl[1] + 1e-12)  # no gaps
        self.assertAlmostEqual(coords[0], 0)
        self.assertAlmostEqual(coords[-1], L)

    def _test_finest_ignores_ghosts(self):
        """
        ghost values are not guaranteed to be valid, so they are poisoned
        with NaN which must not appear in the finest data
        """
        from pyphare.pharesee.hierarchy.hierarchy_utils import flat_finest_field

        r = Run(self.diag_dir())
        time = 0.0
        # GetN drops ghosts, so read the file as is
        N = ScalarField(r._get_hierarchy(time, "ions_pop_protons_density.h5"))
        for lvl in N.levels(time).values():
            for patch in lvl.patches:
                pd = patch["value"]
                self.assertTrue(all(g > 0 for g in pd.ghosts_nbr))
                domain = np.array(pd[pd.box])
                pd.dataset = np.full(pd.dataset.shape, np.nan)
                pd[pd.box] = domain

        data, _ = flat_finest_field(N, "value", time=time)
        self.assertFalse(np.isnan(data).any())

        finest = N.finest(time)
        self.assertFalse(np.isnan(finest[finest.box]).any())

    def test_all(self):
        """
        DO NOT RUN MULTIPLE SIMULATIONS!
        """

        checks = 0
        for test in [method for method in dir(self) if method.startswith("_test_")]:
            getattr(self, test)()
            checks += 1
        self.assertEqual(checks, 23)  # update if you add new tests


if __name__ == "__main__":
    unittest.main()
