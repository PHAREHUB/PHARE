import unittest
from unittest import mock

import numpy as np

from pyphare.core import phare_utilities
from pyphare.pharein import global_vars, simulation, tagging
from pyphare.pharein.initialize import general, hybrid, mhd
from pyphare.pharein.mhd_model import MHDModel

from pyphare import pharein as ph


class TestSimulation(unittest.TestCase):
    def setUp(self):
        self.cells_array = [80, (80, 40), (80, 40, 12)]
        self.dl_array = [0.1, (0.1, 0.2), (0.1, 0.2, 0.3)]
        self.domain_size_array = [100.0, (100.0, 80.0), (100.0, 80.0, 20.0)]
        self.ndim = [1, 2]  # TODO https://github.com/PHAREHUB/PHARE/issues/232
        self.layout = "yee"
        self.time_step = 0.001
        self.time_step_nbr = 1000
        self.final_time = 1.0
        global_vars.sim = None

    def test_dl(self):
        for cells, domain_size, dim in zip(
            self.cells_array, self.domain_size_array, self.ndim
        ):
            j = simulation.Simulation(
                time_step_nbr=self.time_step_nbr,
                cells=cells,
                domain_size=domain_size,
                final_time=self.final_time,
            )

            if phare_utilities.none_iterable(domain_size, cells):
                domain_size = phare_utilities.listify(domain_size)
                cells = phare_utilities.listify(cells)

            for d in np.arange(dim):
                self.assertEqual(j.dl[d], domain_size[d] / float(cells[d]))

            global_vars.sim = None

    def test_boundaries_default_periodic(self):
        j = simulation.Simulation(
            time_step_nbr=1000,
            cells=80,
            domain_size=10,
            final_time=1.0,
        )

        for d in np.arange(j.ndim):
            self.assertTrue(j.periodicities[d])

    def test_boundary_types_kwarg_rejected(self):
        with self.assertRaises(ValueError):
            simulation.Simulation(
                time_step_nbr=1000,
                boundary_types="periodic",
                cells=80,
                domain_size=10,
                final_time=1000,
            )

    def test_time_step(self):
        s = simulation.Simulation(
            time_step_nbr=1000,
            cells=80,
            domain_size=10,
            final_time=10,
        )
        self.assertEqual(0.01, s.time_stepper.time_step)

    # ---- physical outer boundary conditions ----------------------------------

    def _mhd_kwargs(self, **overrides):
        kwargs = dict(
            time_step_nbr=1,
            cells=80,
            domain_size=10,
            final_time=1.0,
            model_options=["MHDModel"],
            mhd_timestepper="SSPRK4_5",
            reconstruction="WENOZ",
            limiter="None",
            riemann="Rusanov",
        )
        kwargs.update(overrides)
        return kwargs

    def test_boundaries_default_none(self):
        # without a domain_boundaries dict every location is periodic and defaults to 'none'
        global_vars.sim = None
        s = simulation.Simulation(**self._mhd_kwargs())
        self.assertEqual([True], s.periodicities)
        for loc in ("xlower", "xupper"):
            self.assertEqual("none", s.domain_boundaries[loc].type)

    def test_physical_direction_requires_both_locations(self):
        # giving only one side of a direction must be rejected
        global_vars.sim = None
        with self.assertRaises(ValueError):
            simulation.Simulation(
                **self._mhd_kwargs(domain_boundaries={"xlower": {"type": "open"}})
            )

    def test_inflow_velocity_scalar_normalized_signed(self):
        # a scalar inflow speed becomes the signed inward-normal component
        global_vars.sim = None
        s = simulation.Simulation(
            **self._mhd_kwargs(
                domain_boundaries={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "velocity": 2.0,
                        "density": 1.0,
                        "pressure": 1.0,
                        "B": [0.5, 1.0, 0.0],
                    },
                    "xupper": {"type": "open"},
                },
            )
        )
        self.assertEqual([False], s.periodicities)
        vx, vy, vz = s.domain_boundaries["xlower"].velocity
        self.assertEqual((2.0, 0.0, 0.0), (vx, vy, vz))  # +x inward at lower

    def test_inflow_scalar_B_rejected(self):
        # the magnetic field of an inflow must be a 3-vector, not a scalar
        global_vars.sim = None
        with self.assertRaises((TypeError, ValueError)):
            simulation.Simulation(
                **self._mhd_kwargs(
                    domain_boundaries={
                        "xlower": {
                            "type": "super-magnetofast-inflow",
                            "velocity": 2.0,
                            "density": 1.0,
                            "pressure": 1.0,
                            "B": 0.5,
                        },
                        "xupper": {"type": "open"},
                    },
                )
            )

    def _inflow_sim(self, velocity=2.0, B=(0.5, 1.0, 0.0), **overrides):
        global_vars.sim = None
        return simulation.Simulation(
            **self._mhd_kwargs(
                domain_boundaries={
                    "xlower": {
                        "type": "super-magnetofast-inflow",
                        "velocity": velocity,
                        "density": 1.0,
                        "pressure": 1.0,
                        "B": list(B),
                    },
                    "xupper": {"type": "open"},
                },
                **overrides,
            )
        )

    def test_inflow_consistent_normal_b_accepted(self):
        s = self._inflow_sim()
        MHDModel(bx=lambda x: 0.5 + 0.1 * np.sin(np.pi * x / 10.0), by=lambda x: 1.0)
        self.assertIs(s, global_vars.sim)
        self.assertIsNotNone(s.model)

    def test_inflow_inconsistent_normal_b_rejected(self):
        self._inflow_sim()
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x: 0.7 + x * 0, by=lambda x: 1.0)

    def test_inflow_oblique_inconsistent_normal_b_rejected(self):
        self._inflow_sim(velocity=[2.0, 0.5, 0.0])
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x: 0.7 + x * 0, by=lambda x: 1.0)

    def test_inflow_normal_b_checked_along_whole_face(self):
        self._inflow_sim(cells=(80, 40), domain_size=(10.0, 8.0))
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x, y: 0.5 + 0.01 * np.sin(2 * np.pi * y / 8.0))

    def test_inflow_normal_b_checked_at_upper_face(self):
        global_vars.sim = None
        simulation.Simulation(
            **self._mhd_kwargs(
                domain_boundaries={
                    "xlower": {"type": "open"},
                    "xupper": {
                        "type": "super-magnetofast-inflow",
                        "velocity": 2.0,
                        "density": 1.0,
                        "pressure": 1.0,
                        "B": [0.5, 1.0, 0.0],
                    },
                },
            )
        )
        with self.assertRaises(ValueError):
            MHDModel(bx=lambda x: np.where(x > 9.0, 0.6, 0.5), by=lambda x: 1.0)

    def tagging_simulation(self, **kwargs):
        return simulation.Simulation(
            time_step_nbr=self.time_step_nbr,
            cells=self.cells_array[0],
            domain_size=self.domain_size_array[0],
            final_time=self.final_time,
            max_nbr_levels=2,
            refinement="tagging",
            **kwargs,
        )

    def test_tagging_unknown_top_level_key_raises(self):
        with self.assertRaises(ValueError):
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={"method": "lohner-per-direction", "quantites": {"rho": 0.4}},
            )

    def test_tagging_without_refinement_raises(self):
        with self.assertRaises(ValueError):
            tagging.resolve_tagging(
                max_nbr_levels=3,
                tagging={"method": "lohner-per-direction", "quantities": {"B": 0.1}},
            )

    def test_tagging_valid_spec_normalizes(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "lohner-per-direction",
                "quantities": {"B": {"threshold": 0.3, "reltol": 0.02}, "rho": 0.4},
            },
        )
        self.assertIsInstance(tag, tagging.LohnerPerDirectionTagging)
        self.assertEqual(tag.method, "lohner-per-direction")
        self.assertEqual(
            tag.quantities,
            [
                ("B", {"threshold": 0.3, "reltol": 0.02, "ref_magnitude": 1e-8}),
                ("rho", {"threshold": 0.4, "reltol": 0.01, "ref_magnitude": 1e-8}),
            ],
        )

    def test_tagging_lohner_default_params(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={"method": "lohner-per-direction", "quantities": {"B": 0.1}},
        )
        self.assertEqual(
            tag.quantities,
            [("B", {"threshold": 0.1, "reltol": 0.01, "ref_magnitude": 1e-8})],
        )

    def test_tagging_lohner_cross_derivatives(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "lohner-cross-derivatives",
                "quantities": {"B": {"threshold": 0.2, "ref_magnitude": 1e-6}},
            },
        )
        self.assertIsInstance(tag, tagging.LohnerCrossDerivativesTagging)
        self.assertEqual(tag.method, "lohner-cross-derivatives")
        self.assertEqual(
            tag.quantities,
            [("B", {"threshold": 0.2, "reltol": 0.01, "ref_magnitude": 1e-6})],
        )

    def test_tagging_list_form(self):
        as_list = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "lohner-per-direction",
                "quantities": [("B", 0.1), ("rho", {"threshold": 0.2, "reltol": 0.0})],
            },
        )
        as_dict = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "lohner-per-direction",
                "quantities": {"B": 0.1, "rho": {"threshold": 0.2, "reltol": 0.0}},
            },
        )
        self.assertEqual(as_list.quantities, as_dict.quantities)

    def test_tagging_unknown_quantity_key_raises(self):
        with self.assertRaises(ValueError) as error:
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={"quantities": {"B": {"threshold": 0.1, "reltol": 0.1}}},
            )
        self.assertIn("reltol", str(error.exception))
        self.assertIn("threshold", str(error.exception))

    def test_tagging_dict_without_threshold_raises(self):
        with self.assertRaises(ValueError) as error:
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={
                    "method": "lohner-per-direction",
                    "quantities": {"B": {"reltol": 0.1}},
                },
            )
        self.assertIn("threshold", str(error.exception))

    def test_tagging_params_key_rejected(self):
        with self.assertRaises(ValueError):
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={"method": "lohner-per-direction", "params": {"reltol": 0.1}},
            )

    def test_tagging_lohner_negative_rejected(self):
        for key, value in (("reltol", -0.1), ("ref_magnitude", -1.0)):
            with self.assertRaises(ValueError) as error:
                tagging.resolve_tagging(
                    refinement="tagging",
                    tagging={
                        "method": "lohner-per-direction",
                        "quantities": {"B": {"threshold": 0.1, key: value}},
                    },
                )
            self.assertIn(key, str(error.exception))

    def test_tagging_non_finite_parameters_rejected(self):
        cases = [
            ("lohner-per-direction", "reltol", float("nan")),
            ("lohner-per-direction", "ref_magnitude", float("nan")),
            ("lohner-cross-derivatives", "ref_magnitude", float("inf")),
            ("gradient-scale", "ref_magnitude", float("nan")),
            ("gradient-scale", "ref_magnitude", float("inf")),
            ("default", "threshold", float("nan")),
        ]
        for method, key, value in cases:
            spec = (
                {"threshold": 0.1, "ref_magnitude": 0.1}
                if method == "gradient-scale"
                else {}
            )
            spec = {"threshold": 0.1, **spec, key: value}
            with self.subTest(method=method, key=key, value=value):
                with self.assertRaises(ValueError) as error:
                    tagging.resolve_tagging(
                        refinement="tagging",
                        tagging={"method": method, "quantities": {"B": spec}},
                    )
                self.assertIn(key, str(error.exception))

    def test_tagging_explicit_none_is_not_a_number(self):
        with self.assertRaises(TypeError) as error:
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={
                    "method": "lohner-per-direction",
                    "quantities": {"B": {"threshold": 0.1, "reltol": None}},
                },
            )
        self.assertIn("reltol must be a number", str(error.exception))

    def test_tagging_malformed_list_entries_rejected(self):
        for entry in ("Bx", ("B",), ("B", 0.1, 0.2), (1, 0.1)):
            with self.subTest(entry=entry):
                with self.assertRaises(TypeError) as error:
                    tagging.resolve_tagging(
                        refinement="tagging", tagging={"quantities": [entry]}
                    )
                self.assertIn("(name, threshold or dict)", str(error.exception))

    def test_tagging_string_thresholds_rejected(self):
        with self.assertRaises(TypeError):
            tagging.resolve_tagging(
                refinement="tagging", tagging={"quantities": {"B": "0.1"}}
            )
        with self.assertRaises(TypeError) as error:
            tagging.resolve_tagging(refinement="tagging", tagging_threshold="0.1")
        self.assertIn("tagging_threshold", str(error.exception))

    def test_tagging_numeric_types_become_float(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "lohner-per-direction",
                "quantities": {
                    "B": 1,
                    "rho": {"threshold": np.float32(0.25), "reltol": 0},
                },
            },
        )
        self.assertEqual(
            tag.quantities,
            [
                ("B", {"threshold": 1.0, "reltol": 0.01, "ref_magnitude": 1e-8}),
                ("rho", {"threshold": 0.25, "reltol": 0.0, "ref_magnitude": 1e-8}),
            ],
        )
        for _, parameters in tag.quantities:
            for value in parameters.values():
                self.assertIs(type(value), float)
        with self.assertRaises((TypeError, ValueError)):
            tagging.resolve_tagging(
                refinement="tagging", tagging={"quantities": {"B": True}}
            )

    def test_gradient_scale_spec(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "gradient-scale",
                "quantities": {
                    "B": {"threshold": 0.1, "ref_magnitude": 0.01},
                    "rho": {"threshold": 0.2, "ref_magnitude": 0.05},
                },
            },
        )
        self.assertIsInstance(tag, tagging.GradientScaleTagging)
        self.assertEqual(tag.method, "gradient-scale")
        self.assertEqual(
            tag.quantities,
            [
                ("B", {"threshold": 0.1, "ref_magnitude": 0.01}),
                ("rho", {"threshold": 0.2, "ref_magnitude": 0.05}),
            ],
        )

    def test_gradient_scale_ref_magnitude_required(self):
        with self.assertRaises(ValueError) as error:
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={"method": "gradient-scale", "quantities": {"B": 0.1}},
            )
        self.assertIn("ref_magnitude", str(error.exception))

    def test_gradient_scale_without_quantities_raises(self):
        with self.assertRaises(ValueError) as error:
            tagging.resolve_tagging(
                refinement="tagging", tagging={"method": "gradient-scale"}
            )
        self.assertIn("ref_magnitude", str(error.exception))
        self.assertIn("quantities", str(error.exception))
        self.assertNotIn("'B'", str(error.exception))

    def test_gradient_scale_ref_magnitude_must_be_positive(self):
        for ref_magnitude in (0.0, -1.0):
            with self.assertRaises(ValueError) as error:
                tagging.resolve_tagging(
                    refinement="tagging",
                    tagging={
                        "method": "gradient-scale",
                        "quantities": {
                            "B": {"threshold": 0.1, "ref_magnitude": ref_magnitude}
                        },
                    },
                )
            self.assertIn("ref_magnitude", str(error.exception))

    def test_gradient_scale_rejects_reltol(self):
        with self.assertRaises(ValueError) as error:
            tagging.resolve_tagging(
                refinement="tagging",
                tagging={
                    "method": "gradient-scale",
                    "quantities": {
                        "B": {"threshold": 0.1, "ref_magnitude": 0.1, "reltol": 0.1}
                    },
                },
            )
        self.assertIn("reltol", str(error.exception))

    def test_tagging_old_lohner_name_rejected(self):
        with self.assertRaises(ValueError) as error:
            tagging.resolve_tagging(refinement="tagging", tagging={"method": "lohner"})
        self.assertIn("lohner-per-direction", str(error.exception))

    def test_tagging_none_when_not_tagging(self):
        self.assertIsNone(tagging.resolve_tagging())

    def test_tagging_unknown_method_raises(self):
        with self.assertRaises(ValueError):
            tagging.resolve_tagging(refinement="tagging", tagging={"method": "wavelet"})

    def test_tagging_not_available_rejected(self):
        for name in ("J", "Ey", "V", "pressure", "massDensity"):
            with self.assertRaises(ValueError) as error:
                tagging.resolve_tagging(
                    refinement="tagging", tagging={"quantities": {name: 0.1}}
                )
            self.assertIn("not available for tagging yet", str(error.exception))

    def test_tagging_duplicate_quantity_rejected(self):
        with self.assertRaises(ValueError) as error:
            tagging.resolve_tagging(
                refinement="tagging", tagging={"quantities": [("B", 0.1), ("B", 0.4)]}
            )
        self.assertIn("given more than once", str(error.exception))

    def test_tagging_fallback_uses_tagging_threshold(self):
        expected = [
            (
                {"method": "lohner-per-direction"},
                {"reltol": 0.01, "ref_magnitude": 1e-8},
            ),
            ({"quantities": {}}, {}),
            (None, {}),
        ]
        for spec, parameters in expected:
            tag = tagging.resolve_tagging(
                refinement="tagging", tagging_threshold=0.3, tagging=spec
            )
            self.assertEqual(tag.quantities, [("B", {"threshold": 0.3, **parameters})])

    def test_re_numpify_does_not_corrupt_tagging_quantities(self):
        quantities = [
            ("B", {"threshold": 0.1, "reltol": 0.01, "ref_magnitude": 1e-8}),
            ("rho", {"threshold": 0.4, "reltol": 0.02, "ref_magnitude": 1e-8}),
        ]
        tag = tagging.LohnerPerDirectionTagging(quantities=[*quantities])
        simulation.re_numpify_recursive(tag, "quantities")
        self.assertEqual(tag.quantities, quantities)
        self.assertIsInstance(tag.quantities, list)
        self.assertIsInstance(tag.quantities[0], tuple)
        self.assertIs(type(tag.quantities[0][1]["threshold"]), float)

    def test_re_numpify_still_numpifies_homogeneous_numeric_lists(self):
        owner = simulation.Simulation.__new__(simulation.Simulation)
        object.__setattr__(owner, "cells", [80, 40, 12])
        simulation.re_numpify_recursive(owner, "cells")
        self.assertIsInstance(owner.cells, np.ndarray)
        np.testing.assert_array_equal(owner.cells, np.array([80, 40, 12]))

    def test_tagging_serialize_deserialize_roundtrip_preserves_quantities(self):
        sim = self.tagging_simulation(
            tagging={
                "method": "lohner-per-direction",
                "quantities": {"B": 0.1, "rho": {"threshold": 0.4, "reltol": 0.02}},
            },
        )
        restored = simulation.deserialize(simulation.serialize(sim))
        self.assertEqual(
            restored.tagging.quantities,
            [
                ("B", {"threshold": 0.1, "reltol": 0.01, "ref_magnitude": 1e-8}),
                ("rho", {"threshold": 0.4, "reltol": 0.02, "ref_magnitude": 1e-8}),
            ],
        )
        for _, parameters in restored.tagging.quantities:
            for value in parameters.values():
                self.assertIs(type(value), float)

    def test_default_path_populates_legacy_keys(self):
        sim = self.tagging_simulation()
        dp = FakeDictPopulator()
        general.populate_amr(dp, sim)

        base = "simulation/AMR/refinement/tagging/"
        tagging_calls = {
            path[len(base) :]: value
            for path, value in dp.calls.items()
            if path.startswith(base)
        }
        self.assertEqual(
            tagging_calls,
            {
                "method": "auto",
                "threshold": 0.1,
                "nbr_quantities": 1,
                "Q0/name": "B",
                "Q0/threshold": 0.1,
            },
        )
        self.assertEqual(sim.tagging.method, "default")

    def model_method(self, module, key):
        dp = FakeDictPopulator()
        with mock.patch.object(module, "DictPopulator", lambda: dp):
            module.populateDict(global_vars.sim)
        return dp.calls[f"simulation/AMR/refinement/tagging/{key}"]

    def test_method_reaches_model_method_keys(self):
        def one(x):
            return 1.0 + 0.0 * x

        for method in (
            "default",
            "lohner-per-direction",
            "lohner-cross-derivatives",
        ):
            global_vars.sim = None
            self.tagging_simulation(tagging={"method": method})
            ph.MaxwellianFluidModel(bx=one, by=one, bz=one, protons={"density": one})
            ph.ElectronModel(closure="isothermal", Te=0.0)
            self.assertEqual(self.model_method(hybrid, "hybrid_method"), method)

            global_vars.sim = None
            self.tagging_simulation(
                tagging={"method": method},
                model_options=["MHDModel"],
                reconstruction="Linear",
                limiter="VanLeer",
                riemann="Rusanov",
                mhd_timestepper="TVDRK2",
            )
            ph.MHDModel(
                density=one, vx=one, vy=one, vz=one, bx=one, by=one, bz=one, p=one
            )
            self.assertEqual(self.model_method(mhd, "mhd_method"), method)


class FakeDictPopulator:
    def __init__(self):
        self.calls = {}

    def _record(self, path, value):
        assert path not in self.calls, f"duplicate path populated: {path}"
        self.calls[path] = value

    def __getattr__(self, name):
        if name in ("add_string", "add_int", "add_double", "add_bool"):
            return self._record
        if name.startswith("add_"):
            return lambda *args: None
        raise AttributeError(name)


class TestTaggingPopulateDict(unittest.TestCase):
    def test_gradient_scale_populate_dict(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "gradient-scale",
                "quantities": {"rho": {"threshold": 0.2, "ref_magnitude": 0.05}},
            },
        )
        dp = FakeDictPopulator()
        tag.populate_dict(dp)

        base = "simulation/AMR/refinement/tagging/"
        self.assertEqual(dp.calls[base + "Q0/name"], "rho")
        self.assertEqual(dp.calls[base + "Q0/threshold"], 0.2)
        self.assertEqual(dp.calls[base + "Q0/ref_magnitude"], 0.05)
        self.assertNotIn(base + "Q0/reltol", dp.calls)

    def test_lohner_multi_quantity_populate_dict(self):
        tag = tagging.resolve_tagging(
            refinement="tagging",
            tagging={
                "method": "lohner-per-direction",
                "quantities": {
                    "B": {"threshold": 0.1, "reltol": 0.02, "ref_magnitude": 1e-20},
                    "rho": 0.4,
                },
            },
        )
        dp = FakeDictPopulator()
        tag.populate_dict(dp)

        base = "simulation/AMR/refinement/tagging/"
        self.assertEqual(dp.calls[base + "nbr_quantities"], 2)
        self.assertEqual(dp.calls[base + "Q0/name"], "B")
        self.assertEqual(dp.calls[base + "Q0/threshold"], 0.1)
        self.assertEqual(dp.calls[base + "Q1/name"], "rho")
        self.assertEqual(dp.calls[base + "Q1/threshold"], 0.4)
        self.assertEqual(dp.calls[base + "Q0/reltol"], 0.02)
        self.assertEqual(dp.calls[base + "Q0/ref_magnitude"], 1e-20)
        self.assertEqual(dp.calls[base + "Q1/reltol"], 0.01)
        self.assertEqual(dp.calls[base + "Q1/ref_magnitude"], 1e-8)
        self.assertFalse(any("params" in key for key in dp.calls))
        self.assertNotIn(base + "method", dp.calls)


if __name__ == "__main__":
    unittest.main()
