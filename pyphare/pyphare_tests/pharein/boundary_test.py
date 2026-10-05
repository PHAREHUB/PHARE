import unittest

from pyphare.pharein import boundary


class TestBoundaryStructural(unittest.TestCase):
    def test_default_all_periodic(self):
        periodicities, resolved = boundary.resolve_boundaries(2)
        self.assertEqual([True, True], periodicities)
        for loc in ("xlower", "xupper", "ylower", "yupper"):
            self.assertIsInstance(resolved[loc], boundary.NoneBoundary)
            self.assertEqual("none", resolved[loc].type)

    def test_one_given_location_makes_direction_physical(self):
        periodicities, resolved = boundary.resolve_boundaries(
            2,
            domain_boundaries={"xlower": {"type": "open"}, "xupper": {"type": "open"}},
            model_options=["MHDModel"],
        )
        self.assertEqual([False, True], periodicities)
        self.assertIsInstance(resolved["ylower"], boundary.NoneBoundary)

    def test_missing_opposite_location_raises(self):
        for given in ("xlower", "xupper"):
            with self.assertRaises(ValueError):
                boundary.resolve_boundaries(
                    2,
                    domain_boundaries={given: {"type": "open"}},
                    model_options=["MHDModel"],
                )

    def test_none_type_rejected(self):
        with self.assertRaises(ValueError):
            boundary.resolve_boundaries(
                2,
                domain_boundaries={"xlower": {"type": "none"}, "xupper": {"type": "open"}},
                model_options=["MHDModel"],
            )

    def test_missing_type_raises(self):
        with self.assertRaises(KeyError):
            boundary.resolve_boundaries(
                2,
                domain_boundaries={"xlower": {}, "xupper": {"type": "open"}},
                model_options=["MHDModel"],
            )

    def test_physical_only_supported_by_mhd_model(self):
        with self.assertRaises(ValueError):
            boundary.resolve_boundaries(
                2,
                domain_boundaries={"xlower": {"type": "open"}, "xupper": {"type": "open"}},
                model_options=["HybridModel"],
            )

    def test_unknown_location_rejected(self):
        with self.assertRaises(ValueError):
            boundary.resolve_boundaries(
                2,
                domain_boundaries={"not_a_location": {"type": "open"}},
                model_options=["MHDModel"],
            )

    def test_location_beyond_dimension_rejected(self):
        with self.assertRaises(ValueError):
            boundary.resolve_boundaries(
                1,
                domain_boundaries={"ylower": {"type": "open"}, "yupper": {"type": "open"}},
                model_options=["MHDModel"],
            )

    def test_open_and_reflective_resolve(self):
        _, resolved = boundary.resolve_boundaries(
            2,
            domain_boundaries={"xlower": {"type": "open"}, "xupper": {"type": "reflective"}},
            model_options=["MHDModel"],
        )
        self.assertIsInstance(resolved["xlower"], boundary.OpenBoundary)
        self.assertIsInstance(resolved["xupper"], boundary.ReflectiveBoundary)


class TestInflowOutflowData(unittest.TestCase):
    def _resolve(self, **bcs):
        return boundary.resolve_boundaries(
            2,
            domain_boundaries=bcs,
            model_options=["MHDModel"],
        )[1]

    def test_inflow_velocity_scalar_normalized_signed(self):
        resolved = self._resolve(
            xlower={
                "type": "super-magnetofast-inflow",
                "velocity": 2.0,
                "density": 1.0,
                "pressure": 1.0,
                "B": [0.5, 1.0, 0.0],
            },
            xupper={"type": "open"},
        )
        self.assertEqual((2.0, 0.0, 0.0), resolved["xlower"].velocity)

    def test_inflow_scalar_B_rejected(self):
        with self.assertRaises((TypeError, ValueError)):
            self._resolve(
                xlower={
                    "type": "super-magnetofast-inflow",
                    "velocity": 2.0,
                    "density": 1.0,
                    "pressure": 1.0,
                    "B": 0.5,
                },
                xupper={"type": "open"},
            )

    def test_inflow_missing_parameter_raises_keyerror(self):
        with self.assertRaises(KeyError):
            self._resolve(
                xlower={
                    "type": "super-magnetofast-inflow",
                    "velocity": 2.0,
                    "density": 1.0,
                    "B": [0.5, 1.0, 0.0],
                },
                xupper={"type": "open"},
            )

    def test_inflow_unknown_parameter_rejected(self):
        with self.assertRaises(ValueError):
            self._resolve(
                xlower={
                    "type": "super-magnetofast-inflow",
                    "velocity": 2.0,
                    "density": 1.0,
                    "pressure": 1.0,
                    "B": [0.5, 1.0, 0.0],
                    "temperature": 1.0,
                },
                xupper={"type": "open"},
            )

    def test_callable_inflow_value_rejected(self):
        with self.assertRaises(ValueError):
            self._resolve(
                xlower={
                    "type": "super-magnetofast-inflow",
                    "velocity": 2.0,
                    "density": 1.0,
                    "pressure": 1.0,
                    "B": [lambda x, y, t: 0.5, 1.0, 0.0],
                },
                xupper={"type": "open"},
            )

    def test_parameter_rejected_on_parameterless_type(self):
        with self.assertRaises(ValueError):
            self._resolve(
                xlower={"type": "open", "density": 1.0},
                xupper={"type": "open"},
            )


class TestBoundaryTypeEnum(unittest.TestCase):
    def test_every_type_maps_to_a_cpp_enum_member(self):
        from pyphare.cpp import cpp_etc_lib

        members = cpp_etc_lib().BoundaryType.__members__
        for boundary_type in boundary._type_to_class:
            self.assertIn(boundary.boundary_type_enum_member(boundary_type), members)


if __name__ == "__main__":
    unittest.main()
