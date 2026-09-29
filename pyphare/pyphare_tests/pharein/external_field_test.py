import unittest

import numpy as np
from ddt import data, ddt, unpack
import pyphare.pharein.global_vars as global_vars
from pyphare.pharein import simulation
from pyphare.pharein.external_field import (
    DipoleExternalField,
    ExternalField,
    ZeroExternalField,
    UserDefinedExternalField,
    resolve_external_field,
)


def az_2d(x, y):
    return x + y


def az_2d_t(x, y, t):
    return (x + y) * t


def dazdt_2d(x, y, t):
    return x + y


def az_3d(x, y, z):
    return x + y + z


def ay_1d(x):
    return x


def az_2d_keyword_only(x, *, y):
    return x


def dazdt_2d_keyword_only(x, y, *, t):
    return x + y


def resolve_user_defined(ndim, potential, **extra):
    """Resolve a user-defined declaration, narrowed to what the tests then read."""
    ef = resolve_external_field(
        ndim,
        external_field={"type": "user-defined", "potential": potential, **extra},
    )
    assert isinstance(ef, UserDefinedExternalField)
    return ef


_dipole = {
    "type": "dipole",
    "position": (0.5, 1.5),
    "moment": (0.0, 2.0),
    "radius": 0.0,
}
_static = {"type": "user-defined", "potential": (None, None, az_2d)}
_timed = {
    "type": "user-defined",
    "potential": (None, None, az_2d_t),
    "potential_time_derivative": (None, None, dazdt_2d),
}

_invalid_declarations = [
    (2, {"type": "quadrupole"}),  # unknown type
    (2, {"position": (0.5, 1.5)}),  # no type
    (2, "dipole"),  # not a dict
    (1, _dipole),  # a dipole makes no sense in 1D
    (3, _dipole),  # both vectors must have ndim components
    (2, {**_dipole, "moment": (0.0, 1.0, 2.0)}),  # three components in 2D
    (2, {**_dipole, "moment": (0.0, 0.0)}),  # a zero moment is no dipole
    (2, {k: v for k, v in _dipole.items() if k != "moment"}),  # missing key
    (2, {k: v for k, v in _dipole.items() if k != "radius"}),  # radius is required
    (2, {**_dipole, "value": 3}),  # unknown key
    (2, {**_dipole, "radius": -1.0}),  # a negative radius
    (2, {**_dipole, "radius": float("nan")}),  # not a number
    (2, {**_dipole, "radius": float("inf")}),  # an infinite radius
    (2, {**_dipole, "radius": (1.0, 1.0)}),  # not a scalar
    (2, {**_dipole, "radius": "1"}),  # not a number
    #
    (2, {"type": "user-defined"}),  # no potential
    (2, {**_static, "potential": (None, None, None)}),  # all components None
    (2, {**_static, "potential": (None, az_2d)}),  # not three components
    (2, {**_static, "potential": (None, None, 1.0)}),  # not a callable
    (2, {**_static, "potential": (az_2d, None, az_2d_t)}),  # mixed signatures
    (3, _static),  # signature does not match ndim
    (1, {**_static, "potential": (ay_1d, None, None)}),  # a_x is unused in 1D
    (2, {**_static, "moment": (0.0, 1.0)}),  # unknown key
    # the time dependence of the potential and of its derivative must agree
    (2, {k: v for k, v in _timed.items() if k != "potential_time_derivative"}),
    (2, {**_static, "potential_time_derivative": (None, None, dazdt_2d)}),
    (2, {**_timed, "potential_time_derivative": (None, None, az_2d)}),
    # every component is called with positional arguments only
    (2, {**_static, "potential": (None, None, az_2d_keyword_only)}),
    (2, {**_timed, "potential_time_derivative": (None, None, dazdt_2d_keyword_only)}),
]

_invalid_instances = [
    (2, ExternalField()),  # the base class maps to no C++ updater
    (1, DipoleExternalField((0.5,), (1.0,), 0.0)),  # a dipole makes no sense in 1D
    (2, DipoleExternalField((0.5, 1.5, 0.5), (0.0, 1.0, 0.0), 0.0)),  # 3 components
    (2, DipoleExternalField((0.5, 1.5), (0.0, 1.0), -1.0)),  # a negative radius
    (2, UserDefinedExternalField(False, (None, None, None), None)),  # all None
    (2, UserDefinedExternalField(False, (None, None, az_3d), None)),  # wrong ndim
    (2, UserDefinedExternalField(True, (None, None, az_2d), None)),  # not time dependent
    (2, UserDefinedExternalField(True, (None, None, az_2d_t), None)),  # no derivative
]


class RecordingPopulator:
    """A dict_populator() stand-in recording what a description writes, and where."""

    def __init__(self):
        self.written = {}

    def add_double(self, path, value):
        self.written[path] = float(value)

    def add_bool(self, path, value):
        self.written[path] = bool(value)

    def add_enum_int(self, path, enum_name, member_name):
        self.written[path] = (enum_name, member_name)

    def add_space_time_function(self, path, fn):
        self.written[path] = fn


@ddt
class TestExternalFieldResolution(unittest.TestCase):
    def test_absent_gives_no_external_field(self):
        self.assertIsInstance(resolve_external_field(2), ZeroExternalField)
        self.assertIsInstance(
            resolve_external_field(2, external_field=None), ZeroExternalField
        )
        self.assertIsInstance(
            resolve_external_field(2, external_field={"type": "zero"}), ZeroExternalField
        )

    def test_dipole_is_resolved(self):
        ef = resolve_external_field(
            2,
            external_field={
                "type": "dipole",
                "position": (0.5, 1.5),
                "moment": (0.0, 2.0),
                "radius": 1,
            },
        )

        assert isinstance(ef, DipoleExternalField)
        self.assertEqual(ef.position, (0.5, 1.5))
        self.assertEqual(ef.moment, (0.0, 2.0))
        self.assertEqual(ef.radius, 1.0)
        self.assertIsInstance(ef.radius, float)

    def test_dipole_accepts_a_zero_radius(self):
        """0 is the point dipole: a valid, explicit choice."""
        ef = resolve_external_field(
            2,
            external_field={
                "type": "dipole",
                "position": (0.5, 1.5),
                "moment": (0.0, 2.0),
                "radius": 0.0,
            },
        )
        self.assertEqual(ef.radius, 0.0)

    def test_user_defined_is_resolved(self):
        ef = resolve_user_defined(2, (None, None, az_2d))

        self.assertIsInstance(ef, UserDefinedExternalField)
        self.assertFalse(ef.is_time_dependent)
        self.assertIsNone(ef.potential_time_derivative)

        # every component is normalized to the space-time signature, the 'None' ones
        # included, so that one C++ SpaceTimeFunction type serves both cases
        self.assertTrue(all(callable(a) for a in ef.potential))
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])
        np.testing.assert_array_equal(ef.potential[2](x, y, 7.0), az_2d(x, y))

    def test_none_components_default_to_zero_shaped_like_the_coordinates(self):
        ef = resolve_user_defined(2, (None, None, az_2d))

        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])

        # an array, not a scalar: the pybind wrapper must not have to broadcast it.
        # the time is accepted and ignored, the potential being static here
        for defaulted in (ef.potential[0], ef.potential[1]):
            np.testing.assert_array_equal(defaulted(x, y, 7.0), np.zeros(x.size))

        np.testing.assert_array_equal(ef.potential[2](x, y, 7.0), x + y)

    def test_user_defined_time_dependence_comes_from_the_signature(self):
        ef = resolve_user_defined(
            2,
            (None, None, az_2d_t),
            potential_time_derivative=(None, None, dazdt_2d),
        )

        dadt = ef.potential_time_derivative
        assert dadt is not None

        self.assertTrue(ef.is_time_dependent)
        self.assertIs(ef.potential[2], az_2d_t)
        self.assertIs(dadt[2], dazdt_2d)

        # the defaulted components take the time as a last argument too
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])
        for defaulted in (ef.potential[0], dadt[0]):
            np.testing.assert_array_equal(defaulted(x, y, 0.5), np.zeros(x.size))

    @data(*_invalid_instances)
    @unpack
    def test_invalid_instances_are_rejected(self, ndim, ef):
        with self.assertRaises(ValueError):
            resolve_external_field(ndim, external_field=ef)

    def test_instances_are_resolved_as_their_dict_declaration(self):
        x, y = np.array([0.0, 1.0]), np.array([2.0, 3.0])

        dipole = resolve_external_field(
            2, external_field=DipoleExternalField([0.5, 1.5], [0.0, 1.0], 0)
        )
        self.assertEqual(dipole, DipoleExternalField((0.5, 1.5), (0.0, 1.0), 0.0))

        ef = resolve_external_field(
            2, external_field=UserDefinedExternalField(False, (None, None, az_2d), None)
        )
        self.assertTrue(ef.resolved)
        np.testing.assert_array_equal(ef.potential[2](x, y, 0.5), az_2d(x, y))
        np.testing.assert_array_equal(ef.potential[0](x, y, 0.5), np.zeros(x.size))

    def test_resolved_instances_are_kept(self):
        ef = resolve_user_defined(2, (None, None, az_2d))
        self.assertIs(resolve_external_field(2, external_field=ef), ef)

    def test_user_defined_is_resolved_in_1d_and_3d(self):
        x = np.array([0.0, 1.0, 2.0])

        # in 1D a_x contributes to no curl term at all, and is expected to be None
        ay = resolve_user_defined(1, (None, ay_1d, None)).potential[1]
        np.testing.assert_array_equal(ay(x, 7.0), ay_1d(x))

        az = resolve_user_defined(3, (None, None, az_3d)).potential[2]
        np.testing.assert_array_equal(az(x, x, x, 7.0), az_3d(x, x, x))

    @data(*_invalid_declarations)
    @unpack
    def test_invalid_declarations_are_rejected(self, ndim, ef):
        with self.assertRaises(ValueError):
            resolve_external_field(ndim, external_field=ef)


class TestExternalFieldPopulateDict(unittest.TestCase):
    def test_zero_writes_only_its_type(self):
        dp = RecordingPopulator()
        ZeroExternalField().populate_dict(dp)

        self.assertEqual(
            dp.written,
            {
                "simulation/external_field/type": (
                    "ExternalFieldUpdaterType",
                    "zero",
                )
            },
        )

    def test_dipole_writes_its_vectors_component_wise(self):
        dp = RecordingPopulator()
        DipoleExternalField((0.5, 1.5), (0.0, 1.0), 0.0).populate_dict(dp)

        self.assertEqual(
            dp.written,
            {
                "simulation/external_field/type": (
                    "ExternalFieldUpdaterType",
                    "dipole",
                ),
                # in 2D both vectors have x and y components only
                "simulation/external_field/position/x": 0.5,
                "simulation/external_field/position/y": 1.5,
                "simulation/external_field/moment/x": 0.0,
                "simulation/external_field/moment/y": 1.0,
                "simulation/external_field/radius": 0.0,
            },
        )

    def test_dipole_writes_its_radius(self):
        dp = RecordingPopulator()
        DipoleExternalField((0.5, 1.5), (0.0, 1.0), 0.25).populate_dict(dp)

        self.assertEqual(dp.written["simulation/external_field/radius"], 0.25)


    def test_user_defined_writes_its_three_components_and_the_time_flag(self):
        dp = RecordingPopulator()
        resolve_user_defined(2, (None, None, az_2d)).populate_dict(dp)

        path = "simulation/external_field"
        self.assertEqual(
            sorted(dp.written),
            sorted(
                [
                    f"{path}/type",
                    f"{path}/is_time_dependent",
                    f"{path}/potential/x",
                    f"{path}/potential/y",
                    f"{path}/potential/z",
                ]
            ),
        )
        self.assertEqual(
            dp.written[f"{path}/type"], ("ExternalFieldUpdaterType", "user-defined")
        )
        self.assertFalse(dp.written[f"{path}/is_time_dependent"])

        # what is written is callable with the space-time signature, whatever the user gave:
        # a_z is theirs, a_x and a_y are the zero defaulters
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])
        np.testing.assert_array_equal(dp.written[f"{path}/potential/z"](x, y, 7.0), x + y)
        np.testing.assert_array_equal(
            dp.written[f"{path}/potential/x"](x, y, 7.0), np.zeros(x.size)
        )

    def test_user_defined_writes_the_derivative_only_when_time_dependent(self):
        dp = RecordingPopulator()
        resolve_user_defined(
            2, (None, None, az_2d_t), potential_time_derivative=(None, None, dazdt_2d)
        ).populate_dict(dp)

        path = "simulation/external_field"
        self.assertTrue(dp.written[f"{path}/is_time_dependent"])
        for axis in "xyz":
            self.assertIn(f"{path}/potential_time_derivative/{axis}", dp.written)

        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])
        np.testing.assert_array_equal(
            dp.written[f"{path}/potential_time_derivative/z"](x, y, 7.0), x + y
        )

        # a static field writes no derivative at all
        static = RecordingPopulator()
        resolve_user_defined(2, (None, None, az_2d)).populate_dict(static)
        self.assertFalse(
            any("potential_time_derivative" in key for key in static.written)
        )


class TestSimulationExternalField(unittest.TestCase):
    def setUp(self):
        global_vars.sim = None

    def tearDown(self):
        global_vars.sim = None

    def test_simulation_defaults_to_no_external_field(self):
        sim = simulation.Simulation(
            time_step=0.001, time_step_nbr=10, cells=(20, 20), dl=(0.1, 0.1)
        )
        self.assertIsInstance(sim.external_field, ZeroExternalField)

    def test_simulation_accepts_a_dipole(self):
        sim = simulation.Simulation(
            time_step=0.001,
            time_step_nbr=10,
            cells=(20, 20),
            dl=(0.1, 0.1),
            external_field={
                "type": "dipole",
                "position": (0.5, 1.5),
                "moment": (0.0, 2.0),
                "radius": 0.0,
            },
        )
        self.assertIsInstance(sim.external_field, DipoleExternalField)

    def test_simulation_accepts_a_user_defined_field(self):
        sim = simulation.Simulation(
            time_step=0.001,
            time_step_nbr=10,
            cells=(20, 20),
            dl=(0.1, 0.1),
            external_field={"type": "user-defined", "potential": (None, None, az_2d)},
        )
        self.assertIsInstance(sim.external_field, UserDefinedExternalField)


if __name__ == "__main__":
    unittest.main()
