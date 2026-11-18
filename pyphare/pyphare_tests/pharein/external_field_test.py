import unittest

import numpy as np
from ddt import data, ddt, unpack

from pyphare.pharein import global_vars, simulation
from pyphare.pharein.diagnostics import ElectromagDiagnostics
from pyphare.pharein.external_field import (
    DipoleExternalField,
    ExternalMagneticField,
    UserDefinedExternalField,
)
from pyphare.pharein.uniform_model import UniformModel


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


def declare_simulation(ndim):
    global_vars.sim = None
    return simulation.Simulation(
        time_step=0.001, time_step_nbr=10, cells=(20,) * ndim, dl=(0.1,) * ndim
    )


_dipole = dict(position=(0.5, 1.5), moment=(0.0, 2.0), radius=0.0)
_dipole_1d = dict(position=(0.5,), moment=(1.0,), radius=0.0)

D, U = DipoleExternalField, UserDefinedExternalField

_invalid_declarations = [
    (1, D, _dipole_1d),  # a dipole makes no sense in 1D
    (3, D, _dipole),  # both vectors must have ndim components
    (2, D, {**_dipole, "moment": (0.0, 1.0, 2.0)}),  # three components in 2D
    (2, D, {**_dipole, "moment": (0.0, 0.0)}),  # a zero moment is no dipole
    (2, D, {**_dipole, "radius": -1.0}),  # a negative radius
    (2, D, {**_dipole, "radius": float("nan")}),  # not a number
    (2, D, {**_dipole, "radius": float("inf")}),  # an infinite radius
    #
    (2, U, {}),  # all components None
    (2, U, dict(a0z=1.0)),  # not a callable
    (2, U, dict(a0x=az_2d, a0z=az_2d_t)),  # mixed signatures
    (3, U, dict(a0z=az_2d)),  # signature does not match ndim
    (1, U, dict(a0x=ay_1d)),  # a0x is unused in 1D
    # the time dependence of the potential and of its derivative must agree
    (2, U, dict(a0z=az_2d_t)),
    (2, U, dict(a0z=az_2d, da0z_dt=dazdt_2d)),
    (2, U, dict(a0z=az_2d_t, da0z_dt=az_2d)),
    (2, U, dict(a0z=az_2d_t, da0z_dt=1.0)),
    # every component is called with positional arguments only
    (2, U, dict(a0z=az_2d_keyword_only)),
    (2, U, dict(a0z=az_2d_t, da0z_dt=dazdt_2d_keyword_only)),
]

_mistyped_declarations = [
    (2, D, {**_dipole, "radius": (1.0, 1.0)}),  # not a scalar
    (2, D, {**_dipole, "radius": "1"}),  # not a number
    (2, D, {k: v for k, v in _dipole.items() if k != "moment"}),  # missing argument
    (2, D, {**_dipole, "value": 3}),  # unknown argument
    (2, U, dict(a0z=az_2d, moment=(0.0, 1.0))),  # unknown argument
    (2, U, dict(potential=(None, None, az_2d))),
]


class RecordingPopulator:
    """A DictPopulator stand-in recording what a description writes, and where."""

    def __init__(self):
        self.written = {}

    def add_double(self, path, value):
        self.written[path] = float(value)

    def add_bool(self, path, value):
        self.written[path] = bool(value)

    def add_enum_int(self, path, enum_name, member_name):
        self.written[path] = (enum_name, member_name)

    def add_space_time_function(self, ndim, path, fn):
        self.written[path] = fn


class SimulationResetting(unittest.TestCase):
    def setUp(self):
        global_vars.sim = None

    def tearDown(self):
        global_vars.sim = None


@ddt
class TestExternalFieldDeclaration(SimulationResetting):
    def test_dipole_is_registered(self):
        sim = declare_simulation(2)
        ef = DipoleExternalField(position=(0.5, 1.5), moment=(0.0, 2.0), radius=1)

        self.assertIs(sim.external_field, ef)
        self.assertEqual(ef.position, (0.5, 1.5))
        self.assertEqual(ef.moment, (0.0, 2.0))
        self.assertEqual(ef.radius, 1.0)
        self.assertIsInstance(ef.radius, float)

    def test_dipole_accepts_a_zero_radius(self):
        """0 is the point dipole: a valid, explicit choice."""
        sim = declare_simulation(2)
        ef = DipoleExternalField(**_dipole)
        self.assertIs(sim.external_field, ef)
        self.assertEqual(ef.radius, 0.0)

    def test_dipole_accepts_lists_and_ints(self):
        sim = declare_simulation(2)
        ef = DipoleExternalField(position=[0.5, 1.5], moment=[0, 2], radius=0)
        self.assertIs(sim.external_field, ef)
        self.assertEqual(ef.position, (0.5, 1.5))
        self.assertEqual(ef.moment, (0.0, 2.0))
        self.assertTrue(all(isinstance(m, float) for m in ef.moment))

    def test_user_defined_is_resolved(self):
        sim = declare_simulation(2)
        ef = UserDefinedExternalField(a0z=az_2d)

        self.assertIs(sim.external_field, ef)
        self.assertFalse(ef.is_time_dependent)
        self.assertIsNone(ef.potential_time_derivative)

        # every component is normalized to the space-time signature, the 'None' ones
        # included, so that one C++ SpaceTimeFunction type serves both cases
        self.assertTrue(all(callable(a) for a in ef.potential))
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])
        np.testing.assert_array_equal(ef.potential[2](x, y, 7.0), az_2d(x, y))

    def test_none_components_default_to_zero_shaped_like_the_coordinates(self):
        declare_simulation(2)
        ef = UserDefinedExternalField(a0z=az_2d)

        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])

        # an array, not a scalar: the pybind wrapper must not have to broadcast it.
        # the time is accepted and ignored, the potential being static here
        for defaulted in (ef.potential[0], ef.potential[1]):
            np.testing.assert_array_equal(defaulted(x, y, 7.0), np.zeros(x.size))

        np.testing.assert_array_equal(ef.potential[2](x, y, 7.0), x + y)

    def test_user_defined_time_dependence_comes_from_the_signature(self):
        declare_simulation(2)
        ef = UserDefinedExternalField(a0z=az_2d_t, da0z_dt=dazdt_2d)

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

    def test_user_defined_is_resolved_in_1d_and_3d(self):
        x = np.array([0.0, 1.0, 2.0])

        # in 1D a0x contributes to no curl term at all, and is expected to be None
        declare_simulation(1)
        ay = UserDefinedExternalField(a0y=ay_1d).potential[1]
        np.testing.assert_array_equal(ay(x, 7.0), ay_1d(x))

        declare_simulation(3)
        az = UserDefinedExternalField(a0z=az_3d).potential[2]
        np.testing.assert_array_equal(az(x, x, x, 7.0), az_3d(x, x, x))

    @data(*_invalid_declarations)
    @unpack
    def test_invalid_declarations_are_rejected(self, ndim, cls, kwargs):
        sim = declare_simulation(ndim)
        with self.assertRaises(ValueError):
            cls(**kwargs)
        self.assertIsNone(sim.external_field)

    @data(*_mistyped_declarations)
    @unpack
    def test_mistyped_declarations_are_rejected(self, ndim, cls, kwargs):
        sim = declare_simulation(ndim)
        with self.assertRaises(TypeError):
            cls(**kwargs)
        self.assertIsNone(sim.external_field)

    def test_error_names_the_offending_component(self):
        declare_simulation(2)
        with self.assertRaisesRegex(ValueError, "da0z_dt"):
            UserDefinedExternalField(a0z=az_2d_t, da0z_dt=az_2d)


class TestExternalFieldOrdering(SimulationResetting):
    def test_declaring_without_a_simulation_raises(self):
        with self.assertRaises(RuntimeError):
            DipoleExternalField(**_dipole)
        with self.assertRaises(RuntimeError):
            UserDefinedExternalField(a0z=az_2d)

    def test_declaring_after_the_model_raises(self):
        sim = declare_simulation(2)
        UniformModel()
        with self.assertRaises(RuntimeError):
            UserDefinedExternalField(a0z=az_2d)
        self.assertIsNone(sim.external_field)

    def test_declaring_twice_raises_and_keeps_the_first(self):
        sim = declare_simulation(2)
        first = UserDefinedExternalField(a0z=az_2d)
        with self.assertRaises(RuntimeError):
            DipoleExternalField(**_dipole)
        self.assertIs(sim.external_field, first)

    def test_the_base_class_is_not_declarable(self):
        sim = declare_simulation(2)
        with self.assertRaises(TypeError):
            ExternalMagneticField()
        self.assertIsNone(sim.external_field)

    def test_simulation_rejects_the_external_field_keyword(self):
        with self.assertRaises(ValueError):
            simulation.Simulation(
                time_step=0.001,
                time_step_nbr=10,
                cells=(20, 20),
                dl=(0.1, 0.1),
                external_field={
                    "type": "user-defined",
                    "potential": (None, None, az_2d),
                },
            )


class TestExternalFieldPopulateDict(SimulationResetting):
    def test_dipole_writes_its_vectors_component_wise(self):
        declare_simulation(2)
        dp = RecordingPopulator()
        DipoleExternalField(
            position=(0.5, 1.5), moment=(0.0, 1.0), radius=0.0
        ).populate_dict(dp, 2)

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
        declare_simulation(2)
        dp = RecordingPopulator()
        DipoleExternalField(
            position=(0.5, 1.5), moment=(0.0, 1.0), radius=0.25
        ).populate_dict(dp, 2)

        self.assertEqual(dp.written["simulation/external_field/radius"], 0.25)

    def test_user_defined_writes_its_three_components_and_the_time_flag(self):
        declare_simulation(2)
        dp = RecordingPopulator()
        UserDefinedExternalField(a0z=az_2d).populate_dict(dp, 2)

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
        # a0z is theirs, a0x and a0y are the zero defaulters
        x = np.array([0.0, 1.0, 2.0])
        y = np.array([3.0, 4.0, 5.0])
        np.testing.assert_array_equal(
            dp.written[f"{path}/potential/z"](x, y, 7.0), x + y
        )
        np.testing.assert_array_equal(
            dp.written[f"{path}/potential/x"](x, y, 7.0), np.zeros(x.size)
        )

    def test_user_defined_writes_the_derivative_only_when_time_dependent(self):
        declare_simulation(2)
        dp = RecordingPopulator()
        UserDefinedExternalField(a0z=az_2d_t, da0z_dt=dazdt_2d).populate_dict(dp, 2)

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
        declare_simulation(2)
        static = RecordingPopulator()
        UserDefinedExternalField(a0z=az_2d).populate_dict(static, 2)
        self.assertFalse(
            any("potential_time_derivative" in key for key in static.written)
        )


class TestSimulationExternalField(SimulationResetting):
    def test_simulation_defaults_to_no_external_field(self):
        sim = declare_simulation(2)
        self.assertIsNone(sim.external_field)

    def test_b0_diagnostics_require_an_external_field(self):
        declare_simulation(2)
        with self.assertRaises(RuntimeError):
            ElectromagDiagnostics(quantity="B0", write_timestamps=[0.0])

    def test_b0_diagnostics_with_an_external_field(self):
        declare_simulation(2)
        UserDefinedExternalField(a0z=az_2d)
        diag = ElectromagDiagnostics(quantity="B0", write_timestamps=[0.0])
        self.assertEqual(diag.quantity, "/EM_B0")


if __name__ == "__main__":
    unittest.main()
