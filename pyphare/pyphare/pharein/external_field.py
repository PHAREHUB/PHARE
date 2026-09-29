"""
External magnetic field resolution and validation for pharein.Simulation.
"""

import inspect
import math
import numbers
from abc import ABC
from dataclasses import dataclass, field

from ..core import phare_utilities

_AXES = ("x", "y", "z")
_CPPDICT_PATH = "simulation/external_field"


def _add_vector(component_adder, path, components):
    """Write a vector as the x/y/z sub-dict shape the C++ side reads component by component."""
    for axis, component in zip(_AXES, components):
        component_adder(f"{path}/{axis}", component)


def _check_components(name, dict, expected):
    values = phare_utilities.listify(dict[name])

    if len(values) != expected:
        raise ValueError(
            f"Error: external_field '{name}' must have {expected} components, got {len(values)}"
        )
    if not all(phare_utilities.is_scalar(v) for v in values):
        raise ValueError(f"Error: external_field '{name}' components must be scalars")

    return tuple(float(v) for v in values)


def _check_radius(radius):
    """A finite, non-negative scalar, 0 standing for a point dipole."""
    if isinstance(radius, bool) or not isinstance(radius, numbers.Real):
        raise ValueError(
            f"Error: dipole external_field 'radius' must be a number, got {radius!r}"
        )
    if not math.isfinite(radius) or radius < 0.0:
        raise ValueError(
            f"Error: dipole external_field 'radius' must be finite and >= 0, got {radius}"
        )
    return float(radius)


def _check_callables(name, dict):
    """
    Read a vector of callables from the external_field dict.
    """
    values = phare_utilities.listify(dict[name])

    if len(values) != len(_AXES):
        raise ValueError(
            f"Error: external_field '{name}' must have {len(_AXES)} components "
            f"(one per axis, whatever the dimensionality), got {len(values)}"
        )
    if not all(v is None or callable(v) for v in values):
        raise ValueError(
            f"Error: external_field '{name}' components must be callables or 'None'"
        )

    return tuple(values)


@dataclass
class ExternalField(ABC):
    """
    Base class of the external field descriptions.

    'type' names the C++ ExternalFieldUpdaterType member this description maps onto.
    """

    type = None

    def populate_dict(self, dp):
        """Mirror the public `external_field` dict shape (type + per-type params) on the C++ side.

        `dp` is a dict_populator() (see pharein.initialize.general): an object exposing
        add_string/add_double/add_int/add_enum_int/... - passed in rather than imported, to
        avoid a circular import between this module and pharein.initialize.general.
        """
        dp.add_enum_int(
            "simulation/external_field/type", "ExternalFieldUpdaterType", self.type
        )


@dataclass
class ZeroExternalField(ExternalField):
    """Zero external field: B0 and its time derivative stay zero."""

    type = "zero"


@dataclass
class DipoleExternalField(ExternalField):
    """
    A static magnetic dipole of moment 'moment' placed at 'position',
    with uniform field below inside 'radius'.
    """

    type = "dipole"

    position: tuple
    moment: tuple
    radius: float

    def populate_dict(self, dp):
        super().populate_dict(dp)
        _add_vector(dp.add_double, "simulation/external_field/position", self.position)
        _add_vector(dp.add_double, "simulation/external_field/moment", self.moment)
        dp.add_double("simulation/external_field/radius", self.radius)


@dataclass
class UserDefinedExternalField(ExternalField):
    """
    A external magnetic field prescribed via pharein
    """

    type = "user-defined"

    is_time_dependent: bool
    potential: tuple
    potential_time_derivative: tuple | None
    resolved: bool = field(default=False, kw_only=True, repr=False, compare=False)

    def populate_dict(self, dp):
        super().populate_dict(dp)
        dp.add_bool(f"{_CPPDICT_PATH}/is_time_dependent", self.is_time_dependent)
        _add_vector(
            dp.add_space_time_function, f"{_CPPDICT_PATH}/potential", self.potential
        )
        if self.is_time_dependent:
            _add_vector(
                dp.add_space_time_function,
                f"{_CPPDICT_PATH}/potential_time_derivative",
                self.potential_time_derivative,
            )


# ------------------------------------------------------------------------------
def _check_keys(keys, allowed, type):
    extra = set(keys) - allowed
    if extra:
        raise ValueError(
            f"Error: invalid external_field keys for type '{type}': {sorted(extra)}, "
            f"allowed {sorted(allowed)}"
        )


def _get_signature_size(f):
    return len(inspect.signature(f).parameters)


def _is_signature_time_dependent(f, ndim):
    return _get_signature_size(f) == ndim + 1


def _contains_only_none(list):
    return all(f is None for f in list)


def _contains_consistent_signature_sizes(list):
    """
    Check that all elements of an iterable containing either
    callables or Nones have the same signature size.
    """
    return len(set([_get_signature_size(f) for f in list if f is not None])) <= 1


def _check_that_vector_of_callable_is_valid(list, name):
    if _contains_only_none(list):
        raise ValueError(f"Error: user-defined '{name}' contains only 'None'")
    if not _contains_consistent_signature_sizes(list):
        raise ValueError(
            f"Error: elements in '{name}' have not a consistent signature."
        )


def _first_callable(list):
    """First non-None element; the caller must have rejected the all-None case."""
    return next(f for f in list if f is not None)


def _check_signature_size(list, name, ndim, is_time_dependent):
    """
    Every provided component takes the 'ndim' coordinates, and the time as a last
    argument when the field is time dependent.
    """
    expected = ndim + 1 if is_time_dependent else ndim
    arguments = "coordinates and time" if is_time_dependent else "coordinates"

    for axis, f in zip(_AXES, list):
        if f is None:
            continue
        size = _get_signature_size(f)
        if size != expected:
            raise ValueError(
                f"Error: external_field '{name}' component 'a{axis}' takes {size} "
                f"argument(s), expected {expected} ({arguments})"
            )
        try:
            inspect.signature(f).bind(*([None] * expected))
        except TypeError as exc:
            raise ValueError(
                f"Error: external_field '{name}' component 'a{axis}' must accept "
                f"{expected} positional argument(s) ({arguments})"
            ) from exc


def _zero_defaulter(ndim):
    if ndim == 1:
        return lambda x, t: x * 0.0
    if ndim == 2:
        return lambda x, y, t: x * 0.0
    if ndim == 3:
        return lambda x, y, z, t: x * 0.0


def _normalized_components(list, ndim, is_time_dependent):
    """
    Do two things:
    - None component -> zero function
    - spatial-only functions -> space-time functions (ex: f(x) -> f(x, t))
    """
    zero = _zero_defaulter(ndim)

    def as_space_time(f):
        if f is None:
            return zero
        return f if is_time_dependent else (lambda *args: f(*args[:-1]))

    return tuple(as_space_time(f) for f in list)


def _resolve_dict_user_defined_external_field(external_field, *, ndim):
    _check_keys(
        external_field,
        {"type", "potential", "potential_time_derivative"},
        "user-defined",
    )
    if "potential" not in external_field:
        raise ValueError("Error: a user-defined external field requires 'potential'")
    potential = _check_callables("potential", external_field)
    _check_that_vector_of_callable_is_valid(potential, "potential")

    # the signature size of the components the user did provide tells whether the
    # potential depends on time; potential[0] may be None, hence _first_callable
    size = _get_signature_size(_first_callable(potential))
    if size not in (ndim, ndim + 1):
        raise ValueError(
            f"Error: external_field 'potential' components must take {ndim} "
            f"coordinate(s), optionally followed by the time, "
            f"got {size} argument(s)"
        )
    is_time_dependent = _is_signature_time_dependent(_first_callable(potential), ndim)

    _check_signature_size(potential, "potential", ndim, is_time_dependent)

    if ndim == 1 and potential[0] is not None:
        raise ValueError(
            "Error: in 1D, first component of the potential cannot"
            " contribute to the external magnetic field. Expected as None."
        )

    potential = _normalized_components(potential, ndim, is_time_dependent)

    has_derivative = "potential_time_derivative" in external_field
    if is_time_dependent and not has_derivative:
        raise ValueError(
            "Error: a time-dependent 'potential' requires "
            "'potential_time_derivative'"
        )
    if has_derivative and not is_time_dependent:
        raise ValueError(
            "Error: 'potential_time_derivative' was given but 'potential' does "
            "not depend on time"
        )

    potential_time_derivative = None
    if is_time_dependent:
        potential_time_derivative = _check_callables(
            "potential_time_derivative", external_field
        )
        _check_that_vector_of_callable_is_valid(
            potential_time_derivative, "potential_time_derivative"
        )
        _check_signature_size(
            potential_time_derivative,
            "potential_time_derivative",
            ndim,
            is_time_dependent,
        )
        potential_time_derivative = _normalized_components(
            potential_time_derivative, ndim, is_time_dependent
        )

    return UserDefinedExternalField(
        is_time_dependent, potential, potential_time_derivative, resolved=True
    )


def _resolve_dict_external_field(external_field, *, ndim):
    valid_types = ("zero", "dipole", "user-defined")
    type_ = external_field.get("type")
    if type_ not in valid_types:
        raise ValueError(
            f"Error: external_field dict requires 'type' in {valid_types}, got {type_!r}"
        )
    if type_ == "zero":
        _check_keys(external_field, {"type"}, type_)
        return ZeroExternalField()
    elif type_ == "dipole":
        _check_keys(external_field, {"type", "position", "moment", "radius"}, type_)
        if ndim == 1:
            raise ValueError("Error: a dipole external_field makes no sense in 1D")
        for key in ("position", "moment", "radius"):
            if key not in external_field:
                raise ValueError(f"Error: dipole external_field requires '{key}'")

        position = _check_components("position", external_field, ndim)
        moment = _check_components("moment", external_field, ndim)

        if all(m == 0.0 for m in moment):
            raise ValueError("Error: dipole external_field 'moment' cannot be zero")

        radius = _check_radius(external_field["radius"])

        return DipoleExternalField(position, moment, radius)
    elif type_ == "user-defined":
        return _resolve_dict_user_defined_external_field(external_field, ndim=ndim)


def _resolve_instance_external_field(external_field, *, ndim):
    if type(external_field) is ZeroExternalField:
        return external_field

    if type(external_field) is DipoleExternalField:
        return _resolve_dict_external_field(
            {
                "type": "dipole",
                "position": external_field.position,
                "moment": external_field.moment,
                "radius": external_field.radius,
            },
            ndim=ndim,
        )

    if type(external_field) is UserDefinedExternalField:
        if external_field.resolved:
            return external_field
        declaration = {"type": "user-defined", "potential": external_field.potential}
        if external_field.potential_time_derivative is not None:
            declaration["potential_time_derivative"] = (
                external_field.potential_time_derivative
            )
        resolved = _resolve_dict_user_defined_external_field(declaration, ndim=ndim)
        if resolved.is_time_dependent != external_field.is_time_dependent:
            raise ValueError(
                "Error: user-defined external_field 'is_time_dependent' is "
                f"{external_field.is_time_dependent}, but the signature of its "
                "'potential' components says otherwise"
            )
        return resolved

    raise ValueError(
        f"Error: unsupported external_field type {type(external_field).__name__}"
    )


def resolve_external_field(ndim, **kwargs):
    """
    Resolve the public 'external_field' Simulation option into a validated ExternalField.
    """
    external_field = kwargs.get("external_field")

    if external_field is None:
        return ZeroExternalField()

    if isinstance(external_field, ExternalField):
        return _resolve_instance_external_field(external_field, ndim=ndim)

    if not isinstance(external_field, dict):
        raise ValueError(
            "Error: external_field must be a dict or an ExternalField, got "
            f"{type(external_field).__name__}"
        )

    return _resolve_dict_external_field(external_field, ndim=ndim)
