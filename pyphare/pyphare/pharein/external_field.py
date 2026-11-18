"""
External magnetic field declarations for pharein.
"""

import inspect
import math
import numbers
from functools import partial

from ..core import phare_utilities
from . import global_vars

_AXES = ("x", "y", "z")
_CPPDICT_PATH = "simulation/external_field"
_POTENTIAL_NAMES = ("a0x", "a0y", "a0z")
_DERIVATIVE_NAMES = ("da0x_dt", "da0y_dt", "da0z_dt")


def _add_vector(component_adder, path, components):
    """Write a vector as the x/y/z sub-dict shape the C++ side reads component by component."""
    for axis, component in zip(_AXES, components):
        component_adder(f"{path}/{axis}", component)


def _check_components(name, values, expected):
    values = phare_utilities.listify(values)

    if len(values) != expected:
        raise ValueError(
            f"Error: external field '{name}' must have {expected} components, got {len(values)}"
        )
    if not all(phare_utilities.is_scalar(v) for v in values):
        raise ValueError(f"Error: external field '{name}' components must be scalars")

    return tuple(float(v) for v in values)


def _check_radius(radius):
    """A finite, non-negative scalar, 0 standing for a point dipole."""
    if isinstance(radius, bool) or not isinstance(radius, numbers.Real):
        raise TypeError(
            f"Error: dipole external field 'radius' must be a number, got {radius!r}"
        )
    if not math.isfinite(radius) or radius < 0.0:
        raise ValueError(
            f"Error: dipole external field 'radius' must be finite and >= 0, got {radius}"
        )
    return float(radius)


def _check_callables(names, values):
    for name, value in zip(names, values):
        if value is not None and not callable(value):
            raise ValueError(
                f"Error: external field '{name}' must be a callable or 'None'"
            )
    return tuple(values)


def _get_signature_size(f):
    return len(inspect.signature(f).parameters)


def _contains_consistent_signature_sizes(list):
    """
    Check that all elements of an iterable containing either
    callables or Nones have the same signature size.
    """
    return len({_get_signature_size(f) for f in list if f is not None}) <= 1


def _check_that_vector_of_callable_is_valid(list, names):
    if phare_utilities.all_none(*list):
        raise ValueError(
            f"Error: a user-defined external field requires one of {', '.join(names)}"
        )
    if not _contains_consistent_signature_sizes(list):
        raise ValueError(
            f"Error: {', '.join(names)} do not have a consistent signature."
        )


def _first_callable(list):
    """First non-None element; the caller must have rejected the all-None case."""
    return next(f for f in list if f is not None)


def _check_signature_size(list, names, ndim, is_time_dependent):
    """
    Every provided component takes the 'ndim' coordinates, and the time as a last
    argument when the field is time dependent.
    """
    expected = ndim + 1 if is_time_dependent else ndim
    arguments = "coordinates and time" if is_time_dependent else "coordinates"

    for name, f in zip(names, list):
        if f is None:
            continue
        size = _get_signature_size(f)
        if size != expected:
            raise ValueError(
                f"Error: external field '{name}' takes {size} "
                f"argument(s), expected {expected} ({arguments})"
            )
        try:
            inspect.signature(f).bind(*([None] * expected))
        except TypeError as exc:
            raise ValueError(
                f"Error: external field '{name}' must accept "
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


class ExternalMagneticField:
    """
    Base class of the external magnetic field declarations, not declarable itself:
    declare a :class:`~pyphare.pharein.DipoleExternalField` or a
    :class:`~pyphare.pharein.UserDefinedExternalField`.

    'type' names the C++ ExternalFieldUpdaterType member a declaration maps onto.
    """

    type = None

    def __init__(self):
        raise TypeError(
            "Error: ExternalMagneticField is not declarable, "
            "declare a DipoleExternalField or a UserDefinedExternalField"
        )

    @staticmethod
    def _declaring_simulation():
        sim = global_vars.sim
        if sim is None:
            raise RuntimeError("A simulation must be declared before an external field")
        if sim.model is not None:
            raise RuntimeError("An external field must be declared before the model")
        if sim.external_field is not None:
            raise RuntimeError("An external field is already declared")
        return sim

    def populate_dict(self, dp, ndim):
        """Mirror the declaration (type + per-type params) on the C++ side.

        `dp` is a DictPopulator (see pharein.initialize): an object exposing
        add_string/add_double/add_int/add_enum_int/... - passed in rather than imported, to
        avoid a circular import between this module and pharein.initialize.

        :meta private:
        """
        dp.add_enum_int(f"{_CPPDICT_PATH}/type", "ExternalFieldUpdaterType", self.type)


class DipoleExternalField(ExternalMagneticField):
    """
    A static magnetic dipole, prescribed as the external (background) magnetic
    field B0 of a Hybrid or MHD simulation, in 2D or 3D. Declare it at most once,
    after the :class:`~pyphare.pharein.Simulation` and before the model.

    **Usage example:**

    .. code-block:: python

        import pyphare.pharein as ph

        ph.Simulation(...)
        ph.DipoleExternalField(position=(10.0, 10.0), moment=(0.0, 1.0), radius=1.0)
        ph.MaxwellianFluidModel(...)

    **Parameters**:

        * **position** (*tuple*), the dipole center, one component per dimension.
        * **moment** (*tuple*), the dipole moment, one component per dimension, not zero.
        * **radius** (*float*), finite and >= 0: inside this radius around the center the
          field is uniform, which removes the singularity at the center. 0 is the point
          dipole.
    """

    type = "dipole"

    def __init__(self, *, position, moment, radius):
        sim = self._declaring_simulation()
        ndim = sim.ndim

        if ndim == 1:
            raise ValueError("Error: a dipole external field makes no sense in 1D")

        self.position = _check_components("position", position, ndim)
        self.moment = _check_components("moment", moment, ndim)

        if all(m == 0.0 for m in self.moment):
            raise ValueError("Error: dipole external field 'moment' cannot be zero")

        self.radius = _check_radius(radius)

        sim.set_external_field(self)

    def populate_dict(self, dp, ndim):
        super().populate_dict(dp, ndim)
        _add_vector(dp.add_double, f"{_CPPDICT_PATH}/position", self.position)
        _add_vector(dp.add_double, f"{_CPPDICT_PATH}/moment", self.moment)
        dp.add_double(f"{_CPPDICT_PATH}/radius", self.radius)


class UserDefinedExternalField(ExternalMagneticField):
    """
    An external (background) magnetic field B0 = curl(A0) given by the components of
    its vector potential A0, for a Hybrid or MHD simulation. Declare it at most once,
    after the :class:`~pyphare.pharein.Simulation` and before the model.

    **Usage example:**

    .. code-block:: python

        import pyphare.pharein as ph

        ph.Simulation(...)

        def az(x, y):
            return x * y

        ph.UserDefinedExternalField(a0z=az)
        ph.MaxwellianFluidModel(...)

    **Parameters**:

        * **a0x**, **a0y**, **a0z** (*function*), default None, the components of A0,
          functions of the spatial coordinates (e.g. ``a0z(x, y)`` in 2D). A ``None``
          component is zero; at least one is given. In 1D, ``a0x`` does not contribute to
          the field and must stay ``None``. The potential is time dependent when its
          components also take the time as a last argument (e.g. ``a0z(x, y, t)``).
        * **da0x_dt**, **da0y_dt**, **da0z_dt** (*function*), default None, the time
          derivative of A0, functions of the spatial coordinates and the time. Required
          (at least one) when the potential is time dependent, forbidden otherwise.
    """

    type = "user-defined"

    def __init__(
        self,
        *,
        a0x=None,
        a0y=None,
        a0z=None,
        da0x_dt=None,
        da0y_dt=None,
        da0z_dt=None,
    ):
        sim = self._declaring_simulation()
        ndim = sim.ndim

        potential = _check_callables(_POTENTIAL_NAMES, (a0x, a0y, a0z))
        _check_that_vector_of_callable_is_valid(potential, _POTENTIAL_NAMES)

        size = _get_signature_size(_first_callable(potential))
        if size not in (ndim, ndim + 1):
            raise ValueError(
                f"Error: external field {', '.join(_POTENTIAL_NAMES)} must take {ndim} "
                f"coordinate(s), optionally followed by the time, "
                f"got {size} argument(s)"
            )
        is_time_dependent = size == ndim + 1

        _check_signature_size(potential, _POTENTIAL_NAMES, ndim, is_time_dependent)

        if ndim == 1 and a0x is not None:
            raise ValueError(
                "Error: in 1D, a0x cannot contribute to the external magnetic field. "
                "Expected as None."
            )

        derivative = (da0x_dt, da0y_dt, da0z_dt)
        has_derivative = not phare_utilities.all_none(*derivative)
        if is_time_dependent and not has_derivative:
            raise ValueError(
                "Error: a time-dependent potential requires one of "
                f"{', '.join(_DERIVATIVE_NAMES)}"
            )
        if has_derivative and not is_time_dependent:
            raise ValueError(
                f"Error: {', '.join(_DERIVATIVE_NAMES)} given but the potential "
                "does not depend on time"
            )

        potential_time_derivative = None
        if is_time_dependent:
            derivative = _check_callables(_DERIVATIVE_NAMES, derivative)
            _check_signature_size(derivative, _DERIVATIVE_NAMES, ndim, True)
            potential_time_derivative = _normalized_components(derivative, ndim, True)

        self.is_time_dependent = is_time_dependent
        self.potential = _normalized_components(potential, ndim, is_time_dependent)
        self.potential_time_derivative = potential_time_derivative

        sim.set_external_field(self)

    def populate_dict(self, dp, ndim):
        super().populate_dict(dp, ndim)
        dp.add_bool(f"{_CPPDICT_PATH}/is_time_dependent", self.is_time_dependent)
        add_space_time_function = partial(dp.add_space_time_function, ndim)
        _add_vector(
            add_space_time_function, f"{_CPPDICT_PATH}/potential", self.potential
        )
        if self.is_time_dependent:
            _add_vector(
                add_space_time_function,
                f"{_CPPDICT_PATH}/potential_time_derivative",
                self.potential_time_derivative,
            )
