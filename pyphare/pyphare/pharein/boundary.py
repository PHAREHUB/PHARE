"""
Boundary-condition resolution and validation for pharein.Simulation.

'domain_boundaries' (the public Simulation constructor option) is a per-location dict:
  {"xlower": {"type": "open"}, "xupper": {"type": "super-magnetofast-inflow", "density": ..., ...}}

Each location dict holds the boundary type and, at the same level, the parameters of that type.

A direction is periodic unless one of its locations appears in 'domain_boundaries', in which case it is
physical and both of its locations must be given.

resolve_boundaries(ndim, **kwargs) is the single entry point used by Simulation's checker()
pipeline; it returns the per-direction periodicity flags, stored as Simulation.periodicities, and
a dict[location] -> Boundary stored as Simulation.domain_boundaries.
"""

import math
import numbers
from abc import ABC
from dataclasses import InitVar, dataclass, fields
from inspect import signature

from ..core import phare_utilities


@dataclass
class Boundary(ABC):
    """Base class holding the serialisation entry point common to every BC type."""

    type = None

    def populate_dict(self, dp, bc_path, ndim):
        dp.add_enum_int(
            f"{bc_path}/type", "BoundaryType", boundary_type_enum_member(self.type)
        )


def boundary_type_enum_member(boundary_type):
    return boundary_type.replace("-", "_")


@dataclass
class NoneBoundary(Boundary):
    type = "none"


@dataclass
class OpenBoundary(Boundary):
    type = "open"


@dataclass
class ReflectiveBoundary(Boundary):
    type = "reflective"


_BOUNDARY_NORMAL_INDEX = {"x": 0, "y": 1, "z": 2}


def _normalize_inflow_scalar(location, key, val, ndim, positive=False):
    """A prescribable inflow scalar: a finite float."""
    if (
        not isinstance(val, numbers.Real)
        or not math.isfinite(val)
        or (positive and not val > 0)
    ):
        raise ValueError(
            f"'{key}' at inflow boundary '{location}' must be a finite "
            f"{'positive ' if positive else ''}scalar, got {val!r}"
        )
    return float(val)


def _normalize_inflow_vector(location, key, vec, ndim):
    """A prescribable inflow 3-vector of finite floats."""
    try:
        comps = list(vec)
    except TypeError:
        raise TypeError(
            f"'{key}' at inflow boundary '{location}' must be a 3-vector of floats, "
            f"got {vec!r}"
        )
    if len(comps) != 3:
        raise ValueError(
            f"'{key}' at inflow boundary '{location}' must be a 3-vector, "
            f"got a {len(comps)}-element sequence"
        )
    return tuple(
        _normalize_inflow_scalar(location, f"{key}[{i}]", c, ndim)
        for i, c in enumerate(comps)
    )


def _normalize_inflow_velocity(location, velocity, ndim):
    """Return velocity as a (vx, vy, vz) tuple of floats.

    A scalar constant is interpreted as the inward-normal speed: it is stored with a
    positive sign for lower boundaries (flow enters in the +direction) and a
    negative sign for upper boundaries (flow enters in the -direction).
    """
    if isinstance(velocity, numbers.Real):
        normal_idx = _BOUNDARY_NORMAL_INDEX[location[0]]
        side = location[1:]
        sign = 1.0 if side == "lower" else -1.0
        v = [0.0, 0.0, 0.0]
        v[normal_idx] = sign * float(velocity)
        return tuple(v)
    return _normalize_inflow_vector(location, "velocity", velocity, ndim)


@dataclass
class SuperMagnetofastInflowBoundary(Boundary):
    type = "super-magnetofast-inflow"
    density: object
    pressure: object
    velocity: object
    B: object
    location: InitVar[str]
    ndim: InitVar[int]

    def __post_init__(self, location, ndim):
        self.density = _normalize_inflow_scalar(
            location, "density", self.density, ndim, positive=True
        )
        self.pressure = _normalize_inflow_scalar(
            location, "pressure", self.pressure, ndim, positive=True
        )
        self.velocity = _normalize_inflow_velocity(location, self.velocity, ndim)
        self.B = _normalize_inflow_vector(location, "B", self.B, ndim)

    def populate_dict(self, dp, bc_path, ndim):
        super().populate_dict(dp, bc_path, ndim)
        dp.add_double(f"{bc_path}/density", self.density)
        dp.add_double(f"{bc_path}/pressure", self.pressure)
        _add_inflow_vector(dp, bc_path, "velocity", self.velocity)
        _add_inflow_vector(dp, bc_path, "B", self.B)


# ------------------------------------------------------------------------------

# this dict fills itself automatically once a Boundary's subclass
# has been declared with it type attribute
_type_to_class = {
    el.type: el
    for el in globals().values()
    if isinstance(el, type) and issubclass(el, Boundary) and el is not Boundary
}


def resolve_boundaries(ndim, **kwargs):
    all_directions = ["x", "y", "z"][:ndim]
    sides = ("lower", "upper")

    all_boundary_locations = [
        f"{direction}{side}" for side in sides for direction in all_directions
    ]
    raw = kwargs.get("domain_boundaries", {})
    if not isinstance(raw, dict):
        raise TypeError("A dict should be passed to argument 'domain_boundaries'")

    for location in raw:
        if location not in all_boundary_locations:
            raise ValueError(
                f"Wrong boundary name {location}: should belong to {all_boundary_locations}"
            )

    periodicities = []
    for direction in all_directions:
        given = [f"{direction}{side}" in raw for side in sides]
        if any(given) and not all(given):
            missing = [f"{direction}{side}" for side, g in zip(sides, given) if not g]
            raise ValueError(
                f"Direction '{direction}' is physical since one of its boundaries is given in "
                f"'domain_boundaries', but {missing[0]} is missing"
            )
        periodicities.append(not all(given))

    model_options = phare_utilities.listify(kwargs.get("model_options", "HybridModel"))
    if raw and "MHDModel" not in model_options:
        raise ValueError(
            "non-periodic domain boundaries are only supported by the MHDModel; "
            f"got model_options={model_options}"
        )

    user_types = tuple(t for t in _type_to_class if t != "none")
    resolved = {}
    for location in all_boundary_locations:
        if location not in raw:
            resolved[location] = NoneBoundary()
            continue

        bc = raw[location]
        if not isinstance(bc, dict):
            raise TypeError(
                f"A dict should be passed to the boundary {location} for specifying a "
                f"boundary condition"
            )
        if "type" not in bc:
            raise KeyError(
                f"No key 'type' found in the domain_boundaries dict passed to {location}"
            )
        boundary_type = bc["type"]
        if boundary_type not in user_types:
            raise ValueError(
                f"Boundary type {boundary_type} is not valid: it should belong to {user_types}"
            )

        ctor = _type_to_class[boundary_type]
        params = {key: val for key, val in bc.items() if key != "type"}
        expected = [field.name for field in fields(ctor)]

        unexpected = [key for key in params if key not in expected]
        if unexpected:
            raise ValueError(
                f"Boundary type '{boundary_type}' at '{location}' does not take "
                f"{unexpected}; expected keys: {expected}"
            )
        for key in expected:
            if key not in params:
                raise KeyError(
                    f"Boundary type '{boundary_type}' at '{location}' requires '{key}'"
                )

        # some Boundary subclass needs to be passed location and ndim at construct, see
        # SuperMagnetofastInflowBoundary for instance
        if "location" in signature(ctor).parameters:
            resolved[location] = ctor(**params, location=location, ndim=ndim)
        else:
            resolved[location] = ctor(**params)

    return periodicities, resolved


def _add_inflow_vector(dp, bc_path, key, vec):
    """Serialise a prescribed inflow 3-vector, one double per component."""
    for axis, c in zip("xyz", list(vec), strict=True):
        dp.add_double(f"{bc_path}/{key}/{axis}", c)
