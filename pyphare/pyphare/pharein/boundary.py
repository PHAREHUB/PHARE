"""
Boundary-condition resolution and validation for pharein.Simulation.

'boundary_conditions' (the public Simulation constructor option) is a per-location dict:
  {"xlower": {"type": "open"}, "xupper": {"type": "super-magnetofast-inflow", "data": {...}}, ...}

resolve_boundary_conditions(ndim, **kwargs) is the single entry point used by Simulation's
checker() pipeline; it returns a dict[location] -> BoundaryCondition stored as
Simulation.boundary_conditions.
"""

import math
import numbers
from abc import ABC
from dataclasses import InitVar, dataclass

from ..core import phare_utilities


@dataclass
class BoundaryCondition(ABC):
    """Base class holding the serialisation entry point common to every BC type."""

    type = None

    def populate_dict(self, bc_path, ndim):
        _add_string(f"{bc_path}/type", self.type)


@dataclass
class NoneBC(BoundaryCondition):
    type = "none"


@dataclass
class OpenBC(BoundaryCondition):
    type = "open"


@dataclass
class ReflectiveBC(BoundaryCondition):
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
class SuperMagnetofastInflowBC(BoundaryCondition):
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

    def populate_dict(self, bc_path, ndim):
        super().populate_dict(bc_path, ndim)
        _add_inflow_scalar(bc_path, "density", self.density, ndim)
        _add_inflow_scalar(bc_path, "pressure", self.pressure, ndim)
        _add_inflow_vector(bc_path, "velocity", self.velocity, ndim)
        _add_inflow_vector(bc_path, "B", self.B, ndim)


# ------------------------------------------------------------------------------

_TYPE_CTORS = {
    "none": NoneBC,
    "open": OpenBC,
    "reflective": ReflectiveBC,
    "super-magnetofast-inflow": SuperMagnetofastInflowBC,
}
_REQUIRED_DATA_KEYS = {
    "super-magnetofast-inflow": ("density", "pressure", "velocity", "B"),
}


def resolve_boundary_conditions(ndim, **kwargs):
    all_directions = ["x", "y", "z"][:ndim]
    sides = ("lower", "upper")
    boundary_types = kwargs["boundary_types"]

    physical_directions = [
        direction
        for direction, boundary_type in zip(all_directions, boundary_types)
        if boundary_type == "physical"
    ]
    physical_boundary_locations = [
        f"{direction}{side}" for direction in physical_directions for side in sides
    ]

    model_options = phare_utilities.listify(kwargs.get("model_options", "HybridModel"))
    if physical_boundary_locations and "MHDModel" not in model_options:
        raise ValueError(
            "'physical' boundary_types are only supported by the MHDModel; "
            f"got model_options={model_options}"
        )

    all_boundary_locations = [
        f"{direction}{side}" for side in sides for direction in all_directions
    ]
    raw = kwargs.get("boundary_conditions", {})
    if not isinstance(raw, dict):
        raise TypeError("A dict should be passed to argument 'boundary_conditions'")

    for location in raw:
        if location not in all_boundary_locations:
            raise ValueError(
                f"Wrong boundary name {location}: should belong to {all_boundary_locations}"
            )

    resolved = {}
    for location in all_boundary_locations:
        bc = raw.get(location, {"type": "none"})
        if not isinstance(bc, dict):
            raise TypeError(
                f"A dict should be passed to the boundary {location} for specifying a "
                f"boundary condition"
            )
        if "type" not in bc:
            raise KeyError(
                f"No key 'type' found in the boundary_condition dict passed to {location}"
            )
        boundary_type = bc["type"]
        if boundary_type not in _TYPE_CTORS:
            raise ValueError(
                f"Boundary type {boundary_type} is not valid: it should belong to "
                f"{tuple(_TYPE_CTORS)}"
            )
        if location in physical_boundary_locations and boundary_type == "none":
            raise KeyError(
                f"{location} is a physical boundary and should be provided with a valid "
                f"type other than 'none'."
            )
        if location not in physical_boundary_locations and boundary_type != "none":
            raise ValueError(
                f"{location} is not a physical boundary (its direction is periodic) and must "
                f"have type 'none', got '{boundary_type}'"
            )

        data = bc.get("data", {})
        if boundary_type in _REQUIRED_DATA_KEYS.keys():
            for key in _REQUIRED_DATA_KEYS[boundary_type]:
                if key not in data:
                    raise KeyError(
                        f"Boundary type '{boundary_type}' at '{location}' requires '{key}' "
                        f"inside 'data'"
                    )
            resolved[location] = _TYPE_CTORS[boundary_type](
                **data, location=location, ndim=ndim
            )
        else:
            if "data" in bc:
                raise ValueError(
                    f"Boundary type '{boundary_type}' at '{location}' takes no 'data' block, "
                    f"but one was provided. Supported data-carrying types: "
                    f"{sorted(_REQUIRED_DATA_KEYS.keys())}"
                )
            resolved[location] = _TYPE_CTORS[boundary_type]()

    return resolved


def _add_string(path, val):
    import pybindlibs.dictator as pp

    pp.add_string(path, val)


def _add_double(path, val):
    import pybindlibs.dictator as pp

    pp.add_double(path, float(val))


def _add_inflow_scalar(bc_path, key, val, ndim):
    """Serialise a prescribed inflow scalar."""
    _add_double(f"{bc_path}/data/{key}", val)


def _add_inflow_vector(bc_path, key, vec, ndim):
    """Serialise a prescribed inflow 3-vector, one double per component."""
    for axis, c in zip("xyz", list(vec), strict=True):
        _add_double(f"{bc_path}/data/{key}/{axis}", c)
