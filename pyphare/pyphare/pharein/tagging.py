import math
import numbers
from dataclasses import dataclass
from typing import ClassVar

_TAGGING_PATH = "simulation/AMR/refinement/tagging/"

_NOT_AVAILABLE = {
    "J",
    "Jx",
    "Jy",
    "Jz",
    "E",
    "Ex",
    "Ey",
    "Ez",
    "V",
    "Vx",
    "Vy",
    "Vz",
    "P",
    "velocity",
    "pressure",
    "bulkVel",
    "massDensity",
    "chargeDensity",
}


@dataclass
class Tagging:
    method = None
    quantity_parameters: ClassVar[dict[str, float | None]] = {}

    quantities: list[tuple[str, dict[str, float]]]

    @classmethod
    def check_parameters(cls, name, parameters):
        pass

    def populate_dict(self, dp):
        dp.add_int(_TAGGING_PATH + "nbr_quantities", len(self.quantities))
        for i, (name, parameters) in enumerate(self.quantities):
            dp.add_string(f"{_TAGGING_PATH}Q{i}/name", name)
            for key, value in parameters.items():
                dp.add_double(f"{_TAGGING_PATH}Q{i}/{key}", value)


@dataclass
class DefaultTagging(Tagging):
    method = "default"


@dataclass
class LohnerTagging(Tagging):
    quantity_parameters: ClassVar[dict[str, float | None]] = {
        "reltol": 0.01,
        "ref_magnitude": 1e-8,
    }

    @classmethod
    def check_parameters(cls, name, parameters):
        for key in ("reltol", "ref_magnitude"):
            if not parameters[key] >= 0:
                raise ValueError(
                    f"Error: tagging quantity '{name}': {key} must be >= 0,"
                    f" got {parameters[key]}"
                )


@dataclass
class LohnerPerDirectionTagging(LohnerTagging):
    method = "lohner-per-direction"


@dataclass
class LohnerCrossDerivativesTagging(LohnerTagging):
    method = "lohner-cross-derivatives"


@dataclass
class GradientScaleTagging(Tagging):
    method = "gradient-scale"
    quantity_parameters: ClassVar[dict[str, float | None]] = {"ref_magnitude": None}

    @classmethod
    def check_parameters(cls, name, parameters):
        if not parameters["ref_magnitude"] > 0:
            raise ValueError(
                f"Error: tagging quantity '{name}': ref_magnitude must be > 0,"
                f" got {parameters['ref_magnitude']}"
            )


_TAGGERS = {
    "default": DefaultTagging,
    "lohner-per-direction": LohnerPerDirectionTagging,
    "lohner-cross-derivatives": LohnerCrossDerivativesTagging,
    "gradient-scale": GradientScaleTagging,
}


def _as_float(name, key, value):
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(
            f"Error: tagging quantity '{name}': {key} must be a number, got {value!r}"
        )
    if not math.isfinite(value):
        raise ValueError(
            f"Error: tagging quantity '{name}': {key} must be finite, got {value!r}"
        )
    return float(value)


def _quantity_parameters(cls, name, spec):
    allowed = {"threshold", *cls.quantity_parameters}
    if not isinstance(spec, dict):
        spec = {"threshold": spec}
    unknown = set(spec) - allowed
    if unknown:
        raise ValueError(
            f"Error: invalid keys {sorted(unknown)} for tagging quantity '{name}'"
            f" with method '{cls.method}', expected keys among {sorted(allowed)}"
        )
    required = {"threshold"} | {
        key for key, default in cls.quantity_parameters.items() if default is None
    }
    missing = sorted(required - set(spec))
    if missing:
        raise ValueError(
            f"Error: tagging quantity '{name}' needs {missing} with method '{cls.method}'"
        )
    parameters = {**cls.quantity_parameters, **spec}
    parameters = {
        key: _as_float(name, key, parameters[key])
        for key in ("threshold", *cls.quantity_parameters)
    }
    cls.check_parameters(name, parameters)
    return parameters


def _list_entry(entry):
    if (
        isinstance(entry, (str, bytes))
        or not isinstance(entry, (list, tuple))
        or len(entry) != 2
        or not isinstance(entry[0], str)
    ):
        raise TypeError(
            "Error: tagging 'quantities' list entries must be"
            f" (name, threshold or dict) pairs, got {entry!r}"
        )
    return entry[0], entry[1]


def resolve_tagging(**kwargs):
    if kwargs.get("refinement", "boxes") != "tagging":
        if kwargs.get("tagging", None) is not None:
            raise ValueError(
                "Error: a 'tagging' configuration was provided but refinement is not"
                " 'tagging'; set refinement='tagging' to enable it (or drop the 'tagging' kwarg)"
            )
        return None

    threshold = kwargs.get("tagging_threshold", 0.1)
    if threshold is None:
        threshold = 0.1
    if isinstance(threshold, bool) or not isinstance(threshold, numbers.Real):
        raise TypeError(f"Error: tagging_threshold must be a number, got {threshold!r}")

    spec = kwargs.get("tagging", None)
    if spec is None:
        spec = {}
    if not isinstance(spec, dict):
        raise TypeError(
            "Error: 'tagging' must be a dict {'method':..., 'quantities':...}"
        )

    valid_spec_keys = {"method", "quantities"}
    unknown_keys = set(spec) - valid_spec_keys
    if unknown_keys:
        raise ValueError(
            f"Error: unknown tagging keys {sorted(unknown_keys)},"
            f" expected among {sorted(valid_spec_keys)}"
        )

    method = spec.get("method", "default")
    if method not in _TAGGERS:
        raise ValueError(
            f"Error: invalid tagging method '{method}', expected one of {tuple(_TAGGERS)}"
        )
    cls = _TAGGERS[method]

    quantities = spec.get("quantities", None)
    if quantities:
        if isinstance(quantities, dict):
            quantities = [(str(name), value) for name, value in quantities.items()]
        else:
            quantities = [_list_entry(entry) for entry in quantities]
        names = [name for name, _ in quantities]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ValueError(
                f"Error: tagging quantities {duplicates} given more than once"
            )
    else:
        required = sorted(k for k, v in cls.quantity_parameters.items() if v is None)
        if required:
            raise ValueError(
                f"Error: tagging method '{method}' needs explicit 'quantities',"
                f" each with its {required}"
            )
        quantities = [("B", threshold)]

    not_available = [name for name, _ in quantities if name in _NOT_AVAILABLE]
    if not_available:
        raise ValueError(
            f"Error: tagging quantities {not_available} are not available for tagging yet"
        )

    return cls(
        quantities=[
            (name, _quantity_parameters(cls, name, value)) for name, value in quantities
        ]
    )
