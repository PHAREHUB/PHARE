"""
Tagging (AMR refinement criterion) resolution and validation for pharein.Simulation.

The public `tagging` Simulation option selects a refinement-tagging criterion and the
per-quantity thresholds it is applied to:

    refinement="tagging",
    tagging={"method": "lohner",
             "quantities": {"B": 0.1, "massDensity": 0.4},
             "params": {"reltol": 0.02}}

where `quantities` maps a field name to its own threshold (a cell is tagged if ANY
quantity's indicator exceeds its threshold). All keys are optional: `method` defaults to
"default"; when `quantities` is omitted (or no `tagging` dict is given at all) the criterion
falls back to B at `tagging_threshold` (kept for backward compatibility, default 0.1).

`params` carries method-specific numbers, one dataclass field per key:
  - lohner: `reltol` (Loehner's eps, default 0.02) weights the noise filter in the
    denominator of the estimator (damps the indicator where the field variation is small
    relative to its magnitude); `abstol` (default 1e-30) is an absolute denominator floor.
  - wavelet: multiresolution detail criterion (Domingues et al. 2019,
    10.1016/j.compfluid.2019.06.025). The indicator is the ABSOLUTE prediction error
    |Q - Qpredicted| (units of Q), so thresholds are not comparable with default/lohner
    ones. By default the threshold follows Harten's level scaling eps_l = eps / 2^{dim (L - l)}
    (L the finest level), refining coarse levels more eagerly; `level_scaling` (bool,
    default True) set to False uses the per-quantity threshold unscaled on every level.

Quantity names are not restricted here: the C++ tagger resolves them against the model's
field tree at setup and throws (listing the available names) if a name matches nothing.
Use the model's own field names -- e.g. the hybrid density is `chargeDensity`/`massDensity`
and the MHD density is `rho`; `B` is a universal alias for the magnetic-field components in
either model.

resolve_tagging(**kwargs) is the single entry point used by Simulation's checker() pipeline;
it returns the Tagging instance stored as Simulation.tagging (or None when refinement is not
"tagging").
"""

from abc import ABC
from dataclasses import dataclass
from typing import List, Tuple


_TAGGING_PATH = "simulation/AMR/refinement/tagging/"


@dataclass
class Tagging(ABC):
    """
    Base class holding what every tagging method shares: the (name, threshold) quantity list.
    """

    method = None

    quantities: List[Tuple[str, float]]

    def populate_dict(self, dp):
        """Mirror the public `tagging` dict shape (method + quantities + per-method params) on
        the C++ side.

        `dp` is a dict populator (see pharein.initialize.general.dict_populator): an object
        exposing add_string/add_double/add_int/add_bool - passed in rather than imported, to
        avoid a circular import between this module and pharein.initialize.general.
        """
        dp.add_string(_TAGGING_PATH + "method", self.method)
        dp.add_int(_TAGGING_PATH + "nbr_quantities", len(self.quantities))
        for i, (name, threshold) in enumerate(self.quantities):
            q_path = f"{_TAGGING_PATH}Q{i}/"
            dp.add_string(q_path + "name", name)
            dp.add_double(q_path + "threshold", threshold)


@dataclass
class DefaultTagging(Tagging):
    method = "default"


@dataclass
class LohnerTagging(Tagging):
    method = "lohner"

    reltol: float = 0.02
    abstol: float = 1e-30

    def populate_dict(self, dp):
        super().populate_dict(dp)
        dp.add_double(_TAGGING_PATH + "params/reltol", self.reltol)
        dp.add_double(_TAGGING_PATH + "params/abstol", self.abstol)


@dataclass
class WaveletTagging(Tagging):
    method = "wavelet"

    level_scaling: bool = True

    def populate_dict(self, dp):
        super().populate_dict(dp)
        dp.add_bool(_TAGGING_PATH + "params/level_scaling", self.level_scaling)


# method name -> (dataclass, {param: coercion}). The coercion map doubles as the set of
# params accepted by each method (used to reject typos with a helpful message).
_TAGGERS = {
    "default": (DefaultTagging, {}),
    "lohner": (LohnerTagging, {"reltol": float, "abstol": float}),
    "wavelet": (WaveletTagging, {"level_scaling": bool}),
}


def resolve_tagging(**kwargs):
    """
    Resolve the public 'tagging' / 'refinement' / 'tagging_threshold' Simulation options into a
    validated Tagging instance, or None when refinement is not "tagging".
    """
    if kwargs.get("refinement", "boxes") != "tagging":
        # a fully specified tagging= kwarg with refinement left at its default is almost
        # always a mistake: fail loudly instead of silently ignoring the configuration.
        if kwargs.get("tagging", None) is not None:
            raise ValueError(
                "Error: a 'tagging' configuration was provided but refinement is not"
                " 'tagging'; set refinement='tagging' to enable it (or drop the 'tagging' kwarg)"
            )
        return None

    threshold = kwargs.get("tagging_threshold", 0.1)
    if threshold is None:
        threshold = 0.1
    threshold = float(threshold)

    spec = kwargs.get("tagging", None)
    if spec is None:
        spec = {}
    if not isinstance(spec, dict):
        raise ValueError("Error: 'tagging' must be a dict {'method':..., 'quantities':...}")

    # reject typos in the top-level keys (e.g. "quantites") instead of silently swallowing
    # them and falling back to the B/tagging_threshold default.
    valid_spec_keys = {"method", "quantities", "params"}
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
    cls, coercions = _TAGGERS[method]

    params = spec.get("params", None) or {}
    if not isinstance(params, dict):
        raise ValueError("Error: tagging 'params' must be a dict {name: value}")
    unknown = set(params) - set(coercions)
    if unknown:
        raise ValueError(
            f"Error: invalid tagging params {sorted(unknown)} for method '{method}',"
            f" expected keys among {tuple(coercions)}"
        )
    # each param becomes a constructor kwarg, coerced to the dataclass field's type.
    param_kwargs = {name: coercions[name](value) for name, value in params.items()}

    quantities = spec.get("quantities", None)
    if quantities:
        # accept the {name: threshold} mapping (preferred) or an iterable of
        # (name, threshold) pairs.
        items = quantities.items() if isinstance(quantities, dict) else quantities
        quantities = [(str(name), float(thr)) for name, thr in items]
    else:
        # no quantities specified -> criterion applies to B at tagging_threshold
        quantities = [("B", threshold)]

    return cls(quantities=quantities, **param_kwargs)
