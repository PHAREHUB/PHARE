from dataclasses import dataclass

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

    quantities: list[tuple[str, float]]

    def populate_dict(self, dp):
        dp.add_int(_TAGGING_PATH + "nbr_quantities", len(self.quantities))
        for i, (name, threshold) in enumerate(self.quantities):
            dp.add_string(f"{_TAGGING_PATH}Q{i}/name", name)
            dp.add_double(f"{_TAGGING_PATH}Q{i}/threshold", threshold)


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


_TAGGER_SPECS = {
    "default": (DefaultTagging, {}),
    "lohner": (LohnerTagging, {"reltol": float, "abstol": float}),
}


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
    threshold = float(threshold)

    spec = kwargs.get("tagging", None)
    if spec is None:
        spec = {}
    if not isinstance(spec, dict):
        raise TypeError(
            "Error: 'tagging' must be a dict {'method':..., 'quantities':...}"
        )

    valid_spec_keys = {"method", "quantities", "params"}
    unknown_keys = set(spec) - valid_spec_keys
    if unknown_keys:
        raise ValueError(
            f"Error: unknown tagging keys {sorted(unknown_keys)},"
            f" expected among {sorted(valid_spec_keys)}"
        )

    method = spec.get("method", "default")
    if method not in _TAGGER_SPECS:
        raise ValueError(
            f"Error: invalid tagging method '{method}', expected one of {tuple(_TAGGER_SPECS)}"
        )
    cls, param_coercions = _TAGGER_SPECS[method]

    params = spec.get("params", None) or {}
    if not isinstance(params, dict):
        raise TypeError("Error: tagging 'params' must be a dict {name: value}")
    unknown = set(params) - set(param_coercions)
    if unknown:
        raise ValueError(
            f"Error: invalid tagging params {sorted(unknown)} for method '{method}',"
            f" expected keys among {tuple(param_coercions)}"
        )
    param_kwargs = {
        name: param_coercions[name](value) for name, value in params.items()
    }

    quantities = spec.get("quantities", None)
    if quantities:
        items = quantities.items() if isinstance(quantities, dict) else quantities
        quantities = [(str(name), float(thr)) for name, thr in items]
    else:
        quantities = [("B", threshold)]

    not_available = [name for name, _ in quantities if name in _NOT_AVAILABLE]
    if not_available:
        raise ValueError(
            f"Error: tagging quantities {not_available} are not available for tagging yet"
        )

    return cls(quantities=quantities, **param_kwargs)  # pyright: ignore[reportArgumentType]
