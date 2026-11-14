"""Load res/ci/physics/decks.yml: the decks the nightly physics report covers, per job."""

import os
import re

import yaml

EXTRACTORS = {"ion_ion_beam", "tangential_discontinuity", "figures_only"}
DECK_KEYS = {"name", "source", "workdir", "figures", "extractor", "theory"}
# optional deck keys: own_check, judge (false: extracted, never shown to the agent) with its
# judge_reason, max_figures, key_number (a key of the deck's numbers, charted on the page)
# own_check "none": the deck asserts nothing, so the report's verdict is the only check it gets
OWN_CHECKS = ("deck asserts", "none")
DEFAULT_MAX_FIGURES = 6  # free flash model, one request: keep the image count small
# the page's branch reaches git: nothing that could pass for an option or a path
PAGE_BRANCH = re.compile(r"(?!-)[A-Za-z0-9._/-]+")
DEFAULT_PATH = os.path.join(os.path.dirname(__file__), "decks.yml")


def load(path=None):
    with open(path or DEFAULT_PATH) as f:
        cfg = yaml.safe_load(f)
    check_page(cfg.get("page"))
    for job_name, job_decks in cfg["jobs"].items():
        for deck in job_decks:
            missing = DECK_KEYS - deck.keys()
            if missing:
                raise ValueError(
                    f"{job_name}/{deck.get('name')}: missing {sorted(missing)}"
                )
            if deck["extractor"] not in EXTRACTORS:
                raise ValueError(
                    f"{job_name}/{deck['name']}: unknown extractor {deck['extractor']}"
                )
            check_optional(f"{job_name}/{deck['name']}", deck)
    return cfg


def check_page(page):
    if type(page) is not dict or set(page) != {"branch", "dots"}:
        raise ValueError("page: needs exactly branch and dots")
    if type(page["branch"]) is not str or not PAGE_BRANCH.fullmatch(page["branch"]):
        raise ValueError(f"page: {page['branch']!r} is not a valid branch name")
    if type(page["dots"]) is not int or page["dots"] < 1:
        raise ValueError(f"page: dots must be a positive integer, not {page['dots']!r}")


def check_optional(where, deck):
    if deck.get("own_check", "deck asserts") not in OWN_CHECKS:
        raise ValueError(f"{where}: own_check must be one of {OWN_CHECKS}")
    if type(deck.get("judge", True)) is not bool:
        raise ValueError(f"{where}: judge must be true or false")
    reason = deck.get("judge_reason")
    if reason is not None and type(reason) is not str:
        raise ValueError(f"{where}: judge_reason must be a string")
    if deck.get("judge", True) is False and not reason:
        raise ValueError(f"{where}: judge: false needs a judge_reason")
    key = deck.get("key_number")
    if key is not None and (type(key) is not str or not key):
        raise ValueError(f"{where}: key_number must be a non-empty string")
    cap = deck.get("max_figures", DEFAULT_MAX_FIGURES)
    if type(cap) is not int or cap < 1:  # type(): a YAML 'true' is an int to isinstance
        raise ValueError(
            f"{where}: max_figures must be a positive integer, not {cap!r}"
        )


def job(cfg, name):
    return cfg["jobs"][name]
