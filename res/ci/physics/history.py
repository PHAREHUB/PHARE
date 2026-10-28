"""Night records and per-test state for the nightly page.

Pure module: no network, no import of publish. A night is
{"date", "run": {"sha", "manual", "time", "builds"}, "tests": {name: record}}; the
history is the list of nights, oldest first, kept in history.json.
"""

import json
import math
import re
import time
from pathlib import Path

DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
SHA = re.compile(r"[0-9a-f]{40}")
BUILD_TYPE = re.compile(r"[A-Za-z0-9_]+")
JOB_NAME = re.compile(r"[a-z0-9_]+")

ABNORMAL = ("own check failed", "suspicious")
NOMINAL = ("looks right",)
NEUTRAL = ("cannot tell", "not judged", "no analysis", "own check unknown", "no data")
NOT_RUN = ("not run",)

VERSION = 1
MAX_STATUS, MAX_REASON = 40, 200


def number(x):
    """x as a finite float, or None: ints accepted, bools and strings refused, and
    ints beyond float range are invalid rather than an error."""
    if type(x) not in (int, float):
        return None
    try:
        f = float(x)
    except OverflowError:
        return None
    return f if math.isfinite(f) else None


def _short(text):
    return None if text is None else str(text)[:MAX_REASON]


def test_record(deck, cfg, part, model):
    """The night's record for one deck: the first matching rule wins."""
    job = part.get("job")
    passed = deck.get("deck_passed") if deck else None
    rec = {
        "status": None,
        "job": job,
        "deck_passed": passed,
        "key": None,
        "reason": None,
    }

    def done(status, reason=None, **flags):
        rec["status"], rec["reason"] = status, _short(reason)
        rec.update({k: True for k, v in flags.items() if v})
        return rec

    # Dormant in production: nothing writes decks_results.json yet, so deck_passed is
    # never False until a per-deck result source exists.
    if passed is False:
        return done("own check failed")
    status = deck.get("status") if deck else None
    if part.get("missing") or deck is None or status == "output missing":
        if part.get("missing"):
            return done("not run", f"the {job} job left no results tonight")
        if deck is None:
            return done("not run")
        return done("not run", status)
    if isinstance(status, str) and (
        status.startswith("extract error") or status == "stale output ignored"
    ):
        return done("no data", status)
    rec["key"] = _key(deck, cfg)
    if cfg.get("judge") is False:
        return done("not judged", cfg.get("judge_reason"))
    summary = part.get("summary") or {}
    rc = summary.get("decks_exit_code")
    steps_failed = rc not in (0, None)
    meta = part.get("meta") or {}
    if meta.get("agent_status") != "ok" or not part.get("agent"):
        return done(
            "no analysis",
            meta.get("agent_reason"),
            own_check_unknown=passed is None and steps_failed,
        )
    verdict = model[0] if model else "cannot tell"
    if verdict not in ("looks right", "suspicious", "cannot tell"):
        verdict = "cannot tell"
    if passed is None and steps_failed and verdict in ("looks right", "cannot tell"):
        return done("own check unknown", f"tests step failed (exit {rc})")
    return done(
        verdict,
        own_check_unknown=passed is None and steps_failed and verdict == "suspicious",
    )


def _key(deck, cfg):
    name = cfg.get("key_number")
    numbers = deck.get("numbers")
    if not name or not isinstance(numbers, dict):
        return None
    return number(numbers.get(name))


def night(date, run, tests):
    return {"date": date, "run": run, "tests": tests}


def add_night(history, n):
    """(history, recorded). A manual night never replaces a scheduled one."""
    history = list(history)
    for i, old in enumerate(history):
        if old["date"] == n["date"]:
            if n["run"].get("manual") and not old["run"].get("manual"):
                return history, False
            history[i] = n
            return history, True
    history.append(n)
    history.sort(key=lambda x: x["date"])
    return history, True


def merge_pending(history, pending):
    for n in pending:
        history, _ = add_night(history, n)
    return history


def _kind(status):
    if status in ABNORMAL:
        return "abnormal"
    if status in NOMINAL:
        return "nominal"
    return "neutral"  # NEUTRAL, NOT_RUN and any unknown string


def state(history, name, judged_tonight):
    """The test's state as of the newest night."""
    empty = {
        "status": "no data",
        "abnormal": False,
        "since": None,
        "new_tonight": False,
        "last_judged": None,
        "nights_since_judged": None,
        "prev_sha": None,
    }
    if not history:
        return empty
    last = len(history) - 1
    newest = history[last]["tests"].get(name)
    status = newest["status"] if newest else "not run"
    shown = (
        status
        if status in ABNORMAL + NOMINAL + NEUTRAL + NOT_RUN
        else (f"unknown status {status}")
    )
    collected = []  # (index, record), newest first
    recovered = False  # a passing own check seen, newer than the nights still to read
    for i in range(last, -1, -1):
        rec = history[i]["tests"].get(name)
        if rec is None:
            continue
        kind = _kind(rec["status"])
        if kind == "nominal":
            break
        if kind == "abnormal":
            if recovered and rec["status"] == "own check failed":
                break  # that run was ended by the passing own check
            recovered = False  # a model-started run needs "looks right" to end
            collected.append((i, rec))
        elif rec.get("deck_passed") is True and all(
            r["status"] == "own check failed" for _, r in collected
        ):
            recovered = True
    newest_kind = _kind(status)
    abnormal = newest_kind == "abnormal" or (
        newest_kind == "neutral" and bool(collected)
    )
    if (
        abnormal
        and newest_kind == "neutral"
        and not judged_tonight
        and collected[-1][1]["status"] == "suspicious"
    ):
        abnormal = False
    out = dict(empty, status=shown, abnormal=abnormal)
    if abnormal:
        out["since"] = history[collected[-1][0]]["date"]
        out["new_tonight"] = out["since"] == history[last]["date"] and len(history) > 1
    for i in range(last, -1, -1):
        rec = history[i]["tests"].get(name)
        if rec and _kind(rec["status"]) != "neutral":
            out["last_judged"] = history[i]["date"]
            out["nights_since_judged"] = last - i
            break
    out["prev_sha"] = next(
        (n["run"]["sha"] for n in reversed(history[:last]) if n["run"].get("sha")),
        None,
    )
    return out


def _valid_record(r):
    if not isinstance(r, dict):
        return False
    s = r.get("status")
    if not isinstance(s, str) or len(s) > MAX_STATUS:
        return False
    job = r.get("job")
    if job is not None and not (isinstance(job, str) and JOB_NAME.fullmatch(job)):
        return False
    passed = r.get("deck_passed")
    if passed is not None and type(passed) is not bool:
        return False
    k = r.get("key")
    if k is not None and number(k) is None:
        return False
    reason = r.get("reason")
    if reason is not None and (not isinstance(reason, str) or len(reason) > MAX_REASON):
        return False
    return all(type(r[f]) is bool for f in ("expired", "own_check_unknown") if f in r)


def _valid_build(b):
    return (
        isinstance(b, dict)
        and type(b.get("id")) is int
        and isinstance(b.get("type"), str)
        and BUILD_TYPE.fullmatch(b["type"]) is not None
        and isinstance(b.get("start"), str)
        and DATE.fullmatch(b["start"]) is not None
    )


def _valid_night(n):
    """(night with only its good records, number of bad records), or None if bad."""
    if not isinstance(n, dict):
        return None
    d, run, tests = n.get("date"), n.get("run"), n.get("tests")
    if not (isinstance(d, str) and DATE.fullmatch(d)):
        return None
    if not isinstance(run, dict) or not isinstance(tests, dict):
        return None
    sha = run.get("sha")
    if sha is not None and not (isinstance(sha, str) and SHA.fullmatch(sha)):
        return None
    if type(run.get("manual")) is not bool:
        return None
    builds = run.get("builds")
    if not isinstance(builds, dict) or not all(
        isinstance(j, str) and _valid_build(b) for j, b in builds.items()
    ):
        return None
    good = {k: r for k, r in tests.items() if _valid_record(r)}
    return dict(n, tests=good), len(tests) - len(good)


def load(path):
    """(nights, notes, clean). Anything dropped makes clean False."""
    path = Path(path)
    if not path.exists():
        return [], [], True
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError) as e:
        return [], [f"{path.name} unreadable ({type(e).__name__})"], False
    if (
        not isinstance(data, dict)
        or data.get("version") != VERSION
        or not isinstance(data.get("nights"), list)
    ):
        return [], [f"{path.name} has the wrong shape"], False
    raw = data["nights"]
    nights, bad_nights, bad_records = [], 0, 0
    for n in raw:
        v = _valid_night(n)
        if v is None:
            bad_nights += 1
        else:
            nights.append(v[0])
            bad_records += v[1]
    if bad_nights * 2 > len(raw):
        return [], [f"{path.name}: {bad_nights} of {len(raw)} nights bad"], False
    notes = []
    if bad_nights:
        notes.append(f"{path.name}: dropped {bad_nights} bad night(s)")
    if bad_records:
        notes.append(f"{path.name}: dropped {bad_records} bad record(s)")
    return nights, notes, not notes


def load_pending(path, date):
    """(nights, notes) from history.pending.json; an unreadable file is renamed
    history.pending.broken-<date>.json so it is not lost and not retried."""
    path = Path(path)
    nights, notes, clean = load(path)
    if not clean and not nights:
        broken = path.with_name(f"{path.stem}.broken-{date}{path.suffix}")
        try:
            path.rename(broken)
            notes.append(f"set aside as {broken.name}")
        except OSError as e:
            notes.append(f"could not set aside {path.name} ({type(e).__name__})")
    return nights, notes


def dump(nights, path):
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps({"version": VERSION, "nights": nights}, indent=1) + "\n")
    tmp.replace(path)


def canary(tonight_builds, exists):
    """True if any (build id, job) of tonight still answers True: TeamCity then shows
    outputs, so a False elsewhere means expired rather than a broken lookup."""
    return any(exists(b, j) is True for b, j in tonight_builds)


def expire(history, exists, canary_ok, drawn, deadline_s=120, clock=time.monotonic):
    """(history, notes). Marks records expired when TeamCity says their build's
    physics-out is gone (exists -> False); None (unknown) never expires anything.
    Only the last `drawn` nights are asked, each (build, job) once overall (a reused
    build appears in many nights), until the deadline. A job that is not a plain
    name is never asked: it would go into a URL."""
    history = [
        dict(n, tests={k: dict(r) for k, r in n["tests"].items()}) for n in history
    ]
    if not canary_ok:
        return history, [
            "retention not checked: TeamCity did not show tonight's outputs"
        ]
    stop = clock() + deadline_s
    asked = {}
    for n in history[-drawn:] if drawn > 0 else []:
        builds = n["run"]["builds"]
        for rec in n["tests"].values():
            job = rec.get("job")
            b = builds.get(job) if isinstance(job, str) else None
            if b is None or rec.get("expired") or not JOB_NAME.fullmatch(job):
                continue
            key = (b["id"], job)
            if key not in asked:
                if clock() > stop:
                    return history, ["retention check stopped at the deadline"]
                asked[key] = exists(*key)
            if asked[key] is False:
                rec["expired"] = True
    return history, []
