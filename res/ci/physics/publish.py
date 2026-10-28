"""Nightly physics report, job 'publish': merge the job reports into the nightly page. No LLM here.

Writes index.html and history.json into a clone of the nightly page repo (--site), and
with --push commits and pushes them there. It posts nothing to PHARE's repo: no comment,
no commit status. The only token is NIGHTLY_PAGE_TOKEN, read only with --push, and it
reaches git through the environment alone.

Tonight's builds (ids, types, start dates, commit, trigger, default branch) all come from
TeamCity's guest API. meta.json and summary.json come from the pod that runs code, so
clean_part checks every field before the page sees it; what they say about the commit or
the build id is only cross-checked, never used. The agent's text only ever becomes escaped
plain text on the page. Never a gate: whatever happens, exit 0, and an existing index.html
is left alone when something goes wrong.
"""

import argparse
import base64
import datetime
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from zoneinfo import ZoneInfo

import decks
import history
import judge
import page
import tc

SAFE_FIGURE = re.compile(r"[A-Za-z0-9._-]+")  # as extract's SAFE_NAME, re-checked here


def sections(body, names):
    """The agent's text as (preamble, {deck: (verdict, rest)}): a heading counts only for a
    known deck, once; the verdict is the line right after it, if it is exactly one of the
    model's three, else 'cannot tell' (and the line stays in the text)."""
    pre, secs, current = [], {}, None
    for line in body.split("\n"):
        line = line.rstrip()
        m = re.fullmatch(r"### (\S+)", line)
        if m and m.group(1) in names and m.group(1) not in secs:
            current = m.group(1)
            secs[current] = [None, []]
        elif current is None:
            pre.append(line)
        elif secs[current][0] is None:
            ok = line.startswith("Verdict: ") and line[9:] in judge.VERDICTS
            secs[current][0] = line[9:] if ok else "cannot tell"
            if not ok:
                secs[current][1].append(line)
        else:
            secs[current][1].append(line)
    return "\n".join(pre), {
        n: (v or "cannot tell", "\n".join(r)) for n, (v, r) in secs.items()
    }


NUMBER_KEY = re.compile(r"[A-Za-z0-9_]{1,40}")
TEXT_CAP = 200
NOT_GIVEN_MAX = 20
FIGURE_SUFFIXES = (".png", ".mp4")
AGENT_STATUSES = ("ok", "unavailable", "not run")


def str_matching(value, pattern):
    return value if isinstance(value, str) and pattern.fullmatch(value) else None


def int_or_none(value):
    return value if type(value) is int else None


def finite(x):
    """True for a bool or a finite int/float (a huge int is not finite: float() overflows)."""
    if type(x) is bool:
        return True
    return history.number(x) is not None


def clean_numbers(numbers, name, dropped):
    if not isinstance(numbers, dict):
        if numbers is not None:
            dropped.append(f"numbers of {name}")
        return {}
    good = {}
    for k, v in numbers.items():
        if isinstance(k, str) and NUMBER_KEY.fullmatch(k) and finite(v):
            good[k] = v
        else:
            dropped.append(f"number {name}.{k}")
    return good


def short_text(value, what, dropped):
    """A str cut to TEXT_CAP characters; anything else is dropped (None stays quiet)."""
    if isinstance(value, str):
        return value[:TEXT_CAP]
    if value is not None:
        dropped.append(what)
    return None


def clean_not_given(value, dropped):
    """The names the agent was not given: at most NOT_GIVEN_MAX of them, TEXT_CAP long."""
    if not isinstance(value, list):
        return []
    names = [n[:TEXT_CAP] for n in value if isinstance(n, str)]
    if len(names) > NOT_GIVEN_MAX:
        dropped.append(f"{len(names) - NOT_GIVEN_MAX} more not_given names")
    return names[:NOT_GIVEN_MAX]


def clean_part(part, caps):
    """The trust boundary: meta.json and summary.json are written in the pod that runs code.
    Every field publish renders or sends is checked; a bad one becomes None/'unknown' or is
    dropped, and what was dropped is listed (part['dropped'], shown as inline code).
    caps: {deck name: max_figures} for the job's decks in publish's own decks.yml."""
    meta, summary, dropped = part["meta"], part["summary"], []
    clean = {
        "build_id": int_or_none(meta.get("build_id")),
        "sha": str_matching(meta.get("sha"), history.SHA),
        "opencode_version": str_matching(meta.get("opencode_version"), judge.VERSION)
        or "unknown",
        "agent_status": meta.get("agent_status")
        if meta.get("agent_status") in AGENT_STATUSES
        else None,
        "agent_reason": short_text(meta.get("agent_reason"), "agent_reason", dropped),
        "unavailable_nights": int_or_none(meta.get("unavailable_nights")) or 0,
        "not_given": clean_not_given(meta.get("not_given"), dropped),
    }
    decks_in = summary.get("decks")
    decks_out, seen = [], set()
    for d in decks_in if isinstance(decks_in, list) else []:
        name = d.get("name") if isinstance(d, dict) else None
        if not isinstance(name, str) or name not in caps or name in seen:
            dropped.append(f"deck {name}")
            continue
        seen.add(name)
        figures = []
        for f in d.get("figures") if isinstance(d.get("figures"), list) else []:
            ok = (
                isinstance(f, str)
                and SAFE_FIGURE.fullmatch(f)
                and f.startswith(f"{name}__")
                and f.endswith(FIGURE_SUFFIXES)
            )
            if ok and len(figures) < caps[name]:
                figures.append(f)
            else:
                dropped.append(f"figure {f}" + (" over max_figures" if ok else ""))
        passed = d.get("deck_passed")
        decks_out.append(
            {
                "name": name,
                "figures": figures,
                "deck_passed": passed if type(passed) is bool else None,
                "freshness": short_text(
                    d.get("freshness"), f"freshness of {name}", dropped
                ),
                "status": short_text(d.get("status"), f"status of {name}", dropped),
                "numbers": clean_numbers(d.get("numbers"), name, dropped),
            }
        )
    rc = int_or_none(summary.get("decks_exit_code"))
    return part | {
        "meta": clean,
        "summary": {"decks_exit_code": rc, "decks": decks_out},
        "dropped": dropped,
    }


def read_part(d, job, caps):
    """One job's meta.json and summary.json (without both, 'no report') and agent.md, the
    agent's text alone (absent when the agent gave none), checked by clean_part."""
    try:
        part = {
            "job": job,
            "meta": json.loads((d / "meta.json").read_text()),
            "summary": json.loads((d / "summary.json").read_text()),
        }
    except (OSError, ValueError):
        return {"job": job, "missing": True}
    if not isinstance(part["meta"], dict) or not isinstance(part["summary"], dict):
        return {"job": job, "missing": True}
    try:
        part["agent"] = (d / "agent.md").read_text()
    except OSError:
        part["agent"] = ""
    return clean_part(part, caps)


def paris():
    """(timezone, notes): Europe/Paris, or UTC with a note when the system has no tzdata."""
    try:
        return ZoneInfo("Europe/Paris"), []
    except Exception:  # noqa: BLE001  ZoneInfoNotFoundError, KeyError, OSError: no tzdata
        return datetime.UTC, ["Europe/Paris unknown on this machine: UTC used"]


def own_build_id(a):
    """The page writer's TeamCity build id: --own-build-id, else the build properties."""
    if a.own_build_id is not None:
        return a.own_build_id
    value = tc.build_properties().get("teamcity.build.id", "")
    return int(value) if value.isdigit() else None


def jobs_tonight(base, cfg, own, a, get, notes):
    """(builds {job: build}, manual, the page writer's own build or None, the head build
    or None), all from TeamCity. A failed read gives no builds, never an error. A night
    whose trigger cannot be read is manual: only a schedule makes a scheduled night.
    A page writer outside the pipeline has no head: its job builds are the ones it
    downloaded, else (TeamCity names none when nothing was downloaded) each job's last
    finished build, and the head is the newest one depending on them."""
    try:
        head, sources = a.head_build_id, []
        if head is None and own is not None:
            head, others = tc.head_of(base, own, get)
            if head is None:
                sources = tc.artifact_sources(base, own, get)
                seen = {b["type"] for b in sources}
                for t in cfg["teamcity"].get("job_types", {}).values():
                    b = None if t in seen else tc.last_finished(base, t, get)
                    sources += [b] if b else []
                if sources:
                    newest = max(b["id"] for b in sources)
                    head, others = tc.head_of(base, newest, get)
            if others:
                notes.append(
                    f"other pipeline heads found: {', '.join(map(str, others))}"
                )
        start_from = head if head is not None else own
        if start_from is None:
            raise ValueError("no build id")
        builds = {}
        for b in sources or tc.upstream_builds(base, start_from, get):
            job = tc.job_of(b["type"], cfg)
            if job and (job not in builds or b["id"] > builds[job]["id"]):
                builds[job] = b
        info = tc.build_info(base, start_from, get)
        if own is None:
            mine = None
        elif own == start_from:
            mine = info
        else:
            mine = tc.build_info(base, own, get)
        manual = info is None or info["triggered"] != "schedule"
        return builds, manual, mine, info
    except Exception as e:  # noqa: BLE001  guest reads: any failure means unknown builds
        notes.append(f"TeamCity read failed ({type(e).__name__})")
        return {}, True, None, None


def common_sha(builds, notes):
    shas = {b["sha"] for b in builds.values()} - {None}
    if len(shas) > 1:
        notes.append("jobs disagree on the commit: commit unknown")
    return shas.pop() if len(shas) == 1 else None


def clean_parts(a, cfg):
    """The cleaned parts of --jobs (a missing one is {"job", "missing": True}); each other
    gets part["sections"] = {deck: (verdict, text)} from the agent's text, if it answered."""
    parts = []
    for job in a.jobs.split(","):
        decks_of_job = cfg["jobs"].get(job, [])
        caps = {
            d["name"]: d.get("max_figures", decks.DEFAULT_MAX_FIGURES)
            for d in decks_of_job
        }
        part = read_part(Path(a.inp, job), job, caps)
        if not part.get("missing"):
            asked = part["meta"]["agent_status"] == "ok"
            part["sections"] = sections(part["agent"], set(caps))[1] if asked else {}
        parts.append(part)
    return parts


def part_notes(part, cfg, builds):
    """What the page says about one job, from its cleaned part and TeamCity's build."""
    job, meta = part["job"], part["meta"]
    notes = []
    pinned = str(cfg.get("opencode_version"))
    if meta["opencode_version"] != pinned:
        notes.append(
            f"{job}: OpenCode {meta['opencode_version']} differs from the pinned {pinned}"
        )
    if meta["agent_status"] == "unavailable" and meta["unavailable_nights"] > 1:
        nights = meta["unavailable_nights"]
        notes.append(f"agent unavailable {nights} nights in a row ({job})")
    if any(d["freshness"] == "unknown" for d in part["summary"]["decks"]):
        notes.append(f"freshness unknown: figures may be from an earlier run ({job})")
    if meta["not_given"]:
        shown = ", ".join(meta["not_given"][:20])
        notes.append(
            f"{job}: not given to the agent (symlink or not a regular file): {shown}"
        )
    if part["dropped"]:
        shown = ", ".join(x[:80] for x in part["dropped"][:20])
        notes.append(
            f"{job}: not valid in meta.json or summary.json, not shown: {shown}"
        )
    build = builds.get(job)
    if build and meta["sha"] and build["sha"] and meta["sha"] != build["sha"]:
        notes.append(
            f"{job}: meta.json commit differs from TeamCity's, TeamCity's is used"
        )
    if build and meta["build_id"] not in (None, build["id"]):
        notes.append(
            f"{job}: meta.json build id differs from TeamCity's, TeamCity's is used"
        )
    return notes


def records(parts, cfg):
    """{deck name: night record} for every deck of decks.yml."""
    by_job = {p["job"]: p for p in parts}
    tests = {}
    for job, job_decks in cfg["jobs"].items():
        part = by_job.get(job) or {"job": job, "missing": True}
        decks_in = [] if part.get("missing") else part["summary"]["decks"]
        for deck_cfg in job_decks:
            name = deck_cfg["name"]
            deck = next((d for d in decks_in if d["name"] == name), None)
            model = part.get("sections", {}).get(name)
            tests[name] = history.test_record(deck, deck_cfg, part, model)
    return tests


def update_history(site, tonight, builds, cfg, base, get, notes):
    """(nights to show, tonight recorded, writes {file name: nights, None to delete}, files
    touched, both sides of a rename included). Nothing is written here."""
    date = tonight["date"]
    hist_path, pend_path = site / "history.json", site / "history.pending.json"
    nights, hnotes, clean = history.load(hist_path)
    notes += hnotes
    files, merged = [], False
    if pend_path.exists():
        pending, pnotes = history.load_pending(pend_path, date)
        notes += pnotes
        if not pend_path.exists():  # unreadable: set aside, both names are touched
            files += [pend_path.name, f"history.pending.broken-{date}.json"]
        elif clean and not pnotes:
            nights, merged = history.merge_pending(nights, pending), True
    nights, recorded = history.add_night(nights, tonight)
    if not recorded:
        notes.append(
            f"tonight is not recorded: a scheduled night already covers {date}"
        )

    asked = {}  # canary and expire ask about the same (build, job) pairs: once each

    def exists(build_id, job):
        if (build_id, job) not in asked:
            asked[build_id, job] = tc.physics_out_exists(base, build_id, job, get)
        return asked[build_id, job]

    canary_ok = history.canary([(b["id"], j) for j, b in builds.items()], exists)
    nights, xnotes = history.expire(nights, exists, canary_ok, cfg["page"]["dots"])
    notes += xnotes
    writes = {}
    if clean:
        writes[hist_path.name] = nights
        if merged:
            writes[pend_path.name] = None
    else:
        kept = history.load(pend_path)[0] if pend_path.exists() else []
        kept, _ = history.add_night(kept, tonight)
        writes[pend_path.name] = kept
        notes.append(
            f"history.json has problems: {len(kept)} nights kept aside in history.pending.json"
        )
    return nights, recorded, writes, list(dict.fromkeys(files + list(writes)))


def figure_url_for(base, figure_base):
    if figure_base is not None:
        return lambda job, build, name: (
            f"{figure_base.rstrip('/')}/{job}/figures/{name}"
        )
    return lambda job, build, name: (
        tc.artifact_url(
            base, build["type"], build["id"], f"physics-out/{job}/figures/{name}"
        )
        if build
        else None
    )


def push_refusal(a, builds, others=()):
    """Why --push is refused, or None: only a night whose every known build (the job builds,
    at least one, and the pipeline head and page writer when known) is on the default
    branch reaches the public page."""
    if a.figure_base is not None:
        return "--figure-base links to local files"
    if not builds:
        return "no job build known from TeamCity"
    if not all(b["default_branch"] for b in [*builds.values(), *others]):
        return "a known build is not on the default branch"
    if not os.environ.get("NIGHTLY_PAGE_TOKEN"):
        return "NIGHTLY_PAGE_TOKEN is not set"
    return None


def header_value(token):
    """The base64 of the credentials: what git sends, so what must never be printed."""
    return base64.b64encode(f"x-access-token:{token}".encode()).decode()


def scrub(text, token):
    return (
        text.replace(token, "***").replace(header_value(token), "***")
        if token
        else text
    )


def push_site(site, token, branch, run, files, date=None):
    """Commit the given files of the site clone (nothing if the tree is clean) and push them
    to origin/<branch>. The token and the commit identity reach git only through the
    environment: never argv, .git/config or a printed line."""
    date = date or datetime.date.today().isoformat()
    env = {k: v for k, v in os.environ.items() if k != "NIGHTLY_PAGE_TOKEN"} | {
        "GIT_TERMINAL_PROMPT": "0",
        "GIT_CONFIG_COUNT": "3",
        "GIT_CONFIG_KEY_0": "http.https://github.com/.extraheader",
        "GIT_CONFIG_VALUE_0": f"AUTHORIZATION: basic {header_value(token)}",
        "GIT_CONFIG_KEY_1": "user.name",
        "GIT_CONFIG_VALUE_1": "phare-nightly",
        "GIT_CONFIG_KEY_2": "user.email",
        "GIT_CONFIG_VALUE_2": "phare-nightly@nightly.invalid",
    }
    git = ["git", "-C", str(site)]

    def step(name, *args):
        r = run([*git, *args], env=env, capture_output=True, text=True)
        shown = scrub(f"{r.stdout}{r.stderr}".strip(), token)
        if r.returncode != 0 and name != "diff":  # main prints the one scrubbed line
            raise RuntimeError(
                f"git {name} failed (exit {r.returncode}): {shown[:200]}"
            )
        if shown and name != "diff":
            print(f"git {name}: {shown[:500]}")
        return r.returncode

    step("add", "add", "-A", "--", *files)
    if step("diff", "diff", "--cached", "--quiet") != 0:  # 1: something is staged
        step("commit", "commit", "-m", f"nightly {date}")
    step("push", "push", "origin", f"HEAD:{branch}")


def build(a, cfg, get, notes):
    """(index.html text, (job builds, head and writer builds), history writes, files
    touched, tonight's date)."""
    base = cfg["teamcity"]["url"]
    tz, tznotes = paris()
    notes += tznotes
    own = own_build_id(a)
    builds, manual, mine, head = jobs_tonight(base, cfg, own, a, get, notes)
    if not builds:
        notes.append("tonight's builds unknown")
    if mine is None:
        notes.append("the page writer's start is unknown: today's date is used")
    start = mine["start"].astimezone(tz) if mine else datetime.datetime.now(tz)
    parts = clean_parts(a, cfg)
    for part in parts:
        if not part.get("missing"):
            notes += part_notes(part, cfg, builds)
    run = {
        "sha": common_sha(builds, notes),
        "manual": manual,
        "time": start.strftime("%H:%M"),
        "builds": {
            j: {
                "id": b["id"],
                "type": b["type"],
                "start": b["start"].astimezone(tz).date().isoformat(),
            }
            for j, b in builds.items()
        },
    }
    tonight = history.night(start.date().isoformat(), run, records(parts, cfg))
    nights, recorded, writes, files = update_history(
        Path(a.site), tonight, builds, cfg, base, get, notes
    )
    html = page.render(
        nights,
        tonight,
        recorded,
        cfg,
        [p for p in parts if not p.get("missing")],
        notes,
        figure_url_for(base, a.figure_base),
    )
    return (
        html,
        (builds, [x for x in (head, mine) if x]),
        writes,
        files,
        tonight["date"],
    )


def write_site(site, html, writes):
    """index.html is replaced last, from a temporary file: a failure before that leaves
    yesterday's page as it was."""
    tmp = site / "index.html.tmp"
    tmp.write_text(html)
    for name, nights in writes.items():
        if nights is None:
            (site / name).unlink(missing_ok=True)
        else:
            history.dump(nights, site / name)
    os.replace(tmp, site / "index.html")


def main(argv=None, get=None, run=subprocess.run):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--jobs", required=True)
    ap.add_argument("--site", required=True)
    ap.add_argument("--config", default=None)
    ap.add_argument("--push", action="store_true")
    ap.add_argument("--figure-base", default=None)
    ap.add_argument("--own-build-id", type=int, default=None)
    ap.add_argument("--head-build-id", type=int, default=None)
    try:
        a = ap.parse_args(argv)
    except SystemExit:
        return 0  # never a gate: argparse has printed its usage
    token = os.environ.get("NIGHTLY_PAGE_TOKEN", "") if a.push else ""
    try:
        cfg, notes = decks.load(a.config), []
        html, known, writes, files, date = build(a, cfg, get, notes)
        write_site(Path(a.site), html, writes)
        if a.push:
            why = push_refusal(a, *known)
            if why:
                print(f"publish: push refused: {why}")
            else:
                touched = ["index.html", *files]
                push_site(a.site, token, cfg["page"]["branch"], run, touched, date)
    except Exception as e:  # noqa: BLE001  never a gate: one line, exit 0
        line = f"publish: page or push not done ({type(e).__name__}: {str(e)[:300]})"
        print(scrub(line.replace("\n", " "), token))
    return 0


if __name__ == "__main__":
    sys.exit(main())
