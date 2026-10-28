"""TeamCity helpers for the physics report: this build's identity, and guest (no token) reads."""

import json
import os
import re
import urllib.error
import urllib.request
from datetime import datetime

import history


def build_properties(path=None):
    path = path or os.environ.get("TEAMCITY_BUILD_PROPERTIES_FILE")
    props = {}
    if path and os.path.exists(path):
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if "=" in line and not line.startswith("#"):
                    k, v = line.rstrip("\n").split("=", 1)
                    props[k.strip()] = v.strip().replace("\\:", ":")
    return props


def http_get(url, timeout=30, headers=None):
    req = urllib.request.Request(
        url, headers={"Accept": "application/json"} | (headers or {})
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def get_json(url, get=None):
    return json.loads((get or http_get)(url))


def recent_builds(base, build_type_id, count, get=None):
    """Ids of the last finished nightly builds, newest first. defaultFilter:false alone would
    include personal, feature-branch, cancelled and never-started builds: excluded here."""
    url = (
        f"{base}/guestAuth/app/rest/builds?locator=buildType:(id:{build_type_id}),"
        "branch:(default:true),personal:false,canceled:false,failedToStart:false,"
        f"state:finished,count:{count},defaultFilter:false&fields=build(id)"
    )
    return [b["id"] for b in get_json(url, get).get("build", [])]


def artifact_url(base, build_type_id, build_id, path):
    return f"{base}/guestAuth/repository/download/{build_type_id}/{build_id}:id/{path}"


FIELDS = "id,buildTypeId,startDate,defaultBranch,triggered(type),revisions(revision(version))"
JOB_ID = re.compile(r"_virtual_([a-z0-9]+)_functional_")


def _build(b):
    """One TeamCity build as {"id", "type", "start", "default_branch", "triggered", "sha"},
    or None when an id, type or start is bad. A bad sha or trigger only gives None."""
    try:
        bid, btype = b["id"], b["buildTypeId"]
        if type(bid) is not int or not history.BUILD_TYPE.fullmatch(btype):
            return None
        start = datetime.strptime(b["startDate"], "%Y%m%dT%H%M%S%z")
    except (KeyError, TypeError, ValueError):
        return None
    sha = None
    try:
        version = b["revisions"]["revision"][0]["version"]
        if isinstance(version, str) and history.SHA.fullmatch(version):
            sha = version
    except (KeyError, IndexError, TypeError):
        pass
    trigger = b.get("triggered")
    triggered = trigger.get("type") if isinstance(trigger, dict) else None
    return {
        "id": bid,
        "type": btype,
        "start": start,
        "default_branch": b.get("defaultBranch") is True,
        "triggered": triggered if isinstance(triggered, str) else None,
        "sha": sha,
    }


def _build_list(url, get):
    """The "build" list of a TeamCity answer; anything else is no builds."""
    data = get_json(url, get)
    items = data.get("build") if isinstance(data, dict) else None
    return [b for b in items if isinstance(b, dict)] if isinstance(items, list) else []


def upstream_builds(base, build_id, get=None):
    """The builds the given build depends on (the job builds of a pipeline head)."""
    url = (
        f"{base}/guestAuth/app/rest/builds?locator=snapshotDependency:(to:(id:{int(build_id)})),"
        f"defaultFilter:false&fields=build({FIELDS})"
    )
    found = (_build(b) for b in _build_list(url, get))
    return [b for b in found if b]


def artifact_sources(base, build_id, get=None):
    """The builds whose files the given build downloaded (a page writer's job builds)."""
    url = (
        f"{base}/guestAuth/app/rest/builds/id:{int(build_id)}"
        f"?fields=artifact-dependencies(build({FIELDS}))"
    )
    data = get_json(url, get)
    deps = data.get("artifact-dependencies") if isinstance(data, dict) else None
    items = deps.get("build") if isinstance(deps, dict) else None
    if not isinstance(items, list):
        return []
    found = (_build(b) for b in items if isinstance(b, dict))
    return [b for b in found if b]


def last_finished(base, build_type_id, get=None):
    """The last finished build of a build type on the default branch, or None."""
    if not isinstance(build_type_id, str) or not history.BUILD_TYPE.fullmatch(
        build_type_id
    ):
        return None
    url = (
        f"{base}/guestAuth/app/rest/builds?locator=buildType:(id:{build_type_id}),"
        "branch:(default:true),state:finished,defaultFilter:false,count:1"
        f"&fields=build({FIELDS})"
    )
    found = [b for b in map(_build, _build_list(url, get)) if b]
    return found[0] if found else None


def build_info(base, build_id, get=None):
    url = f"{base}/guestAuth/app/rest/builds/id:{int(build_id)}?fields={FIELDS}"
    return _build(get_json(url, get))


def head_of(base, own_build_id, get=None):
    """(newest pipeline head id, other head ids) of the builds depending on this one."""
    url = (
        f"{base}/guestAuth/app/rest/builds?locator=snapshotDependency:(from:(id:{int(own_build_id)})),"
        "defaultFilter:false&fields=build(id,buildTypeId,startDate)"
    )
    ids = sorted(
        b["id"]
        for b in _build_list(url, get)
        if type(b.get("id")) is int
        and isinstance(b.get("buildTypeId"), str)
        and history.BUILD_TYPE.fullmatch(b["buildTypeId"])
        and "_virtual_" not in b["buildTypeId"]
    )
    return (ids[-1], ids[:-1][::-1]) if ids else (None, [])


def job_of(type_id, cfg):
    """The decks.yml job of a job build type id, or None."""
    m = JOB_ID.search(type_id)
    return m.group(1) if m and m.group(1) in cfg["jobs"] else None


def physics_out_exists(base, build_id, job, get=None):
    """True if the build still lists physics-out/<job>, False on HTTP 404, None otherwise."""
    url = (
        f"{base}/guestAuth/app/rest/builds/id:{int(build_id)}"
        f"/artifacts/children/physics-out/{job}"
    )
    try:
        (get or http_get)(url, timeout=10)
        return True
    except urllib.error.HTTPError as e:
        return False if e.code == 404 else None
    except (urllib.error.URLError, OSError):
        return None
