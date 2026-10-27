#!/usr/bin/env python3
"""Read-only TeamCity access through guest (anonymous) REST: no token, public data only.

usage: tc.py BUILD_ID      prints the build summary and its failed tests, for a quick look
"""

import json
import re
import sys
import urllib.parse
import urllib.request

SERVER = "https://hephaistos.lpp.polytechnique.fr/teamcity"
REST = f"{SERVER}/guestAuth/app/rest"

# PR configurations, production and shadow: the builds whose reds we explain
PR_CONFIGS = (
    "Phare_Phare_BuildGithubPr",
    "Phare_Phare_BuildGithubPrClang",
    "Phare_Phare_BuildGithubPrGccSs",
    "Phare_Phare_BuildGithubPrClangSs",
    "Phare_Next_PrChecksV1_virtual_gcc_V__1",
    "Phare_Next_PrChecksV1_virtual_clang_asan_V__1",
)

BUILD_FIELDS = "id,buildTypeId,status,state,branchName,startDate,finishDate,statusText,webUrl,revisions(revision(version))"


def get(path, accept="application/json", timeout=60):
    request = urllib.request.Request(f"{REST}/{path}", headers={"Accept": accept})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        body = response.read().decode(errors="replace")
    return json.loads(body) if accept == "application/json" else body


def flatten(build):
    """REST build -> flat dict with the PR number and the commit that was built."""
    revisions = build.get("revisions", {}).get("revision", [])
    branch = build.get("branchName") or ""
    pr = re.fullmatch(r"pull/(\d+)", branch)
    return {
        "id": build["id"],
        "config": build["buildTypeId"],
        "status": build.get("status"),
        "status_text": build.get("statusText") or "",
        "branch": branch,
        "pr": int(pr.group(1)) if pr else None,
        "sha": revisions[0]["version"] if revisions else None,
        "start": build.get("startDate"),
        "finish": build.get("finishDate"),
        "url": build.get("webUrl"),
    }


def red_builds(since, configs=PR_CONFIGS, count=1000):
    """Finished, not-successful builds of 'configs' since 'since' (TeamCity date, e.g. 20261019T000000+0000)."""
    builds = []
    for config in configs:
        locator = f"buildType:{config},state:finished,status:FAILURE,defaultFilter:false,branch:default:any,sinceDate:{since},count:{count}"
        answer = get(
            f"builds?locator={urllib.parse.quote(locator, safe=':,()')}&fields=build({BUILD_FIELDS})"
        )
        builds += [flatten(b) for b in answer.get("build", [])]
    return sorted(builds, key=lambda b: b["id"])


def build(build_id):
    return flatten(get(f"builds/id:{build_id}?fields={BUILD_FIELDS}"))


def builds_of(sha):
    """Every build of commit 'sha' (full sha), any configuration, any state."""
    locator = f"revision:{sha},defaultFilter:false,branch:default:any"
    answer = get(
        f"builds?locator={urllib.parse.quote(locator, safe=':,()')}&fields=build({BUILD_FIELDS})"
    )
    return [{**flatten(b), "state": b.get("state")} for b in answer.get("build", [])]


def log(build_id):
    """The full build log as text."""
    request = urllib.request.Request(
        f"{SERVER}/guestAuth/downloadBuildLog.html?buildId={build_id}"
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read().decode(errors="replace")


def failed_tests(build_id, count=20):
    """[(test name, details)]: details hold the assertion message and the test's stdout."""
    locator = f"build:(id:{build_id}),status:FAILURE,count:{count}"
    answer = get(
        f"testOccurrences?locator={urllib.parse.quote(locator, safe=':,()')}&fields=testOccurrence(name,details)"
    )
    return [
        (t["name"], t.get("details") or "") for t in answer.get("testOccurrence", [])
    ]


def test_history(name, count=100):
    """[(status, branch)] of the last runs of test 'name', any configuration, any branch."""
    locator = f"test:(name:{name}),count:{count}"
    answer = get(
        f"testOccurrences?locator={urllib.parse.quote(locator, safe=':,()')}&fields=testOccurrence(status,build(branchName))"
    )
    return [
        (t["status"], t["build"].get("branchName"))
        for t in answer.get("testOccurrence", [])
    ]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    summary = build(sys.argv[1])
    print(json.dumps(summary, indent=2))
    for name, details in failed_tests(sys.argv[1]):
        print(f"--- {name}\n{details[:500]}")
