#!/usr/bin/env python3
"""The red explainer in GitHub Actions: one PR commit, once all its TeamCity builds are done.

TeamCity reports each PR build to GitHub as a commit status ("gh_pr_gcc (Phare)", ...): pending
when queued, then success or failure. The workflow wakes on each status; this script acts only
when none of TeamCity's statuses on the commit is pending any more.

usage: ci.py check SHA                    prints action=red|green|none, pr=N, builds=ID,ID
       ci.py explain PR SHA BUILDS OUT    writes OUT/comment.md (BUILDS: red build ids, comma-separated)
       ci.py post PR SHA FILE             creates or edits the explainer's comment, if SHA is still the head
       ci.py green PR SHA                 edits it to say the reds are gone

GitHub is read with GITHUB_TOKEN (post and green need pull-requests: write; without a token,
reads only: a local dry run), TeamCity as guest.
EXPLAIN_LABEL: explain only PRs carrying this label (trial); unset or empty: every PR.
"""

import json
import os
import re
import sys
import urllib.request
from pathlib import Path

import agent
import comment
import explain
import tc
import triage

REPO = os.environ.get("GITHUB_REPOSITORY", "PHAREHUB/PHARE")
API = f"https://api.github.com/repos/{REPO}"
TEAMCITY = " (Phare)"  # suffix of the status contexts TeamCity publishes
BOT = "github-actions[bot]"  # author of our comments: a marker in anyone else's is ignored
MAX_FAILURES = 4  # distinct failures explained per commit: up to 10 minutes each


def api(path, method="GET", body=None):
    headers = {
        "Accept": "application/vnd.github+json",
        "Content-Type": "application/json",
    }
    if os.environ.get(
        "GITHUB_TOKEN"
    ):  # without: public reads only, for a local dry run
        headers["Authorization"] = f"Bearer {os.environ['GITHUB_TOKEN']}"
    request = urllib.request.Request(
        f"{API}/{path}",
        method=method,
        data=json.dumps(body).encode() if body is not None else None,
        headers=headers,
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def build_id(url):
    """TeamCity build id at the end of a status link (.../build/100788 or .../Phare_Phare_BuildGithubPr/100788)."""
    found = re.search(r"/(\d+)/?$", url or "")
    return int(found.group(1)) if found else None


def decide(statuses):
    """Latest status per context -> ('red' | 'green' | 'none', red build ids).
    'none' while a TeamCity build is still queued or running, or if TeamCity reported nothing."""
    ours = [s for s in statuses if s["context"].endswith(TEAMCITY)]
    if not ours or any(s["state"] == "pending" for s in ours):
        return "none", []
    red = sorted(
        build_id(s["target_url"]) for s in ours if s["state"] in ("failure", "error")
    )
    return ("red", red) if red else ("green", [])


def our_comment(pr):
    """The explainer's comment on the PR, or None."""
    page = 1
    while True:
        comments = api(f"issues/{pr}/comments?per_page=100&page={page}")
        for c in comments:
            if c["user"]["login"] == BOT and c["body"].startswith(comment.marker(pr)):
                return c
        if len(comments) < 100:
            return None
        page += 1


def done(body, sha, action, red=()):
    """Does the comment already say this about this commit, with exactly these red builds?
    (TeamCity sends some statuses twice; a re-run build gets a new id, a passing re-run drops one.)"""
    head = f"**CI on `{sha[:7]}`: "
    said_green = (head + "no red build") in body
    tag = comment.BUILDS.search(body)
    reported = {int(i) for i in tag[1].split(",") if i} if tag else None
    return head in body and said_green == (action == "green") and reported == set(red)


def current(pr, sha):
    """Is sha still the head of the open PR, with the trial label? Asked again before posting:
    an explanation takes minutes, and a newer commit's comment must not be overwritten."""
    pull = api(f"pulls/{pr}")
    label = os.environ.get("EXPLAIN_LABEL", "")
    return not (
        pull["state"] != "open"
        or pull["head"]["sha"] != sha  # pushed again since: the new commit will come
        or (label and label not in {lab["name"] for lab in pull["labels"]})
    )


def wanted(pr, sha, action, red):
    if not current(pr, sha):
        return False
    old = our_comment(pr)
    if old is None:
        return action == "red"  # green with no comment: nothing to correct
    return not done(old["body"], sha, action, red)


def check(sha):
    statuses = api(f"commits/{sha}/status?per_page=100")["statuses"]
    action, red = decide(statuses)
    pr = None
    if action != "none":
        # the PR number comes from TeamCity (branch pull/N): GitHub can't map a fork's commit to its PR
        any_build = (
            red[0]
            if red
            else next(
                build_id(s["target_url"])
                for s in statuses
                if s["context"].endswith(TEAMCITY)
            )
        )
        pr = tc.build(any_build)["pr"]
    if action == "red" and all(triage.canceled(tc.build(i)) for i in red):
        action = (
            "none"  # canceled on purpose (newer push, CI-only PR): nothing to explain
        )
    if action == "none" or pr is None or not wanted(pr, sha, action, red):
        action = "none"
    print(f"action={action}\npr={pr or ''}\nbuilds={','.join(map(str, red))}")


def explain_commit(pr, sha, ids, out):
    builds = [tc.build(i) for i in ids]
    _, groups, others = explain.sort(builds, tc.log)
    try:
        siblings = tc.builds_of(sha)
    except OSError:  # a clue, not a must: explain without it
        siblings = []
    explained, others = [], others.get((pr, sha), [])
    for n, members in enumerate(
        sorted(groups.values(), key=lambda members: members[0][0]["id"])
    ):
        report = (
            explain.explain(
                members,
                ".",
                "origin",
                out / "work" / str(members[0][0]["id"]),
                siblings=siblings,
            )
            if n < MAX_FAILURES
            else comment.SKIPPED
        )
        explained.append(([b for b, _, _ in members], members[0][1], report))
    if not others and not any(isinstance(r, str) for _, _, r in explained):
        print(
            "no usable explanation, nothing put on the CI: no comment"
        )  # e.g. a wrong model name
        return
    text = comment.render(pr, sha, explained, others, explain.signature(agent.MODEL))
    (out / "comment.md").write_text(text)


def post(pr, sha, body):
    if not current(pr, sha):
        print(
            f"{sha[:7]} is no longer the head of open PR {pr} with its label: nothing posted"
        )
        return
    old = our_comment(pr)
    if old is None:
        api(f"issues/{pr}/comments", "POST", {"body": body})
    else:
        api(f"issues/comments/{old['id']}", "PATCH", {"body": body})


if __name__ == "__main__":
    match sys.argv[1:]:
        case ["check", sha]:
            check(sha)
        case ["explain", pr, sha, ids, out]:
            Path(out).mkdir(parents=True, exist_ok=True)
            explain_commit(int(pr), sha, [int(i) for i in ids.split(",")], Path(out))
        case ["post", pr, sha, path]:
            if Path(path).exists():  # explain may have had nothing to say
                post(int(pr), sha, Path(path).read_text())
        case ["green", pr, sha]:
            post(int(pr), sha, comment.render_green(int(pr), sha))
        case _:
            sys.exit(__doc__)
