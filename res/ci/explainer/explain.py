#!/usr/bin/env python3
"""Shared by replay.py (past reds, offline) and ci.py (one commit, in GitHub Actions).

sort() turns red builds into distinct failures; explain() runs the agent on one of them, in a
folder holding the excerpt, the PR's diff and the code of the built commit.
"""

import os
import re
import subprocess
from collections import defaultdict

import agent
import excerpt
import triage


def git(repo, *args, **kwargs):
    return subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
        **kwargs,
    ).stdout


def have(repo, sha):
    return (
        subprocess.run(
            ["git", "-C", str(repo), "cat-file", "-e", f"{sha}^{{commit}}"],
            capture_output=True,
            check=False,
        ).returncode
        == 0
    )


def fetch(repo, remote, pr, sha):
    """Make 'sha' available locally: the PR head first, the bare commit if the PR moved since."""
    for ref in (f"pull/{pr}/head", sha):
        if have(repo, sha):
            return True
        subprocess.run(
            ["git", "-C", str(repo), "fetch", "-q", remote, ref],
            capture_output=True,
            check=False,
        )
    return have(repo, sha)


def diff(repo, remote, sha):
    base = git(repo, "merge-base", f"{remote}/master", sha).strip()
    return git(repo, "diff", "--stat", base, sha) + "\n" + git(repo, "diff", base, sha)


def snapshot(repo, sha):
    def write(dest):
        dest.mkdir(parents=True)
        archive = subprocess.run(
            ["git", "-C", str(repo), "archive", sha], check=True, capture_output=True
        ).stdout
        subprocess.run(["tar", "-x", "-C", str(dest)], input=archive, check=True)
        # a PR's symlink would let the agent read outside its folder (the sandbox checks paths, not
        # link targets); PHARE has none of its own
        for root, dirs, files in os.walk(dest):
            for name in dirs + files:
                if os.path.islink(path := os.path.join(root, name)):
                    os.unlink(path)

    return write


def failure_key(verdict, text):
    """Same failure on several configurations -> same key."""
    if (
        verdict.kind == "compile"
    ):  # file:line only: gcc and clang word the same error differently
        first = next(
            (line for line in text.splitlines() if triage.COMPILE_ERROR.match(line)), ""
        )
        path, line = (first.split(":") + ["", ""])[:2]
        return ("compile", re.sub(r"^(\.\.?/)+", "", path) + ":" + line)
    # no test name: only the same excerpt is the same failure
    return (verdict.kind, verdict.tests or text)


def sort(builds, log_of):
    """log_of: build id -> log text. Returns (rows, groups, others):
    rows    every build with its verdict, for summaries
    groups  (pr, sha, failure key) -> [(build, verdict, excerpt)]: reds the PR may have caused
    others  (pr, sha) -> [(build, verdict)]: reds of the CI, one line each, no agent"""
    rows, groups, others = [], defaultdict(list), defaultdict(list)
    for build in builds:
        lines = triage.clean(log_of(build["id"]))
        verdict = triage.classify(build, lines)
        rows.append(
            {
                **build,
                "kind": verdict.kind,
                "owner": verdict.owner,
                "tests": verdict.tests,
            }
        )
        if build["pr"] is None:
            continue
        if verdict.owner == "ci":
            others[(build["pr"], build["sha"])].append((build, verdict))
            continue
        text = excerpt.excerpt(verdict, lines)
        groups[(build["pr"], build["sha"], failure_key(verdict, text))].append(
            (build, verdict, text)
        )
    return rows, groups, others


def explain(members, repo, remote, workdir, model=agent.MODEL, siblings=()):
    """One distinct failure ([(build, verdict, excerpt)]) -> the agent's report, or None.
    siblings: the commit's other builds (tc.builds_of). The agent's folder stays in workdir,
    with its raw output in raw.txt."""
    first, verdict, text = members[0]
    if not fetch(repo, remote, first["pr"], first["sha"]):
        return None
    agent.prepare(
        workdir,
        agent.context(first, verdict, text, diff(repo, remote, first["sha"]), siblings),
        snapshot(repo, first["sha"]),
    )
    report, raw = agent.run(workdir, model)
    (workdir / "raw.txt").write_text(raw)
    return report


def opencode_version():
    try:
        return subprocess.run(
            ["opencode", "--version"],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "?"


def signature(model):
    return f"Written by {model} via OpenCode {opencode_version()}."
