#!/usr/bin/env python3
"""Whose fault is a red build: the CI's (one line, no AI) or the PR's (excerpt, then the agent).

Rules come from a week of real reds (2026-10-19 to 10-26, 166 PR builds): registry outage,
60-min timeouts, compile errors, test failures. "Tests passed: N" with ctest exit code 8 is a real
failure: TeamCity counts only googletest results, so a failed Python test doesn't show in N.

usage: triage.py BUILD_ID      prints the verdict for one build
"""

import re
import sys
from dataclasses import dataclass

# a TeamCity log line: "[20:36:35]W:\t [Step 3/5] text"; nested blocks add tabs and "[name] "
TC_PREFIX = re.compile(r"^\[\d\d:\d\d:\d\d\][A-Z ]?:\s*(\[Step \d+/\d+\] ?)?")
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def clean(log):
    """Log text without ANSI colours and TeamCity timestamps / step prefixes."""
    return [TC_PREFIX.sub("", ANSI.sub("", line)) for line in log.splitlines()]


@dataclass
class Verdict:
    kind: str  # canceled upstream test compile infra timeout false_red unknown
    owner: str  # "ci", "pr" or "unknown"
    reason: str  # one line, for people
    tests: tuple = ()  # names of the failed tests, from the log


INFRA = (
    (
        re.compile(r"\[pull [^\]]+\] Process exited with code [1-9]"),
        "the build image could not be downloaded from the registry",
    ),
    (re.compile(r"no space left on device", re.I), "the build machine ran out of disk"),
    (
        re.compile(
            r"connection to the build agent was lost|agent .* (lost|disconnected|unregistered)",
            re.I,
        ),
        "the build machine disconnected",
    ),
    (
        re.compile(r"Error while applying patch|cannot find commit", re.I),
        "TeamCity could not check out the commit (it moved, or was force-pushed)",
    ),
    (
        re.compile(r"Killed signal terminated program|virtual memory exhausted", re.I),
        "the compiler ran out of memory",
    ),
)

COMPILE_ERROR = re.compile(r"^\S+?:\d+:\d+: (fatal )?error: ")
TESTS_FAILED = re.compile(r"Tests failed: \d+")
# ctest: "   29 - py3_init-densities_mpi_n_2 (Failed)" after "The following tests FAILED:"
CTEST_FAILED = re.compile(
    r"^\s*\d+ - (\S+) \((Failed|Timeout|SEGFAULT|Subprocess aborted|Child aborted|Exception|Not Run)\)"
)
# phlop: "mpirun -n 4 ... python3 -um test_divb_refinement Test.test_x FAILED: exited 1"
PHLOP_FAILED = re.compile(r"python3 -u\S* (.+?) FAILED: exited")


def failed_tests(lines):
    names = []
    for line in lines:
        match = CTEST_FAILED.match(line) or PHLOP_FAILED.search(line)
        if match and match.group(1) not in names:
            names.append(match.group(1))
    return tuple(names)


def canceled(build):
    return build["status"] == "UNKNOWN" or bool(
        re.search(r"\bCanceled\b", build["status_text"])
    )


def classify(build, lines):
    """build: dict from tc.flatten; lines: clean(log).

    Order matters: a build can carry a harmless infra warning and still fail on a real test, so
    what the PR broke is looked for before what the CI broke."""
    text = build["status_text"]
    if canceled(build):
        return Verdict(
            "canceled", "ci", "canceled (usually a newer push replaced this build)"
        )
    if "Snapshot dependency failed" in text:
        return Verdict(
            "upstream", "ci", "an earlier job of this chain failed; see that one"
        )
    tests = failed_tests(lines)
    if tests or TESTS_FAILED.search(text):
        return Verdict("test", "pr", "test failure", tests)
    errors = [line for line in lines if COMPILE_ERROR.match(line)]
    if errors:
        if all("/subprojects/" in line for line in errors):
            return Verdict(
                "compile",
                "unknown",
                "compile error inside a third-party library (SAMRAI, caliper, ...)",
            )
        return Verdict("compile", "pr", "compile error")
    for pattern, reason in INFRA:
        if any(pattern.search(line) for line in lines):
            return Verdict("infra", "ci", reason)
    if re.search(r"Execution timeout", text) or any(
        "Execution timeout" in line for line in lines
    ):
        return Verdict(
            "timeout",
            "ci",
            "hit the build time limit; re-run, and tell the CI maintainers if it repeats",
        )
    # a step's exit code means something failed: the rules just didn't find what (agent, tail excerpt)
    if text.startswith("Tests passed:") and "exit code" not in text.lower():
        return Verdict(
            "false_red",
            "ci",
            "all tests passed, but TeamCity marked the build red while reading test results",
        )
    return Verdict("unknown", "unknown", "red for a reason the rules don't know")


if __name__ == "__main__":
    import tc

    if len(sys.argv) != 2:
        sys.exit(__doc__)
    build = tc.build(sys.argv[1])
    print(classify(build, clean(tc.log(sys.argv[1]))))
