#!/usr/bin/env python3
"""The few log lines that explain a red: what the agent reads instead of a 4000-line log.

usage: excerpt.py BUILD_ID      prints the excerpt for one build
"""

import re
import sys

from triage import COMPILE_ERROR

MAX_CHARS = 8000  # per excerpt: free models have small contexts, and people skim
NINJA_FAILED = re.compile(r"^FAILED: ")
CTEST_LINE = re.compile(r"^\s*\d+/\d+ Test +#\d+: (\S+) ")
CHECKOUT = re.compile(
    r"/opt/buildagent/work/[0-9a-f]+/"
)  # agent checkout dir -> repo-relative paths
NOISE = re.compile(
    r"^(\[\d+/\d+\] |ccache |\[ RUN      \]|\[       OK \]|\[----------\] \d+ tests? from .* \(\d+ ms total\)$)"
)


def tidy(lines):
    """Repo-relative paths; no compile commands, no googletest lines of tests that passed."""
    return [CHECKOUT.sub("", line) for line in lines if not NOISE.match(line)]


def cap(lines, head=60, tail=100):
    if len(lines) > head + tail:
        lines = (
            lines[:head]
            + [f"[... {len(lines) - head - tail} lines cut ...]"]
            + lines[-tail:]
        )
    text = "\n".join(lines)
    return (
        text
        if len(text) <= MAX_CHARS
        else text[: MAX_CHARS // 2] + "\n[... cut ...]\n" + text[-MAX_CHARS // 2 :]
    )


def compile_excerpt(lines):
    """Every distinct error line, then the first failing compile command with what follows it."""
    errors = list(dict.fromkeys(line for line in lines if COMPILE_ERROR.match(line)))
    errors.sort(key=lambda line: "/subprojects/" in line)  # PHARE's own errors first
    first = next(i for i, line in enumerate(lines) if COMPILE_ERROR.match(line))
    start = max(
        (i for i in range(first) if NINJA_FAILED.match(lines[i])), default=first
    )
    # up to the next failing file: the same error repeats once per compiled module
    end = next(
        (
            i
            for i in range(first, first + 40)
            if i < len(lines) and NINJA_FAILED.match(lines[i])
        ),
        first + 40,
    )
    return cap(
        tidy(["Errors:"] + errors[:10] + ["", "First failure:"] + lines[start:end])
    )


def test_excerpt(lines, name):
    """Output ctest printed for failed test 'name' (--output-on-failure), or the phlop stdout/stderr."""
    for i, line in enumerate(lines):
        match = CTEST_LINE.match(line)
        if (
            match
            and match.group(1) == name
            and re.search(r"\*\*\*(Failed|Timeout|Exception)|SEGFAULT|aborted", line)
        ):
            end = next(
                (j for j in range(i + 1, len(lines)) if CTEST_LINE.match(lines[j])),
                len(lines),
            )
            return cap(tidy(lines[i:end]))
        if name in line and "FAILED: exited" in line:
            return cap(tidy(lines[i : i + 60]), head=60, tail=0)
    return f"(no output found in the log for {name})"


def tail_excerpt(lines, n=80):
    """The lines before the first failing step: for reds the rules don't know."""
    end = next(
        (
            i
            for i, line in enumerate(lines)
            if re.search(r"Process exited with code [1-9]", line)
        ),
        len(lines),
    )
    return cap(tidy(lines[max(0, end - n) : end + 1]))


def excerpt(verdict, lines, max_tests=3):
    if verdict.kind == "compile":
        return compile_excerpt(lines)
    if verdict.kind == "test" and verdict.tests:
        return "\n\n".join(
            f"=== {name}\n{test_excerpt(lines, name)}"
            for name in verdict.tests[:max_tests]
        )
    return tail_excerpt(lines)


if __name__ == "__main__":
    import tc
    import triage

    if len(sys.argv) != 2:
        sys.exit(__doc__)
    lines = triage.clean(tc.log(sys.argv[1]))
    print(excerpt(triage.classify(tc.build(sys.argv[1]), lines), lines))
