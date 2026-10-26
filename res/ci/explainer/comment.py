#!/usr/bin/env python3
"""One PR comment, about the PR's newest red commit: the reds the agent analysed, then the ones
the rules put on the CI (one line each).

One comment per PR, edited on each new commit: the hidden marker lets ci.py find it again.
"""

import re
import textwrap

MAX_REPORT = 1500  # characters per agent report
WIDTH = 100  # a code block doesn't wrap on GitHub
SKIPPED = object()  # report of a failure the agent was not run on


def marker(pr):
    return f"<!-- ci-explainer:{pr} -->"


# the red builds a comment reports, for ci.done
BUILDS = re.compile(r"<!-- builds:([\d,]*) -->")


def builds_tag(ids):
    return f"<!-- builds:{','.join(map(str, sorted(ids)))} -->"


def as_code(report):
    """The agent's text as one code block, the only safe way to post it: GitHub renders nothing
    inside (no link, HTML, mention or entity), whatever the text holds; escaping it as markdown
    failed on one trick after another. The fence is longer than any backtick run in the text, so
    the text can't close it."""
    lines = []
    for line in report.splitlines():
        if line.startswith("Owner:"):  # for replay's grading, not for people
            continue
        line = re.sub(
            r"^> ?", "", line
        )  # the fix prompt's quote marks: noise in a block
        lines += textwrap.wrap(line, WIDTH, break_on_hyphens=False) or [""]
    text = "\n".join(lines)
    if len(text) > MAX_REPORT:
        text = text[:MAX_REPORT] + " [...]"
    longest = max((len(run) for run in re.findall(r"`+", text)), default=0)
    fence = "`" * max(3, longest + 1)
    return f"{fence}text\n{text}\n{fence}"


def render(pr, sha, explained, others, signature):
    """explained: [(builds, verdict, report)], builds of one distinct failure; report is the
    agent's text, None if unusable, SKIPPED if not run. others: [(build, verdict)] reds the rules
    put on the CI."""
    ids = [b["id"] for builds, _, _ in explained for b in builds] + [
        b["id"] for b, _ in others
    ]
    told = sum(
        len(builds) for builds, _, report in explained if isinstance(report, str)
    )
    counts = (
        (told, "explained below"),
        (len(ids) - told - len(others), "not explained"),
        (len(others), "put on the CI by the rules"),
    )
    lines = [
        marker(pr),
        builds_tag(ids),
        f"**CI on `{sha[:7]}`: {len(ids)} red**: "
        + ", ".join(f"{n} {what}" for n, what in counts if n)
        + ".",
        "",
    ]
    for builds, verdict, report in explained:
        where = ", ".join(
            f"{b['config'].replace('Phare_Phare_Build', '').replace('Phare_Next_PrChecksV1_virtual_', 'next ')} ({b['id']})"
            for b in builds
        )
        lines += [
            f"#### {verdict.reason}: {where}",
            as_code(report)
            if isinstance(report, str)
            else "_(not explained: too many distinct failures on this commit)_"
            if report is SKIPPED
            else "_(no explanation: the agent's answer was unusable)_",
            "",
        ]
    if others:
        lines += (
            ["Not this PR:"]
            + [f"- {b['config']} ({b['id']}): {v.reason}" for b, v in others]
            + [""]
        )
    lines += [f"<sub>{signature} Advisory, never blocks. Rate with 👍 / 👎.</sub>"]
    return "\n".join(lines)


def render_green(pr, sha):
    """Replaces the comment once a later commit has no red build."""
    return "\n".join(
        [
            marker(pr),
            builds_tag([]),
            f"**CI on `{sha[:7]}`: no red build.** The reds explained here before are gone.",
        ]
    )
