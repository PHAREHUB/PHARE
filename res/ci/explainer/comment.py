#!/usr/bin/env python3
"""One PR comment, about the PR's newest red commit: the reds it caused (agent reports), then the
others (one line each).

One comment per PR, edited on each new commit: the hidden marker lets ci.py find it again.
"""

import re

MAX_REPORT = 1500  # characters per agent report, after cleaning
SKIPPED = object()  # report of a failure the agent was not run on


def marker(pr):
    return f"<!-- ci-explainer:{pr} -->"


def clean(report):
    """Agent text made safe to post: no @-mentions, no links, no HTML, bounded length."""
    text = re.sub(r"<[^>\n]+>", "", report)  # HTML tags, comments
    text = re.sub(
        r"\[([^\]]*)\]\([^)]*\)", r"\1", text
    )  # markdown links: keep the words
    text = re.sub(r"https?://\S+", "(link removed)", text)
    text = re.sub(r"@(?=[\w-])", "@\u200b", text)  # zero-width space: no GitHub mention
    text = "\n".join(
        line for line in text.splitlines() if not line.startswith("Owner:")
    )
    return text if len(text) <= MAX_REPORT else text[:MAX_REPORT] + " [...]"


def layout(report):
    """'What: x' -> '**What** x'; the fix prompt stays a quote block."""
    text = re.sub(r"^(What|Why):\s*", r"**\1** ", report, flags=re.M)
    return re.sub(
        r"^Fix prompt:\s*$",
        "**Fix**, to paste into your coding agent:",
        text,
        flags=re.M,
    )


def render(pr, sha, explained, others, signature):
    """explained: [(builds, verdict, report)], builds of one distinct failure; report is the
    agent's text, None if unusable, SKIPPED if not run. others: [(build, verdict)] reds this PR
    didn't cause."""
    reds = sum(len(builds) for builds, _, _ in explained) + len(others)
    lines = [
        marker(pr),
        f"**CI on `{sha[:7]}`: {reds} red**, {reds - len(others)} from this PR, {len(others)} not.",
        "",
    ]
    for builds, verdict, report in explained:
        where = ", ".join(
            f"{b['config'].replace('Phare_Phare_Build', '').replace('Phare_Next_PrChecksV1_virtual_', 'next ')} ({b['id']})"
            for b in builds
        )
        lines += [
            f"#### {verdict.reason}: {where}",
            layout(clean(report))
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
            f"**CI on `{sha[:7]}`: no red build.** The reds explained here before are gone.",
        ]
    )
