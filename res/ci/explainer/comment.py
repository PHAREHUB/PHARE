#!/usr/bin/env python3
"""One PR comment, about the PR's newest red commit: the reds the agent analysed, then the ones
the rules put on the CI (one line each).

One comment per PR, edited on each new commit: the hidden marker lets ci.py find it again.
"""

import re

# code shown as written: only the plain forms, ``` blocks opening and closing a line and single-backtick
# spans not touching a backslash or another backtick. Anything GitHub might read differently
# (\`, ``` mid-line, ~~~) is treated as text and escaped: at worst a '<' shows as &lt;.
CODE = re.compile(
    r"^ {0,3}```[^`\n]*\n.*?^ {0,3}```[ \t]*$|(?<![\\`])`[^`\n]+`(?!`)", re.S | re.M
)
MAX_REPORT = 1500  # characters per agent report, after cleaning
SKIPPED = object()  # report of a failure the agent was not run on


def marker(pr):
    return f"<!-- ci-explainer:{pr} -->"


# the red builds a comment reports, for ci.done
BUILDS = re.compile(r"<!-- builds:([\d,]*) -->")


def builds_tag(ids):
    return f"<!-- builds:{','.join(map(str, sorted(ids)))} -->"


def clean(report):
    """Agent text made safe to post: no @-mentions, no links, no HTML, bounded length.
    '<' is escaped, not removed, outside code: `std::vector<double>` and #include <cstdint> survive."""
    parts, end = [], 0
    for code in CODE.finditer(report):
        parts += [report[end : code.start()].replace("<", "&lt;"), code[0]]
        end = code.end()
    text = "".join(parts) + report[end:].replace("<", "&lt;")
    text = re.sub(
        r"\[([^\]]*)\]\([^)]*\)", r"\1", text
    )  # markdown links: keep the words
    # reference links: without their "[r]: url" line, [words][r] is plain text (and a[i][j] stays)
    text = re.sub(r"^ {0,3}\[[^\]]+\]:.*$", "", text, flags=re.M)
    # GitHub links www. addresses too
    text = re.sub(r"(https?://|www\.)\S+", "(link removed)", text)
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
            builds_tag([]),
            f"**CI on `{sha[:7]}`: no red build.** The reds explained here before are gone.",
        ]
    )
