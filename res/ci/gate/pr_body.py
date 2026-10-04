#!/usr/bin/env python3
"""Check that the PR description fills the required sections of .github/pull_request_template.md.

A section is filled when something other than whitespace is left once the template's
<!-- comments --> are removed. "Issue" may be a written justification rather than #123, and
"closes #123" (fixes, resolves, relates to, ...) anywhere in the description also counts as one.
A missing or empty Tests section is a warning.

usage: pr_body.py BODY_FILE     (the PR description, e.g. from github.event.pull_request.body)
exit:  0 the issue is given, 1 otherwise
"""

import re
import sys

from gitdiff import Finding, report

CHECK = "pr-body"
REQUIRED = {"Issue": "error", "Tests": "warning"}
ISSUE_LINK = re.compile(r"\b(close[sd]?|fix(e[sd])?|resolve[sd]?|relates?(\s+to)?)\b:?\s+(\S+/\S+)?#\d+", re.I)
TEMPLATE = ".github/pull_request_template.md"

COMMENT = re.compile(r"<!--.*?(-->|$)", re.S)  # an unclosed comment hides the rest, as on GitHub
HEADING = re.compile(r"^##\s+(.+?)\s*#*\s*$", re.M)


def sections(body):
    """{heading: text} for each `## heading` of the body, comments removed."""
    text = COMMENT.sub("", body.replace("\r\n", "\n"))
    parts = HEADING.split(text)
    return {parts[i].strip(): parts[i + 1] for i in range(1, len(parts) - 1, 2)}


def check(body):
    found = sections(body)
    if ISSUE_LINK.search(COMMENT.sub("", body)):
        found["Issue"] = "linked"
    findings = []
    for name, level in REQUIRED.items():
        text = found.get(name)
        if text is None:
            findings.append(Finding(TEMPLATE, 0, CHECK, f"section '## {name}' is missing from the PR description", level=level))
        elif not text.strip():
            findings.append(Finding(TEMPLATE, 0, CHECK, f"section '## {name}' is empty: see the hints in the template", level=level))
    return findings


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        sys.exit(report(check(f.read()), CHECK))
