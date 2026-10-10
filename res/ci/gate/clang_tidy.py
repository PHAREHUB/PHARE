#!/usr/bin/env python3
"""Turn clang-tidy output (from clang-tidy-diff, so changed lines only) into PR annotations.

usage: clang_tidy.py LOG ROOT     (LOG: clang-tidy's output, ROOT: the checkout it ran in)
exit:  0 always: warning only for now

Only the checks of .clang-tidy are reported. Compiler errors (clang-diagnostic-*) are listed in
the log, not annotated: a header has no compile command of its own, so clang-tidy guesses one
from a nearby file, and the guess can miss a define (e.g. PHARE_HAS_HIGHFIVE).
"""

import os
import re
import sys

from gitdiff import Finding, report

CHECK = "clang-tidy"

# path:line:col: warning: message [check-name]
LINE = re.compile(r"^(?P<path>.+?):(?P<line>\d+):\d+: (?P<kind>warning|error): (?P<msg>.*) \[(?P<check>[\w.,-]+)\]$")


def parse(log, root):
    root = os.path.join(os.path.abspath(root), "")
    findings, seen, not_analysed = [], set(), set()
    for text in log.splitlines():
        m = LINE.match(text)
        if not m:
            continue
        path = os.path.abspath(m["path"])
        if not path.startswith(root):
            continue  # system or dependency header
        path = path[len(root):]
        if m["check"].startswith("clang-diagnostic"):
            not_analysed.add(path)
            continue
        key = (path, int(m["line"]), m["check"])
        if key in seen:
            continue  # a header seen from several files
        seen.add(key)
        findings.append(Finding(path, int(m["line"]), CHECK, f"{m['msg']} [{m['check']}]", level="warning"))
    return findings, sorted(not_analysed)


def main(argv):
    with open(argv[1], errors="replace") as f:
        findings, not_analysed = parse(f.read(), argv[2])
    for path in not_analysed:
        print(f"compiler error, findings may be incomplete (guessed compile command?): {path}")
    report(findings, CHECK)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
