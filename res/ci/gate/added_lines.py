#!/usr/bin/env python3
"""Flag added lines and added files that break a rule of the repo, from one table of rules.

- tests switched off: new skips, disabled or commented-out tests, .off renames (label tests-skip-ok);
- layering: src/core includes nothing from SAMRAI, amr, simulator or python3 (no override);
- conventions: new .cpp under src/ (code is header-only), `using namespace` at namespace scope
  in a header (label convention-ok).

usage: added_lines.py BASE HEAD     (compares HEAD with its merge base with BASE)
exit:  0 nothing to fix, 1 an error is left after label overrides
"""

import re
import sys
from dataclasses import dataclass

from gitdiff import Finding, changed_files, is_test_file, report

CHECK = "added-lines"

CMAKE = re.compile(r"(^|/)CMakeLists\.txt$|\.cmake$")


def tests_or_cmake(path):
    return is_test_file(path) or bool(CMAKE.search(path))


def under(prefix, suffixes=None):
    return lambda path: path.startswith(prefix) and (suffixes is None or path.endswith(suffixes))


@dataclass
class Rule:
    applies: object  # path -> bool
    pattern: re.Pattern
    message: str
    label: str = None  # None: cannot be overridden


SKIP = "tests-skip-ok"
CONVENTION = "convention-ok"

LINE_RULES = [
    # tests switched off
    Rule(tests_or_cmake, re.compile(r"@(unittest\.)?skip(If|Unless)?\b|\bskipTest\(|\bSkipTest\b"), "Python test skipped", SKIP),
    Rule(tests_or_cmake, re.compile(r"pytest\.mark\.(skip|skipif|xfail)\b|\bpytest\.skip\("), "pytest test skipped", SKIP),
    Rule(tests_or_cmake, re.compile(r"\bGTEST_SKIP\b"), "googletest test skipped", SKIP),
    Rule(tests_or_cmake, re.compile(r"\b(TEST|TEST_F|TEST_P|TYPED_TEST)\s*\(\s*\w*\s*,\s*DISABLED_"), "googletest test disabled", SKIP),
    Rule(tests_or_cmake, re.compile(r"\b(TEST|TEST_F|TEST_P|TYPED_TEST)\s*\(\s*DISABLED_"), "googletest suite disabled", SKIP),
    Rule(tests_or_cmake, re.compile(r"\bDISABLED\s+(ON|TRUE|1)\b", re.I), "ctest test disabled", SKIP),
    Rule(tests_or_cmake, re.compile(r"^\s*#+\s*def\s+test\w*\s*\("), "Python test commented out", SKIP),
    Rule(tests_or_cmake, re.compile(r"^\s*//+\s*(TEST|TEST_F|TEST_P|TYPED_TEST)\s*\("), "googletest test commented out", SKIP),
    Rule(tests_or_cmake, re.compile(r"^\s*#+\s*(add_test|add_\w*test|phare_\w*exec)\s*\("), "test registration commented out", SKIP),
    # layering: core is the bottom layer
    Rule(under("src/core/"), re.compile(r'^\s*#\s*include\s*[<"](SAMRAI|amr|simulator|python3)/'), "src/core must not include SAMRAI, amr, simulator or python3 headers"),
    # conventions
    Rule(under("src/", (".hpp", ".h")), re.compile(r"^using\s+namespace\b"), "`using namespace` at namespace scope in a header leaks into every includer", CONVENTION),
]

DISABLED_SUFFIX = re.compile(r"\.(off|disabled|bak|skip)$", re.I)


def file_rules(change):
    """Findings about the file itself rather than one of its lines."""
    if change.status == "R" and DISABLED_SUFFIX.search(change.path) and tests_or_cmake(change.old_path):
        yield Finding(change.path, 0, CHECK, f"renamed from {change.old_path}: no longer picked up", label=SKIP)
    if change.status == "A" and change.path.startswith("src/") and change.path.endswith(".cpp"):
        yield Finding(change.path, 0, CHECK, "new .cpp under src/: PHARE is header-only", label=CONVENTION)


def check(base, head, cwd=None):
    findings = []
    for change in changed_files(base, head, cwd=cwd):
        findings += file_rules(change)
        rules = [r for r in LINE_RULES if r.applies(change.path)]
        for line_no, text in change.added if rules else []:
            rule = next((r for r in rules if r.pattern.search(text)), None)
            if rule:
                findings.append(Finding(change.path, line_no, CHECK, f"{rule.message}: {text.strip()}", label=rule.label))
    return findings


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(report(check(sys.argv[1], sys.argv[2]), CHECK))
