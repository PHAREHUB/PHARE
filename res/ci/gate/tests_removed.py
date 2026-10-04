#!/usr/bin/env python3
"""Flag tests a PR removes, read statically at the merge base and at the PR head.

- test names: googletest TEST/TEST_F/TEST_P/TYPED_TEST(_P), Python `test*` methods (with their
  class) and module-level `test*` functions. A test moved to another file is not a removal.
- CMake registrations: add_*test, phare_*exec and test add_subdirectory calls, per file. A changed
  call (e.g. a level moved from 9 to 11) counts as removed: it changes what runs, and when.
- new test files that no CMake file mentions (warning: never run by ctest).

usage: tests_removed.py BASE HEAD     (compares HEAD with its merge base with BASE)
exit:  0 nothing to fix, 1 an error is left after label overrides (label tests-removed-ok)
"""

import ast
import os
import re
import subprocess
import sys

from gitdiff import Finding, changed_files, file_at, git, is_test_file, merge_base, report

CHECK = "tests-removed"
LABEL = "tests-removed-ok"

GTEST = re.compile(r"\b(?:TEST|TEST_F|TEST_P|TYPED_TEST|TYPED_TEST_P)\s*\(\s*(\w+)\s*,\s*(\w+)\s*\)")
CPP_COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.S)
CPP_SOURCE = (".cpp", ".hpp", ".cc", ".h", ".ipp")

CMAKE = re.compile(r"(^|/)CMakeLists\.txt$|\.cmake$")
CMAKE_COMMENT = re.compile(r"#[^\n]*")
REGISTRATION = re.compile(
    r"\b(add_\w*test|phare_\w*exec|add_subdirectory)\s*\(", re.I
)  # add_test, add_phare_test, add_no_mpi_phare_test, add_python3_test, add_mpi_python3_test, ...


def gtest_names(source):
    """Suite.name; a DISABLED_ prefix is dropped: disabling is reported by added_lines.py."""

    def strip(s):
        return s.removeprefix("DISABLED_")

    return {f"{strip(suite)}.{strip(name)}" for suite, name in GTEST.findall(CPP_COMMENT.sub("", source))}


def python_names(source):
    """`Class.test_x` for test methods, `test_x` for module functions; None if unparsable."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return None
    names = set()

    def visit(body, prefix):
        for node in body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test"):
                names.add(prefix + node.name)
            elif isinstance(node, ast.ClassDef):
                visit(node.body, f"{prefix}{node.name}.")

    visit(tree.body, "")
    return names


def registrations(source):
    """Normalised CMake calls that register tests or test directories, in file order."""
    text = CMAKE_COMMENT.sub("", source)
    calls = []
    for m in REGISTRATION.finditer(text):
        depth, end = 0, m.end() - 1
        for end in range(m.end() - 1, len(text)):
            depth += {"(": 1, ")": -1}.get(text[end], 0)
            if depth == 0:
                break
        call = " ".join(text[m.start() : end + 1].split())
        if m.group(1).lower() != "add_subdirectory" or "test" in call.lower():
            calls.append(call)
    return calls


def test_names(path, source):
    if source is None:
        return set()
    if path.endswith(".py") and is_test_file(path):
        return python_names(source) or set()
    if path.endswith(CPP_SOURCE) and is_test_file(path):
        return gtest_names(source)
    return set()


def mentioned_in_cmake(path, head, cwd):
    """True if a CMakeLists.txt or .cmake file at head names this file (gtest source or python script)."""
    try:
        out = git("grep", "-l", "-F", os.path.basename(path), head, "--", ":(top,glob)**/CMakeLists.txt", ":(top,glob)**/*.cmake", cwd=cwd)
    except subprocess.CalledProcessError:  # git grep exits 1 when nothing matches
        return False
    return bool(out.strip())


def check(base, head, cwd=None):
    findings = []
    base_rev = merge_base(base, head, cwd=cwd)
    before, after = {}, set()  # test name -> path at base; test names at head
    added_in = {}  # path at base -> number of test names its head version adds
    for change in changed_files(base, head, cwd=cwd):
        old = file_at(base_rev, change.old_path, cwd=cwd) if change.status != "A" else None
        new = file_at(head, change.path, cwd=cwd) if change.status != "D" else None

        for name in test_names(change.old_path or "", old):
            before.setdefault(name, change.old_path)
        base_names = test_names(change.old_path or "", old)
        head_names = test_names(change.path, new)
        after |= head_names
        if change.old_path:
            added_in[change.old_path] = len(head_names - base_names)
        if change.path.endswith(".py") and is_test_file(change.path) and new is not None and python_names(new) is None:
            findings.append(Finding(change.path, 0, CHECK, "could not parse: its tests cannot be listed", label=LABEL))

        if CMAKE.search(change.path) or (change.old_path and CMAKE.search(change.old_path)):
            kept = registrations(new) if new is not None else []
            for call in registrations(old) if old is not None else []:
                if call in kept:
                    kept.remove(call)
                else:
                    findings.append(Finding(change.path, 0, CHECK, f"test registration removed or changed: {call}", label=LABEL))

        if change.status == "A" and head_names and change.path.endswith((".py", ".cpp")) and not mentioned_in_cmake(change.path, head, cwd):
            findings.append(Finding(change.path, 0, CHECK, "new test file not mentioned in any CMake file: ctest will not run it", level="warning"))

    for name in sorted(set(before) - after):
        path = before[name]
        hint = f" (this file also gains {added_in[path]} new test name(s): renamed?)" if added_in.get(path) else ""
        findings.append(Finding(path, 0, CHECK, f"test {name} removed or renamed{hint}", label=LABEL))
    return findings


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(report(check(sys.argv[1], sys.argv[2]), CHECK))
