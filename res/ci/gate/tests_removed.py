#!/usr/bin/env python3
"""Flag tests a PR removes, comparing the merge base with the PR head.

- test names: googletest TEST/TEST_F/TEST_P/TYPED_TEST(_P), Python `test*` methods (with their
  class) and module-level `test*` functions, read statically. A test moved to another file is not a removal.
- CMake registrations: add_*test, phare_*exec and test add_subdirectory calls, per file. A changed
  call (e.g. a level moved from 9 to 11) counts as removed: it changes what runs, and when.
- new test files that no CMake file mentions, by name or built from variables in a loop or a
  function (warning: never run by ctest).

With --listed (list_tests.py run on both trees, as the gate workflow does), what ctest and unittest
would actually run replaces the static reading of CMake calls and of Python scripts:
- a ctest test no longer registered, or run on fewer MPI processes: error. Same name, other
  command: warning (e.g. a wrapper script). A level change is not seen: every level is configured.
- a Python Class.method with fewer cases (ddt @data values count separately): error.
- a Python script that no longer imports: error, and its tests are read statically.

usage: tests_removed.py BASE HEAD [--listed BASE.json HEAD.json]   (HEAD vs its merge base with BASE)
exit:  0 nothing to fix, 1 an error is left after label overrides (label tests-removed-ok)
"""

import ast
import json
import os
import re
import subprocess
import sys
from collections import Counter

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
CMAKE_TOKEN = re.compile(r'[^\s()"]+')
CMAKE_VAR = re.compile(r"\$\{[^}]*\}")


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


def registrations(source, path):
    """Normalised CMake calls that register tests or test directories, in file order.

    Under tests/ every add_subdirectory counts (tests/amr/data/particles adds `copy`, `refine`, ...);
    elsewhere only those whose argument names a test directory.
    """
    text = CMAKE_COMMENT.sub("", source)
    calls = []
    for m in REGISTRATION.finditer(text):
        depth, end = 0, m.end() - 1
        for end in range(m.end() - 1, len(text)):
            depth += {"(": 1, ")": -1}.get(text[end], 0)
            if depth == 0:
                break
        call = " ".join(text[m.start() : end + 1].split())
        call = re.sub(r"(?<![\w.])(\d+),(?=\s)", r"\1", call)  # "11," -> "11": CMake reads the leading number
        if m.group(1).lower() != "add_subdirectory" or "test" in call.lower() or is_test_file(path):
            calls.append(call)
    return calls


def test_names(path, source, listed=()):
    """Statically read test names; none for Python scripts in `listed` (counted by list_tests.py)."""
    if source is None or path in listed:
        return set()
    if path.endswith(".py") and is_test_file(path):
        return python_names(source) or set()
    if path.endswith(CPP_SOURCE) and is_test_file(path):
        return gtest_names(source)
    return set()


def builds_name(cmake, name):
    """True if CMake text builds this file name from variables, as a loop or a function does.

    - its stem as a whole word: `_add_test(test_x)` with `add_executable(${n} ${n}.cpp)` in the function
    - a token with variables and some fixed text besides the extension: `test_x_${suffix}.cpp` in a
      foreach; `${n}.cpp` alone would match any file and is skipped
    """
    text = CMAKE_COMMENT.sub("", cmake)
    stem = os.path.splitext(name)[0]
    if re.search(rf"(?<![\w.-]){re.escape(stem)}(?![\w.-])", text):
        return True
    for token in CMAKE_TOKEN.findall(text):
        literals = CMAKE_VAR.split(token.rsplit("/", 1)[-1])
        fixed = re.sub(r"\.\w+$", "", "".join(literals))  # fixed text without the extension
        if len(literals) > 1 and re.search(r"\w", fixed) and re.fullmatch(r"[\w-]*".join(map(re.escape, literals)), name):
            return True
    return False


def mentioned_in_cmake(path, head, cwd):
    """True if a CMake file at head names this file (gtest source or python script).

    The name anywhere in a CMakeLists.txt or .cmake file, or built from variables in the
    CMakeLists.txt of its directory or of a directory above it.
    """
    try:
        out = git("grep", "-l", "-F", os.path.basename(path), head, "--", ":(top,glob)**/CMakeLists.txt", ":(top,glob)**/*.cmake", cwd=cwd)
        if out.strip():
            return True
    except subprocess.CalledProcessError:  # git grep exits 1 when nothing matches
        pass
    directory = os.path.dirname(path)
    while True:
        cmake = file_at(head, f"{directory}/CMakeLists.txt" if directory else "CMakeLists.txt", cwd=cwd)
        if cmake is not None and builds_name(cmake, os.path.basename(path)):
            return True
        if not directory:
            return False
        directory = os.path.dirname(directory)


def python_cases(listing, usable):
    """Counter of Class.method cases over the usable scripts, and {Class.method: {script: n}}."""
    total, where = Counter(), {}
    for script in usable:
        for name, n in listing.get(script, {}).items():
            total[name] += n
            where.setdefault(name, {})[script] = n
    return total, where


def listed_findings(base_list, head_list, cmake_path):
    """Findings from what ctest and unittest run at base and head; the scripts they counted."""
    findings = []
    before, after = base_list["ctest"], head_list["ctest"]
    commands = []
    for name in sorted(before):
        if name not in after:
            findings.append(Finding(cmake_path, 0, CHECK, f"ctest test {name} is no longer registered", label=LABEL))
        elif after[name]["processors"] < before[name]["processors"]:
            findings.append(Finding(cmake_path, 0, CHECK, f"ctest test {name} runs on {after[name]['processors']} MPI process(es) instead of {before[name]['processors']}", label=LABEL))
        elif after[name]["command"] != before[name]["command"]:
            commands.append(name)
    if commands:
        name = commands[0]
        findings.append(Finding(cmake_path, 0, CHECK, f"{len(commands)} ctest test(s) keep their name but run another command ({', '.join(commands)}), e.g. {name}: `{' '.join(before[name]['command'])}` -> `{' '.join(after[name]['command'])}`", level="warning"))

    base_py, head_py = base_list["python"], head_list["python"]
    for script in sorted(head_py):
        if "error" in head_py[script]:
            findings.append(Finding(script, 0, CHECK, f"cannot be imported, its tests cannot be listed: {head_py[script]['error']}", label=LABEL))
    usable = {s for s in set(base_py) | set(head_py) if "error" not in base_py.get(s, {}) and "error" not in head_py.get(s, {})}
    (was, was_in), (now, now_in) = python_cases(base_py, usable), python_cases(head_py, usable)
    for name in sorted(n for n in was if now[n] < was[n]):
        path = next((s for s, n in was_in[name].items() if now_in.get(name, {}).get(s, 0) < n), next(iter(was_in[name])))
        lost = "removed or renamed" if now[name] == 0 else f"has {now[name]} case(s) instead of {was[name]}"
        findings.append(Finding(path, 0, CHECK, f"test {name} {lost}", label=LABEL))
    return findings, usable


def check(base, head, cwd=None, listed=None):
    """listed: (base listing, head listing) from list_tests.py, or None to read everything statically."""
    findings = []
    base_rev = merge_base(base, head, cwd=cwd)
    changes = changed_files(base, head, cwd=cwd)
    listed_py = set()
    if listed:
        cmake_path = next((c.path for c in changes if CMAKE.search(c.path)), "CMakeLists.txt")
        found, listed_py = listed_findings(*listed, cmake_path)
        findings += found
    # a name can be in several files: count the changed files that have it, so that a move
    # keeps the count and a deletion of one of the copies lowers it
    before, after = Counter(), Counter()
    lost_in = {}  # test name -> first path at base whose head version no longer has it
    added_in = {}  # path at base -> number of test names its head version adds
    for change in changes:
        old = file_at(base_rev, change.old_path, cwd=cwd) if change.status != "A" else None
        new = file_at(head, change.path, cwd=cwd) if change.status != "D" else None

        base_names = test_names(change.old_path or "", old, listed_py)
        head_names = test_names(change.path, new, listed_py)
        before.update(base_names)
        after.update(head_names)
        for name in base_names - head_names:
            lost_in.setdefault(name, change.old_path)
        if change.old_path:
            added_in[change.old_path] = len(head_names - base_names)
        if change.path not in listed_py and change.path.endswith(".py") and is_test_file(change.path) and new is not None and python_names(new) is None:
            findings.append(Finding(change.path, 0, CHECK, "could not parse: its tests cannot be listed", label=LABEL))

        if not listed and (CMAKE.search(change.path) or (change.old_path and CMAKE.search(change.old_path))):
            kept = registrations(new, change.path) if new is not None else []
            for call in registrations(old, change.old_path) if old is not None else []:
                if call in kept:
                    kept.remove(call)
                else:
                    findings.append(Finding(change.path, 0, CHECK, f"test registration removed or changed: {call}", label=LABEL))

        if change.status == "A" and (head_names or (change.path in listed_py and python_names(new))) and change.path.endswith((".py", ".cpp")) and not mentioned_in_cmake(change.path, head, cwd):
            findings.append(Finding(change.path, 0, CHECK, "new test file not mentioned in any CMake file: ctest will not run it", level="warning"))

    for name in sorted(n for n in before if after[n] < before[n]):
        path = lost_in[name]
        hint = f" (this file also gains {added_in[path]} new test name(s): renamed?)" if added_in.get(path) else ""
        findings.append(Finding(path, 0, CHECK, f"test {name} removed or renamed{hint}", label=LABEL))
    return findings


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) == 5 and args[2] == "--listed":
        listed = []
        for path in args[3:]:
            with open(path) as f:
                listed.append(json.load(f))
        sys.exit(report(check(args[0], args[1], listed=listed), CHECK))
    if len(args) != 2:
        sys.exit(__doc__)
    sys.exit(report(check(*args), CHECK))
