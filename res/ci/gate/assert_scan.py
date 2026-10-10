#!/usr/bin/env python3
"""Flag tests a PR makes easier to pass: removed assertions, looser tolerances, changed seeds.

Python test files are compared as syntax trees (ast), test by test.
C++ test files are compared on their EXPECT_/ASSERT_ macros, comments and #if 0 blocks left out.

usage: assert_scan.py BASE HEAD     (compares HEAD with its merge base with BASE)
exit:  0 nothing to fix, 1 an error is left after label overrides (label test-values-ok)

Tests removed or renamed are reported by tests_removed.py.
"""

import ast
import re
import sys
from collections import defaultdict

from gitdiff import Finding, changed_files, file_at, is_test_file, merge_base, report

CHECK = "test-values"
LABEL = "test-values-ok"


def flag(path, line, message):
    return Finding(path, line, CHECK, message, label=LABEL)

# tolerance-like names: a larger value passes more results
LOOSER_IF_LARGER = {"atol", "rtol", "abs_tol", "rel_tol", "delta", "tol", "tolerance", "eps", "epsilon", "threshold", "bound"}
# a smaller value passes more results (digits that must agree)
LOOSER_IF_SMALLER = {"places", "decimal", "significant"}
TOLERANCE_NAME = re.compile(r"^(atol|rtol|tol|tolerance|eps|epsilon|threshold)$|_(atol|rtol|tol|tolerance|eps)$", re.I)

# positional tolerance arguments of common assertion calls: call name -> {index: name}
POSITIONAL = {
    "assertAlmostEqual": {2: "places"},
    "assertNotAlmostEqual": {2: "places"},
    "assert_allclose": {2: "rtol", 3: "atol"},
    "allclose": {2: "rtol", 3: "atol"},
    "isclose": {2: "rtol", 3: "atol"},
    "assert_almost_equal": {2: "decimal"},
    "assert_array_almost_equal": {2: "decimal"},
    "assert_approx_equal": {2: "significant"},
}


def call_name(node):
    f = node.func
    return f.attr if isinstance(f, ast.Attribute) else f.id if isinstance(f, ast.Name) else ""


def number(node):
    try:
        value = ast.literal_eval(node)
    except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
        return None
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else None


class Scope:
    def __init__(self):
        self.asserts = 0
        self.tolerances = defaultdict(list)  # (call, name) -> [(value, line)]
        self.seeds = []  # [(value, line)]
        self.line = 0


def scan_scope(body_nodes, scope):
    for node in body_nodes:
        for sub in ast.walk(node):
            if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and sub is not node:
                continue  # nested definitions get their own scope (approximation: still walked once)
            if isinstance(sub, ast.Assert):
                scope.asserts += 1
                test = sub.test
                if isinstance(test, ast.Compare) and len(test.ops) == 1 and isinstance(test.ops[0], (ast.Lt, ast.LtE)):
                    bound = number(test.comparators[0])
                    if bound is not None:
                        scope.tolerances[("assert", "bound")].append((bound, sub.lineno))
            elif isinstance(sub, ast.Call):
                name = call_name(sub)
                if name.lower().startswith("assert") or name in ("allclose", "isclose"):
                    if name.lower().startswith("assert"):
                        scope.asserts += 1
                    for index, tol in POSITIONAL.get(name, {}).items():
                        if index < len(sub.args) and (value := number(sub.args[index])) is not None:
                            scope.tolerances[(name, tol)].append((value, sub.lineno))
                for kw in sub.keywords:
                    if kw.arg in LOOSER_IF_LARGER | LOOSER_IF_SMALLER and (value := number(kw.value)) is not None:
                        scope.tolerances[(name, kw.arg)].append((value, sub.lineno))
                    if kw.arg in ("seed", "random_state") and (value := number(kw.value)) is not None:
                        scope.seeds.append((value, sub.lineno))
                if name.endswith("seed") and sub.args and (value := number(sub.args[0])) is not None:
                    scope.seeds.append((value, sub.lineno))
            elif isinstance(sub, ast.Assign) and len(sub.targets) == 1 and isinstance(sub.targets[0], ast.Name):
                target = sub.targets[0].id
                if TOLERANCE_NAME.search(target) and (value := number(sub.value)) is not None:
                    key = "places" if "places" in target.lower() else "tol"
                    scope.tolerances[("=" + target, key)].append((value, sub.lineno))


def python_scopes(source):
    """{qualified name: Scope} for the module, each class and each function."""
    tree = ast.parse(source)
    scopes = {}

    def visit(body, prefix):
        here = Scope()
        scan_scope([n for n in body if not isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))], here)
        scopes[prefix or "<module>"] = here
        for n in body:
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
                s = Scope()
                s.line = n.lineno
                scan_scope(n.body, s)
                scopes[f"{prefix}{n.name}"] = s
            elif isinstance(n, ast.ClassDef):
                visit(n.body, f"{prefix}{n.name}.")

    visit(tree.body, "")
    return scopes


def looser(name, old, new):
    if name in LOOSER_IF_SMALLER or name == "places":
        return new < old
    return abs(new) > abs(old)


def compare_python(path, old_src, new_src):
    findings = []
    try:
        old, new = python_scopes(old_src), python_scopes(new_src)
    except SyntaxError as e:
        return [flag(path, e.lineno or 0, "could not parse (inconclusive)")]
    # a renamed class keeps its tests: fall back on the method name when it is unique in the new file
    by_name = defaultdict(list)
    for qual, scope in new.items():
        by_name[qual.split(".")[-1]].append(scope)
    for qual, before in old.items():
        short = qual.split(".")[-1]
        after = new.get(qual)
        if after is None and len(by_name[short]) == 1 and qual != "<module>":
            after = by_name[short][0]
        if after is None:
            continue  # removed or renamed tests: tests_removed.py
        line = after.line
        if after.asserts < before.asserts:
            findings.append(flag(path, line, f"{qual}: {before.asserts - after.asserts} assertion(s) removed ({before.asserts} -> {after.asserts})"))
        for key, values in before.tolerances.items():
            for (old_v, _), (new_v, new_line) in zip(values, after.tolerances.get(key, [])):
                if old_v != new_v and looser(key[1], old_v, new_v):
                    findings.append(flag(path, new_line, f"{qual}: {key[0].lstrip('=')} {key[1]} loosened {old_v:g} -> {new_v:g}"))
        if [v for v, _ in before.seeds] != [v for v, _ in after.seeds] and (before.seeds or after.seeds):
            line_no = after.seeds[0][1] if after.seeds else line
            findings.append(flag(path, line_no, f"{qual}: seed changed {[v for v, _ in before.seeds]} -> {[v for v, _ in after.seeds]}"))
    return findings


CPP_ASSERT = re.compile(r"\b(EXPECT|ASSERT)_[A-Z_]+\s*\(")
CPP_CHECK = re.compile(r"\b(EXPECT|ASSERT)_(NEAR|DOUBLE_EQ|FLOAT_EQ|EQ)\s*\(")
CPP_FLOAT = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?[fFlL]?$")
# string and char literals are matched first so that a "//" inside them is not a comment
CPP_COMMENT_OR_LITERAL = re.compile(r'"(?:\\.|[^"\\\n])*"|\'(?:\\.|[^\'\\\n])*\'|//[^\n]*|/\*.*?\*/', re.S)
PP_IF0 = re.compile(r"^\s*#\s*if\s+0\b")
PP_IF = re.compile(r"^\s*#\s*if")
PP_ELSE = re.compile(r"^\s*#\s*(else|elif)\b")
PP_ENDIF = re.compile(r"^\s*#\s*endif\b")


def cpp_code(source):
    """Lines of source with comments and #if 0 blocks blanked: a commented-out check is no check."""

    def blank(m):
        text = m.group(0)
        return text if text[0] in "\"'" else "\n" * text.count("\n")

    lines, depth = CPP_COMMENT_OR_LITERAL.sub(blank, source).split("\n"), 0
    for i, line in enumerate(lines):
        if depth == 0:
            if PP_IF0.match(line):
                depth, lines[i] = 1, ""
            continue
        if PP_ENDIF.match(line) or (depth == 1 and PP_ELSE.match(line)):
            depth -= 1  # the #else branch of an #if 0 is live code
        elif PP_IF.match(line):
            depth += 1
        lines[i] = ""
    return lines


def macro_args(text, start):
    """Top-level comma-separated arguments of the macro call opening at text[start-1] == '('."""
    depth, args, current = 1, [], ""
    for ch in text[start:]:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
            if depth == 0:
                args.append(current.strip())
                return args
        if ch == "," and depth == 1:
            args.append(current.strip())
            current = ""
        else:
            current += ch
    return None  # call continues on another line


def cpp_checks(lines):
    """{(a, b): (kind, tolerance)} for EXPECT/ASSERT _NEAR/_EQ/_DOUBLE_EQ on single lines."""
    found = {}
    for line_no, text in lines:
        for m in CPP_CHECK.finditer(text):
            args = macro_args(text, m.end())
            if not args or len(args) < 2:
                continue
            tol = None
            if m.group(2) == "NEAR" and len(args) == 3 and CPP_FLOAT.match(args[2]):
                tol = float(args[2].rstrip("fFlL"))
            found[(args[0], args[1])] = (m.group(2), tol, line_no)
    return found


def compare_cpp(path, old_src, new_src, change):
    findings = []
    old_code, new_code = cpp_code(old_src), cpp_code(new_src)
    before, after = (sum(len(CPP_ASSERT.findall(line)) for line in code) for code in (old_code, new_code))
    if after < before:
        findings.append(flag(path, 0, f"{before - after} EXPECT_/ASSERT_ check(s) removed ({before} -> {after})"))
    old_checks = cpp_checks([(n, old_code[n - 1]) for n, _ in change.removed])
    new_checks = cpp_checks([(n, new_code[n - 1]) for n, _ in change.added])
    for key, (kind, tol, line_no) in new_checks.items():
        if key not in old_checks:
            continue
        old_kind, old_tol, _ = old_checks[key]
        if old_kind in ("EQ", "DOUBLE_EQ", "FLOAT_EQ") and kind == "NEAR":
            findings.append(flag(path, line_no, f"exact check on {key[0]} replaced by a tolerance"))
        elif kind == "NEAR" and old_kind == "NEAR" and tol is not None and old_tol is not None and tol > old_tol:
            findings.append(flag(path, line_no, f"tolerance on {key[0]} loosened {old_tol:g} -> {tol:g}"))
    old_seed = [t for _, t in change.removed if re.search(r"seed|mt19937", t, re.I)]
    new_seed = [t for _, t in change.added if re.search(r"seed|mt19937", t, re.I)]
    if old_seed and new_seed and [re.findall(r"\d+", t) for t in old_seed] != [re.findall(r"\d+", t) for t in new_seed]:
        findings.append(flag(path, change.added[0][0], "random seed changed"))
    return findings


def check(base, head, cwd=None):
    findings = []
    base_rev = merge_base(base, head, cwd=cwd)
    for change in changed_files(base, head, cwd=cwd):
        if change.status not in ("M", "R") or not is_test_file(change.path):
            continue
        old_src, new_src = file_at(base_rev, change.old_path, cwd=cwd), file_at(head, change.path, cwd=cwd)
        if old_src is None or new_src is None:
            continue
        if change.path.endswith(".py"):
            findings += compare_python(change.path, old_src, new_src)
        elif change.path.endswith((".cpp", ".hpp", ".cc", ".h")):
            findings += compare_cpp(change.path, old_src, new_src, change)
    return findings


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(report(check(sys.argv[1], sys.argv[2]), CHECK))
