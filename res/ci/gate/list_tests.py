#!/usr/bin/env python3
"""List the tests a configured PHARE tree would run, for `tests_removed.py --listed`.

- ctest: every registered test, with its command and MPI processes (`ctest --show-only=json-v1`).
  Configure with every level and MPI on: a test filtered out at configure time would look removed.
- Python: the unittest cases of every test_*.py script, by importing it (no test runs), counted per
  Class.method: the cases ddt generates count under the method they come from, so inserting a
  @data value (which renumbers the cases after it) changes nothing, and dropping one lowers the count.

usage: list_tests.py TREE BUILD_DIR > listed.json
"""

import importlib.util
import inspect
import json
import os
import subprocess
import sys
import unittest
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

from gitdiff import git, is_test_file

IMPORT_TIMEOUT = 120  # seconds per script


def ctest_tests(tree, build):
    """{name: {"command": [...], "processors": N}}, paths made relative to the tree and build dirs."""
    out = subprocess.run(["ctest", "--show-only=json-v1"], cwd=build, check=True, capture_output=True, text=True).stdout
    build, tree = os.path.realpath(build), os.path.realpath(tree)

    def relative(arg):  # build first: it is inside the tree
        return arg.replace(build, "<build>").replace(tree, "<tree>")

    tests = {}
    for test in json.loads(out)["tests"]:
        props = {p["name"]: p["value"] for p in test.get("properties", [])}
        tests[test["name"]] = {
            "command": [relative(a) for a in test.get("command", [])],
            "processors": int(props.get("PROCESSORS", 1)),
        }
    return tests


def python_scripts(tree):
    return sorted(
        p for p in git("ls-files", "*.py", cwd=tree).splitlines() if is_test_file(p) and os.path.basename(p).startswith("test_")
    )


def python_cases(path):
    """{Class.method: number of cases} for one script, imported but not run (in a child process)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(path)))
    spec = importlib.util.spec_from_file_location("__gate_listing__", path)  # not __main__: no unittest.main()
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    def walk(suite):
        for t in suite:
            yield from walk(t) if isinstance(t, unittest.TestSuite) else [t]

    cases = Counter()
    for test in walk(unittest.TestLoader().loadTestsFromModule(module)):
        method = getattr(test, test._testMethodName, None)
        name = inspect.unwrap(method).__name__ if method else test._testMethodName  # ddt: original method
        cases[f"{type(test).__name__}.{name}"] += 1
    return dict(cases)


def list_script(tree, path):
    """{Class.method: n}, or {"error": message} when the script cannot be imported."""
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([os.path.join(tree, "pyphare"), tree, os.environ.get("PYTHONPATH", "")]))
    try:
        p = subprocess.run(
            [sys.executable, os.path.abspath(__file__), "--script", os.path.join(tree, path)],
            cwd=os.path.join(tree, os.path.dirname(path)), env=env, capture_output=True, text=True, timeout=IMPORT_TIMEOUT, check=False,
        )
    except subprocess.TimeoutExpired:
        return {"error": f"import took more than {IMPORT_TIMEOUT} s"}
    if p.returncode != 0:
        lines = p.stderr.strip().splitlines()
        return {"error": lines[-1] if lines else f"exit code {p.returncode}"}
    return json.loads(p.stdout.strip().splitlines()[-1])  # the script may print while imported


def main(tree, build):
    tree = os.path.abspath(tree)
    scripts = python_scripts(tree)
    with ThreadPoolExecutor(os.cpu_count()) as pool:
        python = dict(zip(scripts, pool.map(lambda s: list_script(tree, s), scripts)))
    json.dump({"ctest": ctest_tests(tree, build), "python": python}, sys.stdout, indent=1, sort_keys=True)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--script":
        print("\n" + json.dumps(python_cases(sys.argv[2])))
    elif len(sys.argv) == 3:
        main(*sys.argv[1:])
    else:
        sys.exit(__doc__)
