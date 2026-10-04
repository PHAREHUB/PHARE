"""Read a git diff between two revisions, as the diff checks need it.

The checks read the PR as text: they never build or run its code.
"""

import json
import os
import re
import subprocess
from dataclasses import dataclass, field


def git(*args, cwd=None):
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout


@dataclass
class FileChange:
    status: str  # A, M, D, R
    path: str  # path at head (path at base for D)
    old_path: str  # path at base (None for A)
    added: list = field(default_factory=list)  # [(line number at head, text)]
    removed: list = field(default_factory=list)  # [(line number at base, text)]


def changed_files(base, head, cwd=None):
    """Files changed between the merge base of base/head and head, with their +/- lines."""
    out = git("diff", "--name-status", "-M", f"{base}...{head}", cwd=cwd)
    changes = {}
    for line in out.splitlines():
        parts = line.split("\t")
        status = parts[0][0]
        if status == "R":
            changes[parts[2]] = FileChange("R", parts[2], parts[1])
        elif status == "D":
            changes[parts[1]] = FileChange("D", parts[1], parts[1])
        elif status == "A":
            changes[parts[1]] = FileChange("A", parts[1], None)
        else:
            changes[parts[1]] = FileChange(status, parts[1], parts[1])

    diff = git("diff", "-M", "-U0", "--no-color", f"{base}...{head}", cwd=cwd)
    current, old_no, new_no = None, 0, 0
    for line in diff.splitlines():
        if line.startswith("diff --git"):
            current = None
        elif line.startswith("--- "):
            path = line[4:]
            current = None
            if path != "/dev/null":
                current = next(
                    (c for c in changes.values() if c.old_path == path[2:]), None
                )
        elif line.startswith("+++ "):
            path = line[4:]
            if path != "/dev/null":  # deleted file: keep the match made on "---"
                current = changes.get(path[2:])
        elif line.startswith("@@"):
            m = re.match(r"@@ -(\d+)(?:,\d+)? \+(\d+)(?:,\d+)? @@", line)
            old_no, new_no = int(m.group(1)), int(m.group(2))
        elif current is None:
            continue
        elif line.startswith("+"):
            current.added.append((new_no, line[1:]))
            new_no += 1
        elif line.startswith("-"):
            current.removed.append((old_no, line[1:]))
            old_no += 1
    return list(changes.values())


def file_at(rev, path, cwd=None):
    """Content of path at rev, or None if it does not exist there."""
    try:
        return git("show", f"{rev}:{path}", cwd=cwd)
    except subprocess.CalledProcessError:
        return None


def merge_base(base, head, cwd=None):
    return git("merge-base", base, head, cwd=cwd).strip()


TEST_PATH = re.compile(r"(^|/)(tests?|pyphare_tests)/|(^|/)test_[^/]*$|_test\.[^/]*$")


def is_test_file(path):
    return bool(path) and bool(TEST_PATH.search(path))


@dataclass
class Finding:
    path: str
    line: int
    check: str
    message: str
    level: str = "error"  # error (fails the check), warning, notice
    label: str = None  # PR label that turns this error into a warning

    def __str__(self):
        return f"{self.path}:{self.line}: [{self.check}] {self.level}: {self.message}"

    def annotation(self):
        """GitHub workflow command; escaping as in actions/toolkit (command.ts)."""

        def data(s):
            return s.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")

        def prop(s):
            return data(s).replace(":", "%3A").replace(",", "%2C")

        where = f"file={prop(self.path)}" + (f",line={self.line}" if self.line else "")
        return f"::{self.level} {where},title={prop(self.check)}::{data(self.message)}"


def pr_labels():
    """Labels on the PR, from GATE_LABELS (a JSON list set by the workflow)."""
    return set(json.loads(os.environ.get("GATE_LABELS") or "[]"))


def report(findings, title):
    """Print findings; return the process exit code: 1 if an error is left after label overrides.

    In GitHub Actions the findings are also printed as annotations, shown on the PR's changed lines.
    """
    labels = pr_labels()
    for f in findings:
        if f.level == "error" and f.label in labels:
            f.level, f.message = "warning", f"{f.message} (accepted: label {f.label})"
    errors = [f for f in findings if f.level == "error"]
    print(f"{title}: {len(errors)} error(s), {len(findings) - len(errors)} other finding(s)")
    for f in findings:
        print(f"  {f}")
    if os.environ.get("GITHUB_ACTIONS") == "true":
        for f in findings:
            print(f.annotation())
    for label in sorted({f.label for f in errors if f.label}):
        print(f"{title}: a maintainer can accept these with the PR label '{label}'")
    return 1 if errors else 0
