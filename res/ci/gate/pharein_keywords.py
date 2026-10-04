#!/usr/bin/env python3
"""Flag keywords that user scripts can no longer pass to pharein.

pharein validates **kwargs against hard-coded lists (`accepted_keywords` in simulation.py's
`checker`, `valid_keys`, `mandatory_keywords`, ...). A string removed from such a list breaks every
user script that passes it, and no signature-based tool (griffe) can see it. Added strings are
listed as notices, for the changelog.

usage: pharein_keywords.py BASE HEAD     (compares HEAD with its merge base with BASE)
exit:  0 nothing to fix, 1 an error is left after label overrides (label api-break-ok)
"""

import ast
import re
import sys
from collections import defaultdict

from gitdiff import Finding, changed_files, file_at, merge_base, report

CHECK = "pharein-keywords"
LABEL = "api-break-ok"
PHAREIN = "pyphare/pyphare/pharein/"
KEYWORD_LIST = re.compile(r"keywords$|^valid_\w*keys$|^valid_\w*(options|names)$")


def strings(node):
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        return {e.value for e in node.elts if isinstance(e, ast.Constant) and isinstance(e.value, str)}
    return set()


def keyword_lists(source):
    """{scope: {keyword: line}}: the keywords of every list built from literals, per function.

    Lists of one function are merged: moving a keyword from mandatory_keywords to
    accepted_keywords relaxes the API, it does not break it.
    """
    found = defaultdict(dict)

    def visit(node, scope):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                visit(child, f"{scope}{child.name}.")
                continue
            if isinstance(child, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
                targets = child.targets if isinstance(child, ast.Assign) else [child.target]
                for target in targets:
                    if isinstance(target, ast.Name) and KEYWORD_LIST.search(target.id) and child.value is not None:
                        for kw in strings(child.value):
                            found[scope.rstrip(".") or "<module>"].setdefault(kw, child.lineno)
            visit(child, scope)

    visit(ast.parse(source), "")
    return found


def check(base, head, cwd=None):
    findings = []
    base_rev = merge_base(base, head, cwd=cwd)
    for change in changed_files(base, head, cwd=cwd):
        old_path = change.old_path or ""
        if not (change.path.startswith(PHAREIN) or old_path.startswith(PHAREIN)) or not change.path.endswith(".py"):
            continue
        old = file_at(base_rev, old_path, cwd=cwd) if change.status != "A" else None
        new = file_at(head, change.path, cwd=cwd) if change.status != "D" else None
        try:
            before = keyword_lists(old) if old else {}
            after = keyword_lists(new) if new else {}
        except SyntaxError as e:
            findings.append(Finding(change.path, e.lineno or 0, CHECK, "could not parse: keywords cannot be compared", label=LABEL))
            continue
        for scope in sorted(set(before) | set(after)):
            old_kw, new_kw = before.get(scope, {}), after.get(scope, {})
            for kw in sorted(set(old_kw) - set(new_kw)):
                findings.append(Finding(change.path, 0, CHECK, f"keyword '{kw}' no longer accepted by {scope}: user scripts passing it now fail", label=LABEL))
            for kw in sorted(set(new_kw) - set(old_kw)):
                findings.append(Finding(change.path, new_kw[kw], CHECK, f"keyword '{kw}' now accepted by {scope}: document it", level="notice"))
    return findings


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(report(check(sys.argv[1], sys.argv[2]), CHECK))
