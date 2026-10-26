#!/usr/bin/env python3
"""Run one OpenCode agent on one red build, in a scratch folder it can read but not change.

The folder holds `context.md` (build, excerpt, PR diff), `src/` (the code of the built commit) and
`opencode.json`, which allows only reading tools: no edits, no shell, no web. The answer is taken
between the REPORT markers; anything else the agent prints is dropped.
"""

import json
import os
import re
import subprocess
from pathlib import Path

# the one place the model is chosen; in CI, the repo variable CI_EXPLAINER_MODEL overrides it
MODEL = os.environ.get("CI_EXPLAINER_MODEL") or "opencode/deepseek-v4.1-flash-free"
PROMPT = Path(__file__).with_name("prompt.md")

# how each configuration builds (TeamCity steps as of 2026-10-26; guest access can't read them)
CONFIGS = {
    "Phare_Phare_BuildGithubPr": "gcc, Debug with -O3, devMode (warnings are errors), installed SAMRAI;"
    " serial tests, then MPI tests rebuilt with -march=native",
    "Phare_Phare_BuildGithubPrClang": "clang with AddressSanitizer (memory-error detector), Debug with -O3 -g3,"
    " devMode, cppcheck; serial tests, then MPI tests",
    "Phare_Phare_BuildGithubPrGccSs": "gcc, SAMRAI built from the submodule subprojects/samrai, devMode OFF,"
    " -O3 -march=native (CPU-specific instructions such as fused multiply-add: float results can differ"
    " in the last bits from the other builds); serial tests only",
    "Phare_Phare_BuildGithubPrClangSs": "clang without sanitizer, SAMRAI built from the submodule"
    " subprojects/samrai, devMode, -O3; serial tests only",
    "Phare_Next_PrChecksV1_virtual_gcc_V__1": "trial pipeline: gcc release build + fast tests, 60-min limit",
    "Phare_Next_PrChecksV1_virtual_clang_asan_V__1": "trial pipeline: clang with AddressSanitizer (memory-error detector), 60-min limit",
}

SANDBOX = {
    "$schema": "https://opencode.ai/config.json",
    "autoupdate": False,
    "share": "disabled",
    # off explicitly: language servers and formatters can run code from the PR's build files
    "lsp": False,
    "formatter": False,
    "permission": {
        **dict.fromkeys(("read", "glob", "grep", "list"), "allow"),
        **dict.fromkeys(
            (
                "edit",
                "bash",
                "webfetch",
                "websearch",
                "task",
                "external_directory",
                "todowrite",
                "question",
                "lsp",
                "skill",
                "doom_loop",  # "ask" by default: a run with nobody to answer would hang
            ),
            "deny",
        ),
    },
}

REPORT = re.compile(r"<<<REPORT\s*\n(.*?)\n\s*REPORT>>>", re.S)
FIELDS = ("Owner:", "What:", "Why:", "Fix prompt:")


def siblings_text(build, siblings):
    """One line per other build of the same commit: which configurations passed is a clue."""
    lines = [
        f"- {b['config']}: {b['status'] if b.get('state') == 'finished' else b.get('state')}"
        f" ({CONFIGS[b['config']]})"
        for b in sorted(siblings, key=lambda b: b["config"])
        if b["id"] != build["id"] and b["config"] in CONFIGS
    ]
    return "\n".join(lines) or "- (none)"


def context(build, verdict, excerpt, diff, siblings=(), max_diff=20000):
    if len(diff) > max_diff:
        diff = (
            diff[:max_diff]
            + f"\n[... diff cut, {len(diff) - max_diff} more characters; read src/ for the rest ...]"
        )
    return f"""# Red build {build["id"]}

- PR: #{build["pr"]}, commit {build["sha"]}
- Configuration: {build["config"]}: {CONFIGS.get(build["config"], "?")}
- TeamCity status: {build["status_text"]}
- Sorted as: {verdict.kind} ({verdict.reason}){", failed tests: " + ", ".join(verdict.tests) if verdict.tests else ""}

## Other builds of the same commit

{siblings_text(build, siblings)}

## Excerpt of the log

```
{excerpt}
```

## Diff of the PR (merge base with master .. built commit)

```diff
{diff}
```
"""


def prepare(workdir, context_text, snapshot):
    """snapshot: callable(dest) that writes the built commit's tree into dest."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "context.md").write_text(context_text)
    (workdir / "opencode.json").write_text(json.dumps(SANDBOX, indent=2))
    if not (workdir / "src").exists():
        snapshot(workdir / "src")
    return workdir


def parse(output):
    """The last REPORT block, if it has the four fields in order; else None."""
    blocks = REPORT.findall(output)
    if not blocks:
        return None
    block = blocks[-1].strip()
    positions = [block.find(field) for field in FIELDS]
    if -1 in positions or positions != sorted(positions):
        return None
    return block


def text_of(events):
    """The agent's text from OpenCode's JSON events (one per line); tool calls and logs dropped."""
    texts = []
    for line in events.splitlines():
        try:
            part = json.loads(line).get("part") or {}
        except (ValueError, AttributeError):
            continue
        if part.get("type") == "text":
            texts.append(part.get("text", ""))
    return "\n".join(texts)


def run(workdir, model=MODEL, timeout=600):
    """(report or None, raw JSON events); events are kept even on timeout, to see where it got stuck."""
    # absolute: it is both --dir and cwd, where a relative path would be read twice
    workdir = Path(workdir).resolve()
    env = (
        {
            k: v
            for k, v in os.environ.items()
            if not re.search(r"TOKEN|SECRET|PASSWORD|KEY", k, re.I)
        }
        | {
            "OPENCODE_DISABLE_LSP_DOWNLOAD": "1",
            "OPENCODE_DISABLE_CLAUDE_CODE": "1",  # no CLAUDE.md from the PR read as instructions
        }
    )
    message = PROMPT.read_text() + "\nStart by reading context.md."
    try:
        done = subprocess.run(
            [
                "opencode",
                "run",
                "--pure",
                "--format",
                "json",
                "--dir",
                str(workdir),
                "-m",
                model,
                message,
            ],
            cwd=workdir,
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,  # a failed run is an unusable answer, not a crash
        )
        raw = done.stdout + done.stderr
    except subprocess.TimeoutExpired as expired:
        partial = expired.stdout or b""
        raw = f"(timeout after {timeout} s)\n" + (
            partial.decode(errors="replace") if isinstance(partial, bytes) else partial
        )
    return parse(text_of(raw)), raw
