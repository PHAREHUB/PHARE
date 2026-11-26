"""Nightly physics report, step 'oc_judge': a free OpenCode agent reads summary.json and figures.

The agent gets read-only tools on a scratch directory holding public material only. Its answer
is taken between <<<REPORT and REPORT>>>; a deck left out or without a valid Verdict line counts
as 'cannot tell'. The agent's text goes to agent.md and the facts to meta.json; publish builds
the page from those and summary.json. This script never fails the job.
"""

import argparse
import contextlib
import json
import os
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

import decks
import tc

VERDICTS = ("looks right", "suspicious", "cannot tell")
# free models are "free for a limited time": say it plainly when the configured one is gone.
# TODO(Task 7): check against the exact message opencode 1.18.32 prints for a withdrawn model.
MODEL_GONE = re.compile(
    r"ProviderModelNotFound|unknown model|model\b[^\n]{0,80}?\bnot (?:found|available)",
    re.I,
)
# what opencode 1.18.32 logs, at level ERROR, once the free quota is spent (seen 2026-10-27);
# it then retries without end and prints nothing else
FREE_LIMIT = re.compile(r"usage limit reached", re.I)
VERSION = re.compile(r"\d+(?:\.\d+)+(?:-[0-9A-Za-z.]+)?")
# what the agent inherits, nothing else (a deny-list ages badly): locale, the lab proxy and the
# certificate paths it needs (paths, not secrets). HOME, XDG_* and OPENCODE_* are set by
# run_agent itself
PASSED_ENV = re.compile(
    r"PATH|LANG|TZ|LC_[A-Za-z_]+|(?i:(https?|no|all)_proxy)"
    r"|SSL_CERT_FILE|SSL_CERT_DIR|NODE_EXTRA_CA_CERTS|REQUESTS_CA_BUNDLE|CURL_CA_BUNDLE"
)


def opencode_config(model):
    deny = "deny"
    return {
        "$schema": "https://opencode.ai/config.json",
        "model": model,
        "autoupdate": False,
        "share": "disabled",
        # off explicitly: language servers and formatters can run code
        "lsp": False,
        "formatter": False,
        "permission": {
            "read": "allow",
            "glob": "allow",
            "grep": "allow",
            "list": "allow",
            "edit": deny,
            "bash": deny,
            "task": deny,
            "external_directory": deny,
            "webfetch": deny,
            "websearch": deny,
            "skill": deny,
            "lsp": deny,
            "question": deny,
            "todowrite": deny,
            "doom_loop": deny,  # "ask" by default: a run with nobody to answer would hang
        },
    }


def fetch_history(base, build_type_id, job, current_id, nights, get=None):
    """Previous nights, newest first (TeamCity's order): each one's summary.json decks and the
    agent_status from its meta.json (None when that file is missing)."""
    try:
        ids = [
            i
            for i in tc.recent_builds(base, build_type_id, nights + 1, get)
            if i != current_id
        ]
    except Exception:  # noqa: BLE001  history is optional: never fail the report
        return []

    def artifact(i, name):
        url = tc.artifact_url(base, build_type_id, i, tc.published(name))
        return tc.get_json(url, get)

    history = []
    for i in ids[:nights]:
        try:
            s = artifact(i, "summary.json")
        except Exception:  # noqa: BLE001  a missing night is skipped
            continue
        try:
            agent_status = artifact(i, "meta.json").get("agent_status")
        except Exception:  # noqa: BLE001  older nights have no agent_status
            agent_status = None
        history.append(
            {"build_id": i, "decks": s.get("decks", []), "agent_status": agent_status}
        )
    return history


def unavailable_streak(history):
    """Nights in a row the agent was unavailable, tonight included."""
    n = 1
    for night in history:
        if night.get("agent_status") != "unavailable":
            break
        n += 1
    return n


def last_block(stdout):
    """The agent's text between the last <<<REPORT and REPORT>>>, or None."""
    blocks = re.findall(r"<<<REPORT\s*\n(.*?)\n\s*REPORT>>>", stdout or "", re.S)
    return (
        blocks[-1] if blocks else None
    )  # the last wins: the agent may echo the template


def regular_file(path):
    """A regular file, not a symlink: copyfile follows links, so a link would hand the agent its
    target, from anywhere on the disk."""
    try:
        return stat.S_ISREG(os.lstat(path).st_mode)
    except OSError:
        return False


def prepare_workdir(root, out, summary, history, deck_list, prompt_path):
    """(work dir, names not copied because they are symlinks or not regular files)."""
    work = Path(tempfile.mkdtemp(prefix="oc_judge_"))
    (work / "summary.json").write_text(json.dumps(summary, indent=2))
    (work / "history.json").write_text(json.dumps(history, indent=2))
    (work / "sources").mkdir()
    skipped = []
    for d in deck_list:
        src = Path(root) / d["source"]
        if not os.path.lexists(src):
            continue
        if not regular_file(src):
            skipped.append(src.name)
            continue
        shutil.copyfile(src, work / "sources" / src.name)
    figures = Path(out) / "figures"
    prefixes = tuple(f"{d['name']}__" for d in deck_list)
    if figures.is_dir():
        (work / "figures").mkdir()
        for p in sorted(figures.iterdir()):
            if not p.name.startswith(prefixes):
                continue  # a deck the agent does not judge
            if not regular_file(p):
                skipped.append(p.name)
                continue
            shutil.copyfile(p, work / "figures" / p.name)
    shutil.copyfile(prompt_path, work / "PROMPT.md")
    return work, skipped


KILL_GRACE_S = 5


def run_bounded(cmd, cwd, env, timeout, stop=None):
    """(exit code, stdout, stderr), or (None, stdout, stderr so far) after a timeout or as soon
    as a stderr line matches `stop`. The agent runs in its own process group, and the whole
    group is killed: a worker it started would otherwise hold stdout open and outlive it. Every
    wait is bounded."""
    p = subprocess.Popen(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        errors="replace",
        start_new_session=True,
    )
    out, err = [], []
    stopped = threading.Event()

    def read(f, sink, watch):
        with contextlib.suppress(OSError, ValueError):
            for line in f:
                sink.append(line)
                if watch is not None and watch.search(line):
                    stopped.set()

    readers = [
        threading.Thread(target=read, args=(p.stdout, out, None), daemon=True),
        threading.Thread(target=read, args=(p.stderr, err, stop), daemon=True),
    ]

    def closed(rc):
        # a process that left the group may still hold the pipes: stop reading, never wait on it
        for f in (p.stdout, p.stderr):
            with contextlib.suppress(OSError):
                f.close()
        return rc, "".join(out), "".join(err)

    for t in readers:
        t.start()
    deadline = time.monotonic() + timeout
    while not stopped.is_set() and time.monotonic() < deadline:
        if p.poll() is not None and not any(t.is_alive() for t in readers):
            return closed(p.returncode)
        stopped.wait(0.1)
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(p.pid, sig)
        except ProcessLookupError:
            break
        with contextlib.suppress(subprocess.TimeoutExpired):
            p.wait(timeout=KILL_GRACE_S)
            break
    for t in readers:
        t.join(KILL_GRACE_S)
    return closed(None)


def run_agent(exe, work, model, timeout):
    """(stdout, None) or (None, why); never raises: the report must not make the job red."""
    # the agent's HOME is outside `work`: its own data and logs must not land where it reads
    home = Path(tempfile.mkdtemp(prefix="oc_home_"))
    try:
        (work / "opencode.json").write_text(
            json.dumps(opencode_config(model), indent=2)
        )
        env = {k: v for k, v in os.environ.items() if PASSED_ENV.fullmatch(k)}
        env.update(
            HOME=str(home),
            XDG_CONFIG_HOME=str(home / ".config"),
            XDG_DATA_HOME=str(home / ".local" / "share"),
            XDG_CACHE_HOME=str(home / ".cache"),
            XDG_STATE_HOME=str(home / ".local" / "state"),
            OPENCODE_CONFIG=str(work / "opencode.json"),
            OPENCODE_DISABLE_LSP_DOWNLOAD="1",
            OPENCODE_DISABLE_CLAUDE_CODE="1",
        )
        pngs = (
            sorted(str(p) for p in (work / "figures").glob("*.png"))
            if (work / "figures").exists()
            else []
        )
        # the message comes first and -f last: opencode takes every argument after -f as a file
        cmd = [
            exe,
            "run",
            "--pure",
            "--print-logs",
            "--log-level",
            "ERROR",
            (work / "PROMPT.md").read_text(),
            "-m",
            model,
            "--dir",
            str(work),
        ]
        if pngs:
            cmd += ["-f", *pngs]
        try:
            rc, stdout, stderr = run_bounded(cmd, work, env, timeout, stop=FREE_LIMIT)
        except (OSError, ValueError) as e:
            return None, f"could not run the agent: {type(e).__name__}"
        if rc is None and FREE_LIMIT.search(stderr):
            return None, f"the free usage limit of {model} is reached"
        if rc is None:
            return None, f"timeout after {timeout} s"
        if "<<<REPORT" not in stdout and MODEL_GONE.search(f"{stdout}\n{stderr}"):
            return None, f"model {model} is no longer available on OpenCode"
        if rc != 0:
            return None, f"opencode exit code {rc}"
        return stdout, None
    except (OSError, ValueError) as e:
        return None, f"could not prepare the agent: {type(e).__name__}"
    finally:
        shutil.rmtree(home, ignore_errors=True)


def parse_version(printed):
    """The version number in `opencode --version`'s output ("opencode 1.18.32" -> "1.18.32"):
    the header must hold one word, or publish fences it as untrusted."""
    m = VERSION.search(printed or "")
    return m.group(0) if m else "unknown"


def opencode_version(exe):
    try:
        return parse_version(
            subprocess.run(
                [exe, "--version"],
                capture_output=True,
                text=True,
                errors="replace",
                timeout=60,
                check=False,
            ).stdout
        )
    except Exception:  # noqa: BLE001  the version is only a label
        return "unknown"


def main(argv=None, get=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--job", required=True)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--opencode", default="opencode")
    ap.add_argument("--timeout", type=int, default=None)
    ap.add_argument("--config", default=None)
    ap.add_argument("--no-history", action="store_true")
    a = ap.parse_args(argv)
    cfg = decks.load(a.config)
    deck_list = decks.job(cfg, a.job)
    names = [d["name"] for d in deck_list]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    props = tc.build_properties()
    try:
        build_id = int(props.get("teamcity.build.id", 0)) or None
    except ValueError:
        build_id = None
    build_type_id = props.get("teamcity.buildType.id")
    model = cfg["model"]
    # whatever goes wrong below, meta.json is written and the exit code is 0
    version = "unknown"
    progress = {}  # judged_body notes when the agent has run, for the reason below
    try:
        version = opencode_version(a.opencode)
        info = judged_body(
            a, cfg, deck_list, names, out, build_id, build_type_id, get, progress
        )
    except Exception as e:  # noqa: BLE001  the judge is advice and never makes a job red
        when = "after" if progress.get("agent_ran") else "before"
        reason = f"{type(e).__name__} {when} the agent ran"
        info = {"agent_status": "unavailable", "agent_reason": reason}
    # publish posts the agent's text alone, fenced, and writes every fact itself from
    # meta.json and summary.json: our lines and the model's never share a block
    agent = out / "agent.md"
    agent.unlink(missing_ok=True)
    if info.get("agent_text") is not None:
        agent.write_text(info["agent_text"])
    (out / "meta.json").write_text(
        json.dumps(
            {
                "job": a.job,
                "build_id": build_id,
                "build_type_id": build_type_id,
                # the commit this job tested: publish reports on it, not on its own checkout
                "sha": props.get("build.vcs.number"),
                "model": model,
                "opencode_version": version,
                # read back by the next nights, to count the nights in a row without an agent
                "agent_status": info["agent_status"],
                "agent_reason": info["agent_reason"],
                "unavailable_nights": info.get("unavailable_nights", 0),
                "not_given": info.get("not_given", []),
            }
        )
    )
    return 0


def only_decks(summary, names):
    return summary | {
        "decks": [d for d in summary.get("decks", []) if d.get("name") in names]
    }


def ask_agent(a, cfg, out, judged, summary, history):
    """(stdout or None, why, names not copied) for the decks the agent judges."""
    names = {d["name"] for d in judged}
    prompt = Path(__file__).with_name("prompt.md")
    work, skipped = prepare_workdir(
        a.root,
        out,
        only_decks(summary, names),
        [only_decks(h, names) for h in history],
        judged,
        prompt,
    )
    try:
        stdout, why = run_agent(
            a.opencode, work, cfg["model"], a.timeout or cfg.get("judge_timeout_s", 900)
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return stdout, why, skipped


def judged_body(
    a, cfg, deck_list, names, out, build_id, build_type_id, get=None, progress=None
):
    """info: agent_status ("ok" | "unavailable" | "not run"), agent_reason, unavailable_nights,
    not_given and agent_text (the agent's block, or None)."""
    summary = json.loads((out / "summary.json").read_text())
    if build_type_id and not a.no_history:
        # visible in the build log: a new pipeline version would change the id, and reset history
        print(f"history: build type {build_type_id}")
    history = (
        []
        if a.no_history or not build_type_id
        else fetch_history(
            cfg["teamcity"]["url"],
            build_type_id,
            a.job,
            build_id,
            cfg.get("history_nights", 7),
            get,
        )
    )
    # a deck with judge: false is not shown to the agent at all
    judged = [d for d in deck_list if d.get("judge", True)]
    skipped = []
    info = {"agent_status": "not run", "agent_reason": "no deck to judge"}
    if judged:
        stdout, why, skipped = ask_agent(a, cfg, out, judged, summary, history)
        if progress is not None:
            progress["agent_ran"] = True
        text = last_block(stdout) if stdout else None
        if text is not None:
            info = {"agent_status": "ok", "agent_reason": None, "agent_text": text}
        else:
            info = {
                "agent_status": "unavailable",
                "agent_reason": why or "no <<<REPORT block in the agent's answer",
                "unavailable_nights": unavailable_streak(history),
            }
    info["not_given"] = skipped
    return info


if __name__ == "__main__":
    sys.exit(main())
