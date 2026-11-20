"""Nightly physics report, step 'extract': read each pilot deck's output, write summary.json.

Deterministic. For each deck: its status, the numbers that matter, its figures, the end of its log.
A deck that crashed or left no output is reported, never raised: the other decks still count.
Freshness fails closed: the tests step touches <root>/.t0 before running the decks; any figure
or h5 older than it is from an earlier run and ignored. Without .t0, freshness is "unknown".
"""

import argparse
import glob
import importlib.util
import json
import os
import re
import shutil
import stat
import sys
import traceback
from pathlib import Path

import decks
import metrics

LOG_TAIL_LINES = 60
# the run directory, under the deck's workdir, of the h5 files each numeric extractor reads
H5_DIRS = {"ion_ion_beam": "ion_ion_beam1d", "tangential_discontinuity": "td_noflow"}
# figure names go raw into markdown on GitHub: nothing that could close or open a link
SAFE_NAME = re.compile(r"[A-Za-z0-9._-]+")


def load_module(path):
    spec = importlib.util.spec_from_file_location(Path(path).stem, path)
    if spec is None or not Path(path).exists():
        raise FileNotFoundError(path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def numbers_ion_ion_beam(deck, root, workdir, out):
    mod = load_module(root / deck["source"])
    # the deck's omega is left out: |fft[1]| of a circularly polarised mode has no frequency
    _, _, _, gamma, _, _ = mod.growth_b_right_hand(
        str(workdir / H5_DIRS["ion_ion_beam"]), 10.0
    )
    t = deck["theory"]
    return metrics.growth_rate_summary(gamma, t["value"], t["deck_tolerance"])


def numbers_tangential_discontinuity(deck, root, workdir, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pyphare.pharesee.hierarchy.fromh5 import get_times_from_h5
    from pyphare.pharesee.run import Run

    path = workdir / H5_DIRS["tangential_discontinuity"]
    times = get_times_from_h5(str(path / "EM_B.h5"))
    t0, t1 = times[0], times[-1]
    r = Run(str(path))
    b0, b1 = r.GetB(t0, merged=True), r.GetB(t1, merged=True)
    x = b0["By"][1][0]
    by0, by1 = b0["By"][0](x), b1["By"][0](x)
    length = float(r.GetDomainSize()[0])
    # no flow: the profile should not move
    result = {
        "td_profile_error": metrics.shifted_profile_error(x, by0, by1, 0.0, length),
        "final_time": float(t1),
    }
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, by0, "k", label=f"t = {t0:g}")
    ax.plot(x, by1, label=f"t = {t1:g}")
    ax.set_xlabel("x")
    ax.set_ylabel("By")
    ax.legend()
    fig.tight_layout()
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures / f"{deck['name']}__by.png", dpi=120)
    plt.close(fig)
    return result


def numbers_figures_only(deck, root, workdir, out):
    return {}


NUMBERS = {
    "ion_ion_beam": numbers_ion_ion_beam,
    "tangential_discontinuity": numbers_tangential_discontinuity,
    "figures_only": numbers_figures_only,
}


def t0_mtime(root):
    """When the tests step started (mtime of <root>/.t0), or None when it left no .t0."""
    try:
        return os.stat(Path(root) / ".t0").st_mtime
    except OSError:
        return None


def stale(path, t0):
    return t0 is not None and os.lstat(path).st_mtime < t0


def skip_reason(src, t0=None):
    """Why a figure is not copied, or None. A symlink is never copied: copyfile would follow it
    and put its target, from anywhere on the disk, in front of the agent."""
    if not SAFE_NAME.fullmatch(Path(src).name):
        return "unsafe file name"
    if os.path.islink(src):
        return "symlink"
    if not stat.S_ISREG(os.lstat(src).st_mode):
        return "not a regular file"
    if stale(src, t0):
        return "older than .t0"
    return None


def stale_h5(deck, workdir, t0):
    """'dir/name: older than .t0' for each h5 file the deck's extractor would read."""
    sub = H5_DIRS.get(deck["extractor"])
    if sub is None or t0 is None:
        return []
    return [
        f"{sub}/{p.name}: older than .t0"
        for p in sorted((workdir / sub).glob("*.h5"))
        if stale(p, t0)
    ]


def copy_figures(deck, workdir, out, t0=None):
    """(copied names, skipped 'name: reason' lines)."""
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    cap = deck.get("max_figures", decks.DEFAULT_MAX_FIGURES)
    names, skipped = [], []
    for pattern in deck["figures"]:
        for src in sorted(glob.glob(str(workdir / pattern))):
            why = skip_reason(src, t0) or (
                f"over max_figures {cap}" if len(names) >= cap else None
            )
            if why:
                skipped.append(f"{Path(src).name}: {why}")
                continue
            name = f"{deck['name']}__{Path(src).name}"
            shutil.copyfile(src, figures / name)
            names.append(name)
    return names, skipped


def find_logs(root, workdir, stem):
    # phlop writes a test class's logs under .phlop/<workdir>/<stem>/<Class>/, and a plain
    # script's next to it, named after its directory: .phlop/<workdir>.stdout (and .stderr)
    base = Path(root) / ".phlop" / workdir
    per_dir = [base.with_name(base.name + ext) for ext in (".stdout", ".stderr")]
    return sorted(
        p for p in base.rglob("*") if p.is_file() and stem in str(p.relative_to(base))
    ) + [p for p in per_dir if p.is_file()]


def log_tail(paths):
    lines = []
    for p in paths:
        lines += [f"== {p.name}"] + p.read_text(errors="replace").splitlines()[
            -LOG_TAIL_LINES:
        ]
    return "\n".join(lines[-LOG_TAIL_LINES:])


def error_status(e):
    """A one-line status for a deck's failure, paths cut to the file's name: h5py's message
    for a missing file holds the whole path and runs over a line."""
    if isinstance(e, FileNotFoundError):
        name = e.filename or next(iter(re.findall(r"name = '([^']+)'", str(e))), None)
        if name:
            return f"extract error: {os.path.basename(name)} is missing"
    return f"extract error: {type(e).__name__}: " + re.sub(r"/\S*/", "", str(e))


def run_deck(deck, root, out):
    root, out = Path(root), Path(out)
    workdir = root / deck["workdir"]
    t = deck["theory"]
    result = {
        "name": deck["name"],
        "quantity": t["quantity"],
        "theory": t["value"],
        "deck_tolerance": t["deck_tolerance"],
        "own_check": deck.get("own_check", "deck asserts"),
        "numbers": {},
        "figures": [],
        "skipped": [],
        "log_tail": "",
    }
    t0 = t0_mtime(root)
    result["freshness"] = "unknown" if t0 is None else "ok"
    stem = Path(deck["source"]).stem
    result["log_tail"] = log_tail(find_logs(root, deck["workdir"], stem))
    if not workdir.exists():
        result["status"] = "output missing"
        return result
    try:
        result["figures"], result["skipped"] = copy_figures(deck, workdir, out, t0)
        old_h5 = stale_h5(deck, workdir, t0)
        if old_h5:
            # an earlier run's numbers are not tonight's: not read at all
            result["skipped"] += old_h5
            result["status"] = "stale output ignored"
            return result
        result["numbers"] = NUMBERS[deck["extractor"]](deck, root, workdir, out)
        result["figures"] = sorted(
            set(result["figures"])
            | {p.name for p in (out / "figures").glob(f"{deck['name']}__*")}
        )
    except Exception as e:  # noqa: BLE001  one deck's failure must not hide the others
        result["status"] = error_status(e)
        result["traceback"] = traceback.format_exc(limit=3)
        return result
    if deck["extractor"] == "figures_only" and not result["figures"]:
        old = any(x.endswith(": older than .t0") for x in result["skipped"])
        result["status"] = "stale output ignored" if old else "output missing"
    else:
        result["status"] = "ok"
    return result


def read_deck_results(out):
    """{deck: True|False} from <out>/decks_results.json, written by the tests step: each deck's
    own check. Anything absent, unreadable or not a boolean is None (unknown)."""
    try:
        results = json.loads((Path(out) / "decks_results.json").read_text())
    except (OSError, ValueError):
        return {}
    if not isinstance(results, dict):
        return {}
    return {k: v for k, v in results.items() if isinstance(v, bool)}


def read_rc(out):
    p = Path(out) / "decks_rc"
    return int(p.read_text().strip()) if p.exists() else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--job", required=True)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default=None)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = decks.load(a.config)
    passed = read_deck_results(out)
    summary = {
        "job": a.job,
        "decks_exit_code": read_rc(out),
        "decks": [
            run_deck(d, a.root, out) | {"deck_passed": passed.get(d["name"])}
            for d in decks.job(cfg, a.job)
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
