#
# per launch max-over-ranks inclusive time (s) of the given scopes, from a permutor bench tarball
#   PYTHONPATH=subprojects/phlop python3 tools/bench/functional/scope_times.py \
#     .phare_bench/<date>.tar.gz "Simulator::advance" "SolverPPC::moveIons" ...
#

import sys
import ast
import glob
import tarfile
import tempfile
import collections
from pathlib import Path

from phlop.timing.scope_timer import file_parser


def extract(tar_path, into):
    with tarfile.open(tar_path) as tar:
        tar.extractall(into, filter="data")
    for nested in Path(into).glob("*/*.tar.gz"):
        with tarfile.open(nested) as tar:
            tar.extractall(nested.parent, filter="data")


def scope_totals(timer_file):
    stf = file_parser(timer_file)
    totals = collections.Counter()

    def walk(node):
        totals[stf(node.k)] += node.t
        for child in node.c:
            walk(child)

    for root in stf.roots:
        walk(root)
    return totals


def rows_for(base, scopes):
    for run in sorted(glob.glob(f"{base}/*/*/")):
        summary = Path(run) / ".phare" / "summary.txt"
        if not summary.exists():
            continue
        kwargs = ast.literal_eval(summary.read_text().splitlines()[0])
        per_rank = [scope_totals(f) for f in glob.glob(f"{run}.phare/timings/rank.*.txt")]
        if not per_rank:
            continue
        launch = "x".join(str(kwargs.get(k, "?")) for k in ("mpirun", "pools", "threads"))
        label = f"{kwargs['particle_layout']:9s} {launch}"
        yield label, [max(t[s] for t in per_rank) / 1e9 for s in scopes]


def main():
    tar_path, scopes = sys.argv[1], sys.argv[2:]
    with tempfile.TemporaryDirectory() as tmp:
        extract(tar_path, tmp)
        rows = sorted(rows_for(tmp, scopes))
    print(f"{'':16s}" + "".join(f"{s.split('::')[-1][:22]:>24s}" for s in scopes))
    for label, vals in rows:
        print(f"{label:16s}" + "".join(f"{v:24.2f}" for v in vals))


if __name__ == "__main__":
    main()
