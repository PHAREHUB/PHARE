#
# runs every permutation of a config's kwargs, a directory per permutation
#  env PHARE_BENCH_OUT: output directory, default .phare_bench/<datetime>
#  env PHARE_BENCH_TAR=0: keep run directories, don't tar them
#  env PHARE_BENCH_MON: resource monitoring interval seconds, 0 is off
#
import os
import shutil
import datetime
import itertools
import subprocess
from pathlib import Path

from pyphare import cpp
import pyphare.pharein as ph
from pyphare.simulator.monitoring import MonitoringOptions
from pyphare.simulator.simulator import Simulator


def datetime_now():
    return datetime.datetime.now().replace(microsecond=0).strftime("%Y%m%d-%H%M%S")


MAKE_TAR_FILE = os.getenv("PHARE_BENCH_TAR", "1") != "0"
log_dir = Path(".log")
local_dir = Path(".phare")
out_dir = Path(os.getenv("PHARE_BENCH_OUT", f".phare_bench/{datetime_now()}"))
monitoring_options = MonitoringOptions(
    interval=int(os.getenv("PHARE_BENCH_MON", "5")), rank_modulo=1
)


def sort_summaries(path):
    import re

    with open(path) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    def key(line):
        tb = float(re.search(r"'tag_buffer': (\d+)", line).group(1))
        tt = float(re.search(r"'tagging_threshold': ([0-9.]+)", line).group(1))
        return (tb, tt)

    lines.sort(key=key)

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def mkdir(dir):
    dir.mkdir(parents=True, exist_ok=True)


def clean_dir(dir):
    if dir.exists():
        shutil.rmtree(str(dir))


def make_tarfile(outdir, outfile=None):
    if MAKE_TAR_FILE and cpp.mpi_rank() == 0:
        path = outdir.parent
        outfile = outfile or outdir.name
        args = ["tar", "czf", f"{path/outfile}.tar.gz", "-C", str(path), outdir.name]
        ec = subprocess.call(args)
        assert ec == 0, f"{ec}"
        clean_dir(outdir)


def post_sim(summary, kwargs):
    log = local_dir / "summary.txt"
    with open(log, "w") as f:
        f.write(f"{kwargs}\n")
        f.write(summary)
    to = out_dir / run_name(kwargs)
    shutil.copytree(log_dir, to / log_dir.name)
    shutil.copytree(local_dir, to / local_dir.name, ignore=_other_ranks_timers)
    make_tarfile(to)


def run_name(kwargs):
    """k=v,... of the kwargs, "repeat" a subdirectory so repeats of a config
    are the runs of one directory, see phlop.timing.scope_compare"""

    def value(v):
        return "x".join(map(str, v)) if isinstance(v, (list, tuple)) else str(v)

    config = {k: v for k, v in kwargs.items() if k != "repeat"}
    name = Path(",".join(f"{k}={value(v)}" for k, v in config.items()) or "default")
    return name / f"repeat={kwargs['repeat']}" if "repeat" in kwargs else name


def execute_permutation(config_fn, keys, permutation):
    from pyphare import cpp

    def run():
        ph.global_vars.sim = None
        kwargs = {keys[i]: permutation[i] for i in range(len(keys))}
        sim = Simulator(config_fn(**kwargs)).run(monitoring=monitoring_options)
        summary = sim.summary()
        sim.reset()
        cpp.mpi_barrier()
        if cpp.mpi_rank() == 0:
            try:
                post_sim(summary, kwargs)
            except Exception:  # the other ranks wait in the next barrier
                import traceback
                from mpi4py import MPI

                traceback.print_exc()
                MPI.COMM_WORLD.Abort(1)

    run()
    cpp.mpi_barrier()


def execute(config_fn, permutables):
    from pyphare import cpp

    if cpp.mpi_rank() == 0:
        mkdir(out_dir)  # should be unique
        mkdir(log_dir)
        mkdir(local_dir / "timings")
    cpp.mpi_barrier()
    keys = tuple(e[0] for e in permutables)
    for permutation in itertools.product(*(e[1] for e in permutables)):
        execute_permutation(config_fn, keys, permutation)
    make_tarfile(out_dir)


def _other_ranks_timers(dir, names):
    """timer files of ranks not in this run, left from earlier runs"""
    if Path(dir).resolve() != (local_dir / "timings").resolve():
        return []
    ours = {f"rank.{rank}.txt" for rank in range(cpp.mpi_size())}
    return [n for n in names if n.startswith("rank.") and n not in ours]

