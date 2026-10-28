#
#
#
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


MAKE_TAR_FILE = True
log_dir = Path(".log")
local_dir = Path(".phare")
out_dir = Path(f".phare_bench/{datetime_now()}")
monitoring_options = MonitoringOptions(interval=5, rank_modulo=1)


def sort_summaries(path):
    import re

    with open(path) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    def key(line):
        tb = float(re.search(r"'tag_buffer': (\d+)", line).group(1))
        tt = float(re.search(r"'tagging_threshold': ([0-9.]+)", line).group(1))
        ts = float(re.search(r"'tile_size': ([0-9.]+)", line).group(1))
        return (tb, tt, ts)

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
    to = out_dir / datetime_now()
    shutil.copytree(log_dir, to / log_dir.name)
    shutil.copytree(local_dir, to / local_dir.name)
    make_tarfile(to)


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
            post_sim(summary, kwargs)

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
