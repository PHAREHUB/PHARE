#
#
#
import shutil
import tempfile
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
summary_every = 1  # coarse steps between Simulator.summary() samples


def sort_summaries(path, keys=("tag_buffer", "tagging_threshold")):
    """sorts the JSON lines of bench.py print_summary by kwargs, missing keys last"""
    import json

    with open(path) as f:
        lines = [line.rstrip("\n") for line in f if line.strip()]

    def key(line):
        kwargs = json.loads(line)["kwargs"]
        return tuple((k not in kwargs, kwargs.get(k)) for k in keys)

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


def post_sim(samples, kwargs):
    import json

    log = local_dir / "summary.json"
    with open(log, "w") as f:
        json.dump(dict(kwargs=kwargs, samples=samples), f, indent=1)
    # unique even for permutations finishing in the same second
    to = Path(tempfile.mkdtemp(prefix=f"{datetime_now()}_", dir=out_dir))
    shutil.copytree(log_dir, to / log_dir.name)
    shutil.copytree(local_dir, to / local_dir.name)
    make_tarfile(to)


def execute_permutation(config_fn, keys, permutation):
    from pyphare import cpp

    def run():
        ph.global_vars.sim = None
        kwargs = {keys[i]: permutation[i] for i in range(len(keys))}
        samples, steps = [], [0]  # summaries over time, rank 0 only

        def sample():  # collective, every rank calls it
            if summary := sim.summary():  # None on ranks > 0
                samples.append(dict(time=sim.currentTime(), **summary))

        def post_advance(time):  # after every coarse step, on every rank
            steps[0] += 1
            # the final hierarchy is sampled after run(), even if there was no step
            if time < sim.cpp_sim.endTime() and steps[0] % summary_every == 0:
                sample()

        sim = Simulator(config_fn(**kwargs), post_advance=post_advance)
        sim.run(monitoring=monitoring_options)
        sample()
        sim.reset()
        cpp.mpi_barrier()
        if cpp.mpi_rank() == 0:
            try:
                post_sim(samples, kwargs)
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
